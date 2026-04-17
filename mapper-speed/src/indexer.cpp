#include <algorithm>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <sys/mman.h>

#include "common.hpp"
#include "index.hpp"
#include "nucleotide.hpp"
#include "reference.hpp"

namespace {

struct Config {
    std::string reference_path;
    std::string index_path;
};

void usage() {
    std::cerr << "Usage: indexer -R genome.fa -I genome.idx\n";
}

Config parse_args(int argc, char* argv[]) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-R" && i + 1 < argc) {
            cfg.reference_path = argv[++i];
        } else if (arg == "-I" && i + 1 < argc) {
            cfg.index_path = argv[++i];
        } else if (arg == "-h" || arg == "--help") {
            usage();
            std::exit(0);
        } else {
            usage();
            throw std::runtime_error("unknown argument: " + arg);
        }
    }
    if (cfg.reference_path.empty() || cfg.index_path.empty()) {
        usage();
        throw std::runtime_error("missing required arguments");
    }
    return cfg;
}

struct DenseCounts {
    uint32_t* values = nullptr;
    std::size_t bytes = 0;

    explicit DenseCounts(std::size_t count) : bytes(count * sizeof(uint32_t)) {
        values = static_cast<uint32_t*>(mmap(nullptr, bytes, PROT_READ | PROT_WRITE,
                                             MAP_PRIVATE | MAP_ANON, -1, 0));
        if (values == MAP_FAILED) {
            values = nullptr;
            throw std::runtime_error("failed to allocate dense count table");
        }
        std::memset(values, 0, bytes);
    }

    ~DenseCounts() {
        if (values != nullptr) {
            munmap(values, bytes);
        }
    }
};

std::vector<uint32_t> collect_valid_positions(const mapper_speed::ReferenceData& reference) {
    std::vector<uint32_t> positions;
    if (reference.genome_length < mapper_speed::kSeedLength) {
        return positions;
    }
    positions.reserve(reference.genome_length);

    for (const auto& chrom : reference.chromosomes) {
        if (chrom.length < mapper_speed::kSeedLength) {
            continue;
        }
        uint32_t valid_run = 0;
        uint32_t key = 0;
        for (uint32_t i = 0; i < chrom.length; ++i) {
            const uint32_t global = chrom.start + i;
            const char base = mapper_speed::bitset_get(reference.n_mask_words.data(), global)
                ? 'N'
                : mapper_speed::kCodeToBase[mapper_speed::packed_get(reference.packed_bases.data(), global)];
            const uint8_t code = mapper_speed::base_to_code(base);
            if (code == mapper_speed::kBaseCodeInvalid) {
                valid_run = 0;
                key = 0;
                continue;
            }
            key = (key << 2u) | code;
            if (valid_run < mapper_speed::kSeedLength) {
                ++valid_run;
            }
            if (valid_run >= mapper_speed::kSeedLength) {
                positions.push_back(global - mapper_speed::kSeedLength + 1u);
            }
        }
    }
    return positions;
}

uint32_t seed_key_at(const mapper_speed::ReferenceData& reference, uint32_t global_pos) {
    uint32_t key = 0;
    for (uint32_t i = 0; i < mapper_speed::kSeedLength; ++i) {
        const uint32_t pos = global_pos + i;
        const uint8_t code = mapper_speed::bitset_get(reference.n_mask_words.data(), pos)
            ? mapper_speed::kBaseCodeInvalid
            : mapper_speed::packed_get(reference.packed_bases.data(), pos);
        if (code == mapper_speed::kBaseCodeInvalid) {
            return 0;
        }
        key = (key << 2u) | code;
    }
    return key;
}

void build_sparse(const mapper_speed::ReferenceData& reference,
                  std::vector<mapper_speed::OffsetPageMeta>& page_meta,
                  std::vector<uint8_t>& dense_pages,
                  std::vector<uint32_t>& positions) {
    std::vector<uint32_t> valid_positions = collect_valid_positions(reference);
    std::vector<uint64_t> records;
    records.reserve(valid_positions.size());
    for (uint32_t pos : valid_positions) {
        const uint32_t key = seed_key_at(reference, pos);
        records.push_back((uint64_t{key} << 32u) | pos);
    }

    std::sort(records.begin(), records.end());
    positions.resize(records.size());

    struct UniqueCount {
        uint32_t key = 0;
        uint32_t count = 0;
    };
    std::vector<UniqueCount> unique_counts;
    unique_counts.reserve(records.size());

    uint32_t out = 0;
    std::size_t i = 0;
    while (i < records.size()) {
        const uint32_t key = static_cast<uint32_t>(records[i] >> 32u);
        uint32_t count = 0;
        while (i < records.size() && static_cast<uint32_t>(records[i] >> 32u) == key) {
            positions[out++] = static_cast<uint32_t>(records[i]);
            ++count;
            ++i;
        }
        unique_counts.push_back({key, count});
    }

    const uint64_t page_count = (mapper_speed::kOffsetEntryCount + mapper_speed::kOffsetPageSize - 1u) / mapper_speed::kOffsetPageSize;
    page_meta.assign(static_cast<std::size_t>(page_count), mapper_speed::OffsetPageMeta{});

    std::size_t unique_idx = 0;
    uint32_t running = 0;
    for (uint64_t page = 0; page < page_count; ++page) {
        const uint64_t start_key = page * mapper_speed::kOffsetPageSize;
        const uint64_t remaining = mapper_speed::kOffsetEntryCount - start_key;
        const uint32_t value_count = static_cast<uint32_t>(std::min<uint64_t>(mapper_speed::kOffsetPageSize, remaining));
        mapper_speed::OffsetPageMeta meta{};
        meta.value_count = value_count;
        meta.base_value = running;

        const uint64_t page_key_end = std::min<uint64_t>(uint64_t{1} << (2 * mapper_speed::kSeedLength), start_key + value_count);
        std::size_t page_unique_end = unique_idx;
        while (page_unique_end < unique_counts.size() && unique_counts[page_unique_end].key < page_key_end) {
            ++page_unique_end;
        }

        const bool needs_materialized_page =
            (unique_idx != page_unique_end) || (start_key + value_count == mapper_speed::kOffsetEntryCount);
        if (needs_materialized_page) {
            std::vector<uint32_t> values(value_count, running);
            std::size_t page_pos = unique_idx;
            uint32_t page_running = running;
            for (uint32_t local = 0; local < value_count; ++local) {
                const uint64_t entry = start_key + local;
                values[local] = page_running;
                if (entry >= (uint64_t{1} << (2 * mapper_speed::kSeedLength))) {
                    continue;
                }
                if (page_pos < page_unique_end && unique_counts[page_pos].key == entry) {
                    page_running += unique_counts[page_pos].count;
                    ++page_pos;
                }
            }
            running = page_running;
            unique_idx = page_unique_end;

            const bool constant = std::all_of(values.begin(), values.end(), [&](uint32_t v) { return v == values.front(); });
            meta.base_value = values.front();
            if (!constant) {
                meta.flags = mapper_speed::kOffsetPageDense;
                meta.data_offset = dense_pages.size();
                const std::size_t bytes = values.size() * sizeof(uint32_t);
                const std::size_t old = dense_pages.size();
                dense_pages.resize(old + bytes);
                std::memcpy(dense_pages.data() + old, values.data(), bytes);
            }
        }
        page_meta[static_cast<std::size_t>(page)] = meta;
    }
}

void build_dense(const mapper_speed::ReferenceData& reference,
                 std::vector<mapper_speed::OffsetPageMeta>& page_meta,
                 std::vector<uint8_t>& dense_pages,
                 std::vector<uint32_t>& positions) {
    DenseCounts counts(static_cast<std::size_t>(uint64_t{1} << (2 * mapper_speed::kSeedLength)));
    std::vector<uint32_t> valid_positions = collect_valid_positions(reference);
    for (uint32_t pos : valid_positions) {
        ++counts.values[seed_key_at(reference, pos)];
    }

    uint32_t running = 0;
    for (uint64_t key = 0; key < (uint64_t{1} << (2 * mapper_speed::kSeedLength)); ++key) {
        const uint32_t count = counts.values[key];
        counts.values[key] = running;
        running += count;
    }

    positions.assign(valid_positions.size(), 0);
    const uint64_t page_count = (mapper_speed::kOffsetEntryCount + mapper_speed::kOffsetPageSize - 1u) / mapper_speed::kOffsetPageSize;
    page_meta.assign(static_cast<std::size_t>(page_count), mapper_speed::OffsetPageMeta{});

    for (uint64_t page = 0; page < page_count; ++page) {
        const uint64_t start_key = page * mapper_speed::kOffsetPageSize;
        const uint64_t remaining = mapper_speed::kOffsetEntryCount - start_key;
        const uint32_t value_count = static_cast<uint32_t>(std::min<uint64_t>(mapper_speed::kOffsetPageSize, remaining));
        std::vector<uint32_t> values(value_count, 0);
        for (uint32_t local = 0; local < value_count; ++local) {
            const uint64_t entry = start_key + local;
            values[local] = (entry == (uint64_t{1} << (2 * mapper_speed::kSeedLength))) ? running : counts.values[entry];
        }

        const bool constant = std::all_of(values.begin(), values.end(), [&](uint32_t v) { return v == values.front(); });
        mapper_speed::OffsetPageMeta meta{};
        meta.base_value = values.front();
        meta.value_count = value_count;
        if (!constant) {
            meta.flags = mapper_speed::kOffsetPageDense;
            meta.data_offset = dense_pages.size();
            const std::size_t bytes = values.size() * sizeof(uint32_t);
            const std::size_t old = dense_pages.size();
            dense_pages.resize(old + bytes);
            std::memcpy(dense_pages.data() + old, values.data(), bytes);
        }
        page_meta[static_cast<std::size_t>(page)] = meta;
    }

    for (uint32_t pos : valid_positions) {
        const uint32_t key = seed_key_at(reference, pos);
        positions[counts.values[key]++] = pos;
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const Config cfg = parse_args(argc, argv);
        mapper_speed::ReferenceData reference = mapper_speed::load_reference(cfg.reference_path);
        std::vector<mapper_speed::OffsetPageMeta> page_meta;
        std::vector<uint8_t> dense_pages;
        std::vector<uint32_t> positions;

        const uint64_t phys_mem = mapper_speed::physical_memory_bytes();
        const bool use_dense = phys_mem >= (uint64_t{40} << 30u);

        std::cerr << "[indexer] genome length: " << reference.genome_length << " bp across "
                  << reference.chromosomes.size() << " chromosomes\n";
        std::cerr << "[indexer] build mode: " << (use_dense ? "dense" : "sparse") << "\n";

        if (use_dense) {
            build_dense(reference, page_meta, dense_pages, positions);
        } else {
            build_sparse(reference, page_meta, dense_pages, positions);
        }

        const std::size_t bytes = mapper_speed::write_index(cfg.index_path, reference, page_meta, dense_pages, positions);
        std::cerr << "[indexer] valid " << mapper_speed::kSeedLength << "-mers: " << positions.size() << "\n";
        std::cerr << "[indexer] index written: " << mapper_speed::bytes_to_mebibytes(bytes) << " MiB\n";
        std::cerr << "[indexer] peak RSS: " << mapper_speed::bytes_to_mebibytes(mapper_speed::peak_rss_bytes()) << " MiB\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "indexer error: " << e.what() << '\n';
        return 1;
    }
}
