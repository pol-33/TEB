#include <algorithm>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include "common.hpp"
#include "index.hpp"
#include "nucleotide.hpp"
#include "reference.hpp"

namespace {

struct Config {
    std::string reference_path;
    std::string index_path;
};

struct UniqueCount {
    uint32_t key = 0;
    uint32_t count = 0;
};

constexpr std::size_t kDensePageThreshold =
    (mapper_speed::kOffsetPageSize * sizeof(uint32_t)) / sizeof(mapper_speed::OffsetTransition);

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
        for (uint32_t i = 0; i < chrom.length; ++i) {
            const uint32_t global = chrom.start + i;
            const char base = mapper_speed::bitset_get(reference.n_mask_words.data(), global)
                ? 'N'
                : mapper_speed::kCodeToBase[mapper_speed::packed_get(reference.packed_bases.data(), global)];
            const uint8_t code = mapper_speed::base_to_code(base);
            if (code == mapper_speed::kBaseCodeInvalid) {
                valid_run = 0;
                continue;
            }
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

void append_sparse_records(const std::vector<mapper_speed::OffsetTransition>& transitions,
                           std::vector<uint8_t>& page_data,
                           mapper_speed::OffsetPageMeta& meta) {
    meta.flags = mapper_speed::kOffsetPageSparse;
    meta.record_count = static_cast<uint16_t>(transitions.size());
    meta.data_offset = page_data.size();
    const std::size_t bytes = transitions.size() * sizeof(mapper_speed::OffsetTransition);
    const std::size_t old_size = page_data.size();
    page_data.resize(old_size + bytes);
    std::memcpy(page_data.data() + old_size, transitions.data(), bytes);
}

void append_dense_page(uint32_t base_value,
                       const std::vector<mapper_speed::OffsetTransition>& transitions,
                       uint32_t value_count,
                       std::vector<uint8_t>& page_data,
                       mapper_speed::OffsetPageMeta& meta) {
    meta.flags = mapper_speed::kOffsetPageDense;
    meta.record_count = static_cast<uint16_t>(value_count);
    meta.data_offset = page_data.size();

    std::vector<uint32_t> values(value_count, base_value);
    std::size_t transition_idx = 0;
    uint32_t current = base_value;
    for (uint32_t local = 0; local < value_count; ++local) {
        while (transition_idx < transitions.size() && transitions[transition_idx].local_entry <= local) {
            current = transitions[transition_idx].value_after;
            ++transition_idx;
        }
        values[local] = current;
    }

    const std::size_t bytes = values.size() * sizeof(uint32_t);
    const std::size_t old_size = page_data.size();
    page_data.resize(old_size + bytes);
    std::memcpy(page_data.data() + old_size, values.data(), bytes);
}

void build_compact_index(const mapper_speed::ReferenceData& reference,
                         std::vector<mapper_speed::OffsetPageMeta>& page_meta,
                         std::vector<uint8_t>& page_data,
                         std::vector<uint32_t>& positions) {
    std::vector<uint32_t> valid_positions = collect_valid_positions(reference);
    std::vector<uint64_t> records;
    records.reserve(valid_positions.size());
    for (uint32_t pos : valid_positions) {
        records.push_back((uint64_t{seed_key_at(reference, pos)} << 32u) | pos);
    }

    std::sort(records.begin(), records.end());
    positions.resize(records.size());

    std::vector<UniqueCount> unique_counts;
    unique_counts.reserve(records.size());

    std::size_t i = 0;
    uint32_t out = 0;
    while (i < records.size()) {
        const uint32_t key = static_cast<uint32_t>(records[i] >> 32u);
        uint32_t count = 0;
        while (i < records.size() && static_cast<uint32_t>(records[i] >> 32u) == key) {
            positions[out++] = static_cast<uint32_t>(records[i]);
            ++count;
            ++i;
        }
        unique_counts.push_back(UniqueCount{key, count});
    }

    const uint64_t page_count =
        (mapper_speed::kOffsetEntryCount + mapper_speed::kOffsetPageSize - 1u) / mapper_speed::kOffsetPageSize;
    page_meta.assign(static_cast<std::size_t>(page_count), mapper_speed::OffsetPageMeta{});

    std::size_t unique_idx = 0;
    uint32_t running = 0;
    for (uint64_t page = 0; page < page_count; ++page) {
        const uint64_t start_key = page * mapper_speed::kOffsetPageSize;
        const uint64_t remaining = mapper_speed::kOffsetEntryCount - start_key;
        const uint32_t value_count =
            static_cast<uint32_t>(std::min<uint64_t>(mapper_speed::kOffsetPageSize, remaining));

        mapper_speed::OffsetPageMeta meta{};
        meta.base_value = running;

        const uint64_t page_key_end = std::min<uint64_t>(uint64_t{1} << (2 * mapper_speed::kSeedLength),
                                                         start_key + value_count);
        std::size_t page_unique_end = unique_idx;
        while (page_unique_end < unique_counts.size() && unique_counts[page_unique_end].key < page_key_end) {
            ++page_unique_end;
        }

        if (unique_idx != page_unique_end) {
            std::vector<mapper_speed::OffsetTransition> transitions;
            transitions.reserve(page_unique_end - unique_idx);

            uint32_t page_running = running;
            for (std::size_t idx = unique_idx; idx < page_unique_end; ++idx) {
                const uint32_t local_key = unique_counts[idx].key - static_cast<uint32_t>(start_key);
                page_running += unique_counts[idx].count;
                if (local_key + 1u < value_count) {
                    transitions.push_back(mapper_speed::OffsetTransition{
                        static_cast<uint16_t>(local_key + 1u),
                        page_running
                    });
                }
            }

            if (!transitions.empty()) {
                if (transitions.size() >= kDensePageThreshold) {
                    append_dense_page(running, transitions, value_count, page_data, meta);
                } else {
                    append_sparse_records(transitions, page_data, meta);
                }
            }

            running = page_running;
            unique_idx = page_unique_end;
        }

        page_meta[static_cast<std::size_t>(page)] = meta;
    }
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const Config cfg = parse_args(argc, argv);
        mapper_speed::ReferenceData reference = mapper_speed::load_reference(cfg.reference_path);
        std::vector<mapper_speed::OffsetPageMeta> page_meta;
        std::vector<uint8_t> page_data;
        std::vector<uint32_t> positions;

        std::cerr << "[indexer] genome length: " << reference.genome_length << " bp across "
                  << reference.chromosomes.size() << " chromosomes\n";
        std::cerr << "[indexer] build mode: compact-page-transitions\n";

        build_compact_index(reference, page_meta, page_data, positions);

        const std::size_t bytes =
            mapper_speed::write_index(cfg.index_path, reference, page_meta, page_data, positions);
        std::cerr << "[indexer] valid " << mapper_speed::kSeedLength << "-mers: " << positions.size() << "\n";
        std::cerr << "[indexer] index written: " << mapper_speed::bytes_to_mebibytes(bytes) << " MiB\n";
        std::cerr << "[indexer] peak RSS: " << mapper_speed::bytes_to_mebibytes(mapper_speed::peak_rss_bytes()) << " MiB\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "indexer error: " << e.what() << '\n';
        return 1;
    }
}
