#include <algorithm>
#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <limits>
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
    std::string mode = "auto";
};

constexpr std::size_t kDensePageThreshold =
    (mapper_speed::kOffsetPageSize * sizeof(uint32_t)) / sizeof(mapper_speed::OffsetTransition);

void usage() {
    std::cerr << "Usage: indexer -R genome.fa -I genome.idx [--mode auto|dense|compact]\n";
}

Config parse_args(int argc, char* argv[]) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-R" && i + 1 < argc) {
            cfg.reference_path = argv[++i];
        } else if (arg == "-I" && i + 1 < argc) {
            cfg.index_path = argv[++i];
        } else if (arg == "--mode" && i + 1 < argc) {
            cfg.mode = argv[++i];
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
    if (cfg.mode != "auto" && cfg.mode != "dense" && cfg.mode != "compact") {
        throw std::runtime_error("invalid --mode, expected auto|dense|compact");
    }
    return cfg;
}

template <typename Callback>
void for_each_valid_seed(const mapper_speed::ReferenceData& reference,
                         uint32_t index_stride,
                         Callback&& callback) {
    constexpr uint32_t kRollingMask = 0xFFFFFFFFu;

    for (const auto& chrom : reference.chromosomes) {
        uint32_t rolling_key = 0;
        uint32_t valid_run = 0;
        for (uint32_t i = 0; i < chrom.length; ++i) {
            const uint32_t global = chrom.start + i;
            const uint8_t code = mapper_speed::bitset_get(reference.n_mask_words.data(), global)
                ? mapper_speed::kBaseCodeInvalid
                : mapper_speed::packed_get(reference.packed_bases.data(), global);
            if (code == mapper_speed::kBaseCodeInvalid) {
                rolling_key = 0;
                valid_run = 0;
                continue;
            }

            rolling_key = ((rolling_key << 2u) | code) & kRollingMask;
            if (valid_run < mapper_speed::kSeedLength) {
                ++valid_run;
            }
            if (valid_run >= mapper_speed::kSeedLength) {
                const uint32_t start = global - mapper_speed::kSeedLength + 1u;
                if ((start % index_stride) == 0u) {
                    callback(rolling_key, start);
                }
            }
        }
    }
}

std::size_t count_valid_positions(const mapper_speed::ReferenceData& reference, uint32_t index_stride) {
    std::size_t count = 0;
    for_each_valid_seed(reference, index_stride, [&](uint32_t, uint32_t) {
        ++count;
    });
    return count;
}

uint64_t estimate_dense_index_bytes(const mapper_speed::ReferenceData& reference,
                                    std::size_t dense_positions) {
    const uint64_t page_count =
        (mapper_speed::kOffsetEntryCount + mapper_speed::kOffsetPageSize - 1u) / mapper_speed::kOffsetPageSize;
    const uint64_t page_meta_bytes = page_count * sizeof(mapper_speed::OffsetPageMeta);
    const uint64_t chromosome_bytes = static_cast<uint64_t>(reference.chromosomes.size()) *
        sizeof(mapper_speed::StoredChromosome);
    const uint64_t positions_bytes = static_cast<uint64_t>(dense_positions) * sizeof(uint32_t);
    const uint64_t reference_bytes = reference.packed_bases.size() +
        reference.n_mask_words.size() * sizeof(uint64_t);
    return sizeof(mapper_speed::IndexHeader) + chromosome_bytes + page_meta_bytes +
        reference_bytes + positions_bytes + (256u << 20u);
}

uint32_t choose_index_stride(const mapper_speed::ReferenceData& reference,
                             const Config& cfg) {
    if (cfg.mode == "dense") {
        return mapper_speed::kDenseIndexStride;
    }
    if (cfg.mode == "compact") {
        return mapper_speed::kCompactIndexStride;
    }

    const uint64_t ram = mapper_speed::physical_memory_bytes();
    if (ram == 0u) {
        return mapper_speed::kCompactIndexStride;
    }

    const std::size_t dense_positions = count_valid_positions(reference, mapper_speed::kDenseIndexStride);
    const uint64_t estimated_dense_bytes = estimate_dense_index_bytes(reference, dense_positions);
    return (estimated_dense_bytes <= (ram * 3u) / 5u)
        ? mapper_speed::kDenseIndexStride
        : mapper_speed::kCompactIndexStride;
}

unsigned choose_bucket_bits(std::size_t valid_positions) {
    const uint64_t ram = mapper_speed::physical_memory_bytes();
    const uint64_t target_bucket_bytes = ram == 0u
        ? (uint64_t{1} << 30u)
        : std::max<uint64_t>(uint64_t{512} << 20u, ram / 10u);

    unsigned bits = 0;
    while (bits < 12u) {
        const uint64_t bucket_count = uint64_t{1} << bits;
        const uint64_t avg_bytes =
            (static_cast<uint64_t>(valid_positions) * sizeof(uint64_t) + bucket_count - 1u) / bucket_count;
        if (avg_bytes <= target_bucket_bytes) {
            break;
        }
        ++bits;
    }
    return bits;
}

void collect_bucket_records(const mapper_speed::ReferenceData& reference,
                            uint32_t index_stride,
                            unsigned bucket_bits,
                            uint32_t bucket_id,
                            std::size_t reserve_hint,
                            std::vector<uint64_t>& records) {
    records.clear();
    if (reserve_hint > 0u) {
        records.reserve(reserve_hint);
    }
    const unsigned shift = 32u - bucket_bits;
    for_each_valid_seed(reference, index_stride, [&](uint32_t key, uint32_t pos) {
        if (bucket_bits == 0u || (key >> shift) == bucket_id) {
            records.push_back((uint64_t{key} << 32u) | pos);
        }
    });
}

void append_padding(std::ofstream& out, uint64_t target_offset) {
    const uint64_t current = static_cast<uint64_t>(out.tellp());
    if (current >= target_offset) {
        return;
    }
    std::vector<char> padding(static_cast<std::size_t>(target_offset - current), 0);
    out.write(padding.data(), static_cast<std::streamsize>(padding.size()));
    if (!out) {
        throw std::runtime_error("failed to write index padding");
    }
}

template <typename T>
void write_value(std::ofstream& out, const T& value) {
    out.write(reinterpret_cast<const char*>(&value), sizeof(T));
    if (!out) {
        throw std::runtime_error("failed to write index file");
    }
}

class DirectIndexWriter {
public:
    DirectIndexWriter(const std::string& path,
                      const mapper_speed::ReferenceData& reference,
                      uint32_t index_stride,
                      std::size_t page_count,
                      std::size_t positions_count)
        : out_(path, std::ios::binary | std::ios::trunc) {
        if (!out_) {
            throw std::runtime_error("failed to open index output: " + path);
        }

        std::memcpy(header_.magic, "MSPDIDX1", 8u);
        header_.version = mapper_speed::kIndexVersion;
        header_.flags = index_stride;
        header_.checksum = reference.checksum;
        header_.genome_length = reference.genome_length;
        header_.chromosome_count = static_cast<uint32_t>(reference.chromosomes.size());
        header_.packed_reference_bytes = reference.packed_bases.size();
        header_.n_mask_bytes = reference.n_mask_words.size() * sizeof(uint64_t);
        header_.offset_page_count = page_count;
        header_.positions_count = positions_count;

        uint64_t offset = sizeof(mapper_speed::IndexHeader);
        header_.chromosome_table_offset = offset;
        for (const auto& chrom : reference.chromosomes) {
            offset += sizeof(mapper_speed::StoredChromosome) + chrom.name.size();
        }
        offset = mapper_speed::align_up(offset, 64u);

        header_.packed_reference_offset = offset;
        offset += header_.packed_reference_bytes;
        offset = mapper_speed::align_up(offset, 64u);

        header_.n_mask_offset = offset;
        offset += header_.n_mask_bytes;
        offset = mapper_speed::align_up(offset, 64u);

        header_.offset_page_meta_offset = offset;
        offset += static_cast<uint64_t>(page_count) * sizeof(mapper_speed::OffsetPageMeta);
        offset = mapper_speed::align_up(offset, 64u);

        header_.positions_offset = offset;
        offset += static_cast<uint64_t>(positions_count) * sizeof(uint32_t);
        offset = mapper_speed::align_up(offset, 64u);

        page_data_offset_ = offset;

        write_value(out_, header_);

        for (const auto& chrom : reference.chromosomes) {
            mapper_speed::StoredChromosome stored{};
            stored.name_len = static_cast<uint32_t>(chrom.name.size());
            stored.start = chrom.start;
            stored.length = chrom.length;
            write_value(out_, stored);
            out_.write(chrom.name.data(), static_cast<std::streamsize>(chrom.name.size()));
            if (!out_) {
                throw std::runtime_error("failed to write chromosome name");
            }
        }

        append_padding(out_, header_.packed_reference_offset);
        out_.write(reinterpret_cast<const char*>(reference.packed_bases.data()),
                   static_cast<std::streamsize>(reference.packed_bases.size()));

        append_padding(out_, header_.n_mask_offset);
        out_.write(reinterpret_cast<const char*>(reference.n_mask_words.data()),
                   static_cast<std::streamsize>(reference.n_mask_words.size() * sizeof(uint64_t)));

        append_padding(out_, header_.offset_page_meta_offset);
        std::vector<char> zero_meta(static_cast<std::size_t>(page_count * sizeof(mapper_speed::OffsetPageMeta)), 0);
        out_.write(zero_meta.data(), static_cast<std::streamsize>(zero_meta.size()));

        append_padding(out_, header_.positions_offset);
        positions_cursor_ = header_.positions_offset;
        page_data_cursor_ = page_data_offset_;
    }

    void append_positions(const uint32_t* data, std::size_t count) {
        out_.seekp(static_cast<std::streamoff>(positions_cursor_));
        out_.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(count * sizeof(uint32_t)));
        if (!out_) {
            throw std::runtime_error("failed to append positions");
        }
        positions_cursor_ += count * sizeof(uint32_t);
    }

    uint64_t append_page_data(const void* data, std::size_t bytes) {
        out_.seekp(static_cast<std::streamoff>(page_data_cursor_));
        const uint64_t offset = page_data_cursor_;
        out_.write(reinterpret_cast<const char*>(data), static_cast<std::streamsize>(bytes));
        if (!out_) {
            throw std::runtime_error("failed to append page data");
        }
        page_data_cursor_ += bytes;
        return offset;
    }

    std::size_t finalize(const std::vector<mapper_speed::OffsetPageMeta>& page_meta) {
        out_.seekp(static_cast<std::streamoff>(header_.offset_page_meta_offset));
        out_.write(reinterpret_cast<const char*>(page_meta.data()),
                   static_cast<std::streamsize>(page_meta.size() * sizeof(mapper_speed::OffsetPageMeta)));
        if (!out_) {
            throw std::runtime_error("failed to write page metadata");
        }

        out_.seekp(0);
        write_value(out_, header_);
        out_.seekp(0, std::ios::end);
        out_.flush();
        if (!out_) {
            throw std::runtime_error("failed to finalize index");
        }
        return static_cast<std::size_t>(out_.tellp());
    }

private:
    std::ofstream out_;
    mapper_speed::IndexHeader header_{};
    uint64_t page_data_offset_ = 0;
    uint64_t positions_cursor_ = 0;
    uint64_t page_data_cursor_ = 0;
};

void append_sparse_records(const std::vector<mapper_speed::OffsetTransition>& transitions,
                           DirectIndexWriter& writer,
                           mapper_speed::OffsetPageMeta& meta) {
    meta.flags = mapper_speed::kOffsetPageSparse;
    meta.record_count = static_cast<uint16_t>(transitions.size());
    meta.data_offset = writer.append_page_data(transitions.data(),
                                               transitions.size() * sizeof(mapper_speed::OffsetTransition));
}

void append_dense_page(uint32_t base_value,
                       const std::vector<mapper_speed::OffsetTransition>& transitions,
                       uint32_t value_count,
                       DirectIndexWriter& writer,
                       mapper_speed::OffsetPageMeta& meta) {
    meta.flags = mapper_speed::kOffsetPageDense;
    meta.record_count = static_cast<uint16_t>(value_count);

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

    meta.data_offset = writer.append_page_data(values.data(), values.size() * sizeof(uint32_t));
}

std::size_t build_multipass_index(const mapper_speed::ReferenceData& reference,
                                  const std::string& index_path,
                                  uint32_t index_stride,
                                  std::size_t valid_positions,
                                  unsigned bucket_bits) {
    const uint64_t page_count =
        (mapper_speed::kOffsetEntryCount + mapper_speed::kOffsetPageSize - 1u) / mapper_speed::kOffsetPageSize;
    std::vector<mapper_speed::OffsetPageMeta> page_meta(
        static_cast<std::size_t>(page_count), mapper_speed::OffsetPageMeta{});

    DirectIndexWriter writer(index_path, reference, index_stride, page_meta.size(), valid_positions);
    const uint32_t bucket_count = uint32_t{1} << bucket_bits;
    const std::size_t reserve_hint = (valid_positions + bucket_count - 1u) / bucket_count;

    std::vector<uint64_t> records;
    std::vector<mapper_speed::OffsetTransition> transitions;
    std::vector<uint32_t> grouped_positions;
    std::size_t positions_written = 0;

    for (uint32_t bucket = 0; bucket < bucket_count; ++bucket) {
        std::cerr << "[indexer] bucket " << (bucket + 1u) << "/" << bucket_count << '\n';
        collect_bucket_records(reference, index_stride, bucket_bits, bucket, reserve_hint, records);
        std::sort(records.begin(), records.end());

        const uint64_t bucket_start_key = static_cast<uint64_t>(bucket) << (32u - bucket_bits);
        const uint64_t bucket_end_key = (bucket + 1u == bucket_count)
            ? (uint64_t{1} << 32u)
            : (static_cast<uint64_t>(bucket + 1u) << (32u - bucket_bits));
        const uint64_t start_page = bucket_start_key >> mapper_speed::kOffsetPageShift;
        const uint64_t end_page =
            (bucket_end_key + mapper_speed::kOffsetPageSize - 1u) >> mapper_speed::kOffsetPageShift;

        std::size_t record_idx = 0;
        for (uint64_t page = start_page; page < end_page; ++page) {
            const uint64_t start_key = page * mapper_speed::kOffsetPageSize;
            const uint64_t remaining = mapper_speed::kOffsetEntryCount - start_key;
            const uint32_t value_count =
                static_cast<uint32_t>(std::min<uint64_t>(mapper_speed::kOffsetPageSize, remaining));
            const uint64_t page_key_end = std::min<uint64_t>(uint64_t{1} << 32u, start_key + value_count);

            mapper_speed::OffsetPageMeta meta{};
            meta.base_value = static_cast<uint32_t>(positions_written);
            transitions.clear();

            while (record_idx < records.size() && static_cast<uint32_t>(records[record_idx] >> 32u) < page_key_end) {
                const uint32_t key = static_cast<uint32_t>(records[record_idx] >> 32u);
                const std::size_t group_begin = record_idx;
                while (record_idx < records.size() && static_cast<uint32_t>(records[record_idx] >> 32u) == key) {
                    ++record_idx;
                }
                const std::size_t group_size = record_idx - group_begin;
                grouped_positions.resize(group_size);
                for (std::size_t i = 0; i < group_size; ++i) {
                    grouped_positions[i] = static_cast<uint32_t>(records[group_begin + i]);
                }
                writer.append_positions(grouped_positions.data(), group_size);
                positions_written += group_size;

                const uint32_t local_key = key - static_cast<uint32_t>(start_key);
                if (local_key + 1u < value_count) {
                    transitions.push_back(mapper_speed::OffsetTransition{
                        static_cast<uint16_t>(local_key + 1u),
                        static_cast<uint32_t>(positions_written)
                    });
                }
            }

            if (!transitions.empty()) {
                if (transitions.size() >= kDensePageThreshold) {
                    append_dense_page(meta.base_value, transitions, value_count, writer, meta);
                } else {
                    append_sparse_records(transitions, writer, meta);
                }
            }
            page_meta[static_cast<std::size_t>(page)] = meta;
        }
    }

    page_meta.back().base_value = static_cast<uint32_t>(positions_written);
    return writer.finalize(page_meta);
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const Config cfg = parse_args(argc, argv);
        mapper_speed::ReferenceData reference = mapper_speed::load_reference(cfg.reference_path);
        const uint32_t index_stride = choose_index_stride(reference, cfg);

        std::cerr << "[indexer] genome length: " << reference.genome_length << " bp across "
                  << reference.chromosomes.size() << " chromosomes\n";

        const std::size_t valid_positions = count_valid_positions(reference, index_stride);
        const unsigned bucket_bits = choose_bucket_bits(valid_positions);
        const uint32_t bucket_count = uint32_t{1} << bucket_bits;
        const char* mode_label = (index_stride == mapper_speed::kDenseIndexStride) ? "dense" : "compact";

        std::cerr << "[indexer] build mode: compact-page-transitions/" << mode_label
                  << " (multipass-" << bucket_count << "-bucket, stride-" << index_stride << ")\n";

        const std::size_t bytes =
            build_multipass_index(reference, cfg.index_path, index_stride, valid_positions, bucket_bits);

        std::cerr << "[indexer] valid " << mapper_speed::kSeedLength << "-mers: " << valid_positions << "\n";
        std::cerr << "[indexer] index written: " << mapper_speed::bytes_to_mebibytes(bytes) << " MiB\n";
        std::cerr << "[indexer] peak RSS: " << mapper_speed::bytes_to_mebibytes(mapper_speed::peak_rss_bytes()) << " MiB\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "indexer error: " << e.what() << '\n';
        return 1;
    }
}
