#include <cstdint>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>

#include "buffered_io.hpp"
#include "common.hpp"
#include "index.hpp"
#include "search.hpp"

namespace {

struct Config {
    std::string index_path;
    std::string reads_path;
    std::string output_path;
    int max_errors = -1;
};

void usage() {
    std::cerr << "Usage: mapper -I genome.idx -i reads.fastq -o output.sam -k <0..3>\n";
}

Config parse_args(int argc, char* argv[]) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-I" && i + 1 < argc) {
            cfg.index_path = argv[++i];
        } else if (arg == "-i" && i + 1 < argc) {
            cfg.reads_path = argv[++i];
        } else if (arg == "-o" && i + 1 < argc) {
            cfg.output_path = argv[++i];
        } else if (arg == "-k" && i + 1 < argc) {
            cfg.max_errors = std::stoi(argv[++i]);
        } else if (arg == "-h" || arg == "--help") {
            usage();
            std::exit(0);
        } else {
            usage();
            throw std::runtime_error("unknown argument: " + arg);
        }
    }
    if (cfg.index_path.empty() || cfg.reads_path.empty() || cfg.output_path.empty() ||
        cfg.max_errors < 0 || cfg.max_errors > 3) {
        usage();
        throw std::runtime_error("missing required arguments or invalid -k");
    }
    return cfg;
}

}  // namespace

int main(int argc, char* argv[]) {
    try {
        const Config cfg = parse_args(argc, argv);
        mapper_speed::IndexView index(cfg.index_path);
        mapper_speed::FastqReader reader(cfg.reads_path);
        mapper_speed::BufferedWriter writer(cfg.output_path);
        mapper_speed::MapperEngine engine(index);

        std::cerr << "[mapper] loaded index with " << index.chromosome_count() << " chromosomes, "
                  << index.positions_count() << " indexed seeds\n";
        std::cerr << "[mapper] single-thread verifier dispatch ready\n";

        mapper_speed::FastqRecord record;
        uint64_t processed = 0;
        while (reader.next(record)) {
            writer.write(engine.map_record(record, cfg.max_errors));
            ++processed;
            if ((processed % 10000u) == 0u) {
                std::cerr << "[mapper] processed " << processed << " reads\n";
            }
        }
        writer.flush();

        std::cerr << "[mapper] processed " << processed << " reads total\n";
        std::cerr << "[mapper] peak RSS: " << mapper_speed::bytes_to_mebibytes(mapper_speed::peak_rss_bytes()) << " MiB\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "mapper error: " << e.what() << '\n';
        return 1;
    }
}
