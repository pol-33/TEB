#include "utils.hpp"



void print_kmer_table(const kmer_table_t& table) {
    for (const auto& pair : table) {
        const string& kmer = pair.first;
        const list<int>& positions = pair.second;

        printf("%s : ", kmer.c_str());

        for (int pos : positions)
            printf("%d ", pos);

        printf("\n");
    }
}