#include "fasta_parser.hpp"

#include <time.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdexcept>
#include <iomanip>
#include <algorithm>

using namespace std;

int main(int argc, char* argv[]) {
    clock_t t = clock();
    string text, pattern;
    cin >> text;
    cin >> pattern;
    naive();

    t = clock() -t;
    printf("Total time std (micro-seconds): %ld \n", t);
    return 0;
}