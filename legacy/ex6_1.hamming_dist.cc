#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// Task 1.1: Hamming distance
// Computes the number of mismatched positions between two equal-length sequences
int hamming_distance(const string& X, const string& Y, int& comps) {
  if (X.length() != Y.length()) {
    return -1; // Error: sequences must have equal length
  }
  
  int dist = 0;
  for (size_t i = 0; i < X.length(); ++i) {
    comps++;
    if (X[i] != Y[i]) {
      dist++;
    }
  }
  return dist;
}

// Task 1.2: Naive approximate search
// Finds all occurrences of pattern P in text T with at most k mismatches
// Scans every possible alignment and verifies using Hamming distance
vector<int> naive_approximate_search(const string& T, const string& P, int k, int& comps) {
  vector<int> matches;
  
  // Scan every possible alignment
  for (size_t i = 0; i + P.length() <= T.length(); ++i) {
    // Extract substring of length |P| starting at position i
    string substring = T.substr(i, P.length());
    
    // Compute Hamming distance
    int dist = hamming_distance(substring, P, comps);
    
    // If distance <= k, it's a valid approximate match
    if (dist <= k) {
      matches.push_back(i);
    }
  }
  
  return matches;
}

int main() {
  string T = "ACGTACGTGACG";

  vector<string> patterns = {"ACG", "CGT", "TGA", "GATTACA", "ACGTG", "GAC", "TTT"};

  cout << "T = " << T << "\n\n";

  for (int idx = 0; idx < (int)patterns.size(); ++idx) {
    const string& P = patterns[idx];

    int comps = 0;
    // Search with k=0 (exact matches only)
    vector<int> matches = naive_approximate_search(T, P, 1, comps);

    cout << "P = " << P << "\n";
    cout << "  occurrences (positions in T): ";
    
    if (matches.empty()) {
      cout << "(none)";
    } else {
      for (size_t j = 0; j < matches.size(); ++j) {
        if (j > 0) cout << ' ';
        cout << matches[j];
      }
    }
    
    cout << "\n";
    cout << "  character comparisons = " << comps << "\n\n";
  }

  return 0;
}