#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

// Compare suffixes T[i..] and T[j..] lexicographically
bool suffix_less(const string& T, int i, int j) {
  int n = (int)T.size();
  while (i < n && j < n) {
    if (T[i] < T[j]) return true;
    if (T[i] > T[j]) return false;
    ++i;
    ++j;
  }
  return i == n && j != n;
}

// Naive suffix array construction
vector<int> build_suffix_array(const string& T) {
  int n = (int)T.size();
  vector<int> sa(n);
  for (int i = 0; i < n; ++i) sa[i] = i;

  sort(sa.begin(), sa.end(), [&](int a, int b) {
    return suffix_less(T, a, b);
  });

  return sa;
}

// Compare pattern P with suffix T[i..] lexicographically.
// Return -1 if P < suffix, 0 if P is a prefix of suffix, +1 if P > suffix.
int compare_pattern_to_suffix(const string& T, int i, const string& P) {
  int n = (int)T.size();
  int m = (int)P.size();

  int k = 0;
  while (k < m && i + k < n && P[k] == T[i + k]) ++k;

  if (k == m) return 0;                 // P fully matched, P is a prefix
  if (i + k == n) return 1;             // suffix ended, suffix < P, so P > suffix
  if (P[k] < T[i + k]) return -1;
  return 1;
}

// Binary search over suffix array to test if P occurs in T.
// comparisons counts how many SA midpoints were evaluated.
bool contains_pattern_sa(const string& T, const vector<int>& sa,
                         const string& P, int& comparisons) {
  comparisons = 0;
  int lo = 0;
  int hi = (int)sa.size() - 1;

  while (lo <= hi) {
    int mid = lo + (hi - lo) / 2;
    ++comparisons;

    int c = compare_pattern_to_suffix(T, sa[mid], P);
    if (c == 0) return true;
    if (c < 0) hi = mid - 1;
    else lo = mid + 1;
  }
  return false;
}

int main() {
  string T = "ACGTACGTGACG";
  vector<int> sa = build_suffix_array(T);

  vector<string> patterns = {"ACG", "CGT", "TGA", "GATTACA", "ACGTG", "GAC"};

  cout << "T = " << T << "\n\n";
  cout << "Testing patterns using binary search on suffix array:\n";

  for (int idx = 0; idx < (int)patterns.size(); ++idx) {
    const string& P = patterns[idx];
    int comps = 0;
    bool ok = contains_pattern_sa(T, sa, P, comps);
    cout << "P = " << P << "\toccurs = " << (ok ? "yes" : "no")
         << "\tSA comparisons = " << comps << "\n";
  }

  return 0;
}