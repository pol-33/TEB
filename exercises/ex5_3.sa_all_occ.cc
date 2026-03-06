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

  if (k == m) return 0;      // P fully matched, P is a prefix of the suffix
  if (i + k == n) return 1;  // suffix ended, suffix < P, so P > suffix
  if (P[k] < T[i + k]) return -1;
  return 1;
}

// Returns true if suffix < pattern in lexicographic order.
bool suffix_less_than_pattern(const string& T, int suf_pos, const string& P) {
  // suffix < P  <=>  P > suffix  <=> compare_pattern_to_suffix returns +1
  return compare_pattern_to_suffix(T, suf_pos, P) > 0;
}

// Small helper for DNA strings: return the next lexicographic string after P,
// assuming alphabet order A < C < G < T.
// If P is all 'T's, there is no next string; we return empty and handle it.
string next_dna_string(const string& P) {
  string s = P;
  for (int i = (int)s.size() - 1; i >= 0; --i) {
    if (s[i] == 'A') { s[i] = 'C'; return s; }
    if (s[i] == 'C') { s[i] = 'G'; return s; }
    if (s[i] == 'G') { s[i] = 'T'; return s; }
    // s[i] == 'T': carry
    s[i] = 'A';
  }
  return ""; // overflow (P was all 'T')
}

// Find [L, R) in SA such that P is a prefix of T[SA[j]..] for all j in [L, R).
// comparisons counts how many midpoint comparisons were evaluated across both searches.
pair<int,int> find_occurrence_interval(const string& T, const vector<int>& sa,
                                       const string& P, int& comparisons) {
  comparisons = 0;

  // Lower bound for P: first suffix >= P (in suffix-vs-string order)
  int lo = 0, hi = (int)sa.size();
  while (lo < hi) {
    int mid = lo + (hi - lo) / 2;
    ++comparisons;
    // if suffix < P, go right; else go left
    if (suffix_less_than_pattern(T, sa[mid], P)) lo = mid + 1;
    else hi = mid;
  }
  int L = lo;

  // Upper bound: first suffix >= next(P). For DNA, occurrences are exactly [L, ub(next(P))).
  string Pnext = next_dna_string(P);
  int R = (int)sa.size();

  if (Pnext != "") {
    lo = 0; hi = (int)sa.size();
    while (lo < hi) {
      int mid = lo + (hi - lo) / 2;
      ++comparisons;
      if (suffix_less_than_pattern(T, sa[mid], Pnext)) lo = mid + 1;
      else hi = mid;
    }
    R = lo;
  }

  // If L is out of range or SA[L] does not have P as prefix, then interval is empty.
  if (L >= (int)sa.size()) return {L, L};
  if (compare_pattern_to_suffix(T, sa[L], P) != 0) return {L, L};

  return {L, R};
}

int main() {
  string T = "ACGTACGTGACG";
  vector<int> sa = build_suffix_array(T);

  vector<string> patterns = {"ACG", "CGT", "TGA", "GATTACA", "ACGTG", "GAC", "TTT"};

  cout << "T = " << T << "\n\n";

  for (int idx = 0; idx < (int)patterns.size(); ++idx) {
    const string& P = patterns[idx];

    int comps = 0;
    pair<int,int> interval = find_occurrence_interval(T, sa, P, comps);
    int L = interval.first;
    int R = interval.second;

    cout << "P = " << P << "\n";
    cout << "  interval [L,R) = [" << L << "," << R << "), size = " << (R - L)
         << ", SA comparisons = " << comps << "\n";

    if (R == L) {
      cout << "  occurrences: (none)\n\n";
      continue;
    }

    cout << "  occurrences (positions in T): ";
    for (int j = L; j < R; ++j) {
      if (j > L) cout << ' ';
      cout << sa[j];
    }
    cout << "\n\n";
  }

  return 0;
}