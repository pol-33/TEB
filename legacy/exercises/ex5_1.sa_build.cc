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

// Print: rank, SA[rank], suffix text (tab-separated)
void print_suffix_array_view(const string& T, const vector<int>& sa) {
  cout << "rank\tSA\tsuffix\n";
  for (int r = 0; r < (int)sa.size(); ++r) {
    int i = sa[r];
    cout << r << '\t' << i << '\t';
    for (int j = i; j < (int)T.size(); ++j) cout << T[j];
    cout << '\n';
  }
}

int main() {
  string T = "ACGTACGTACGTACGT";

  vector<int> sa = build_suffix_array(T);

  // Print raw suffix array
  cout << "SA: ";
  for (int i = 0; i < (int)sa.size(); ++i) {
    if (i) cout << ' ';
    cout << sa[i];
  }
  cout << endl;

  // Pretty print
  print_suffix_array_view(T, sa);

  return 0;
}