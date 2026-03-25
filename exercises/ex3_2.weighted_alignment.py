"""
Global alignment with weighted scoring and linear gap penalty.

Given two DNA sequences x and y, compute:
  - the full dynamic programming matrix
  - one optimal global alignment by backtrace
  - the optimal path in the DP matrix

Scoring:
  - match: +2
  - transition mismatch: -1
  - transversion mismatch: -2
  - gap: -2

Time complexity:
  - DP construction: O(n * m)
  - Backtrace: O(n + m)

Space complexity:
  - DP matrix: O(n * m)
"""

from typing import Dict, List, Tuple
from dp_plot import plot_dp_matrix

GAP_PENALTY = -2

SUBSTITUTION_MATRIX: Dict[str, Dict[str, int]] = {
  "A": {"A":  2, "C": -2, "G": -1, "T": -2},
  "C": {"A": -2, "C":  2, "G": -2, "T": -1},
  "G": {"A": -1, "C": -2, "G":  2, "T": -2},
  "T": {"A": -2, "C": -1, "G": -2, "T":  2},
}


def substitution_score(a: str, b: str) -> int:
  """
  Return the substitution score between two nucleotides.

  Time: O(1)
  Space: O(1)
  """
  return SUBSTITUTION_MATRIX[a][b]


def compute_dp_matrix(x: str, y: str) -> List[List[int]]:
  """
  Compute the full DP matrix for weighted global alignment.

  dp[i][j] = best alignment score between x[:i] and y[:j].

  Time: O(n * m)
  Space: O(n * m)
  """
  n = len(x)
  m = len(y)

  dp = [[0] * (m + 1) for _ in range(n + 1)]

  for i in range(1, n + 1):
    dp[i][0] = dp[i - 1][0] + GAP_PENALTY

  for j in range(1, m + 1):
    dp[0][j] = dp[0][j - 1] + GAP_PENALTY

  for i in range(1, n + 1):
    for j in range(1, m + 1):
      dp[i][j] = max(
        dp[i - 1][j - 1] + substitution_score(x[i - 1], y[j - 1]),
        dp[i - 1][j] + GAP_PENALTY,
        dp[i][j - 1] + GAP_PENALTY
      )

  return dp


def backtrace_alignment(
  x: str,
  y: str,
  dp: List[List[int]]
) -> Tuple[str, str, List[Tuple[int, int]]]:
  """
  Reconstruct one optimal global alignment from the DP matrix.

  Returns:
    aligned_x: aligned version of x
    aligned_y: aligned version of y
    path: list of visited DP cells from (0, 0) to (n, m)

  Time: O(n + m)
  Space: O(n + m)
  """
  i = len(x)
  j = len(y)

  aligned_x = []
  aligned_y = []
  path = [(i, j)]

  while i > 0 or j > 0:
    if i > 0 and j > 0:
      score_diag = dp[i - 1][j - 1] + substitution_score(x[i - 1], y[j - 1])
      if dp[i][j] == score_diag:
        aligned_x.append(x[i - 1])
        aligned_y.append(y[j - 1])
        i -= 1
        j -= 1
        path.append((i, j))
        continue

    if i > 0:
      score_up = dp[i - 1][j] + GAP_PENALTY
      if dp[i][j] == score_up:
        aligned_x.append(x[i - 1])
        aligned_y.append("-")
        i -= 1
        path.append((i, j))
        continue

    if j > 0:
      score_left = dp[i][j - 1] + GAP_PENALTY
      if dp[i][j] == score_left:
        aligned_x.append("-")
        aligned_y.append(y[j - 1])
        j -= 1
        path.append((i, j))
        continue

    raise RuntimeError("Backtrace failed: inconsistent DP matrix.")

  aligned_x.reverse()
  aligned_y.reverse()
  path.reverse()

  return "".join(aligned_x), "".join(aligned_y), path


def global_alignment(
  x: str,
  y: str
) -> Tuple[int, List[List[int]], str, str, List[Tuple[int, int]]]:
  """
  Compute the best global alignment score and one optimal alignment.

  Time: O(n * m)
  Space: O(n * m)
  """
  dp = compute_dp_matrix(x, y)
  aligned_x, aligned_y, path = backtrace_alignment(x, y, dp)
  score = dp[len(x)][len(y)]
  return score, dp, aligned_x, aligned_y, path


def print_dp_matrix(dp: List[List[int]], x: str, y: str) -> None:
  """
  Pretty print of the DP matrix with fixed-width columns.

  Time: O(n * m)
  Space: O(1) extra
  """
  width = 4

  print(" " * (width + 2), end="")
  print(f"{'-':>{width}}", end=" ")
  for c in y:
    print(f"{c:>{width}}", end=" ")
  print()

  for i in range(len(dp)):
    row_label = "-" if i == 0 else x[i - 1]
    print(f"{row_label:>2} ", end="")
    for value in dp[i]:
      print(f"{value:>{width}}", end=" ")
    print()


def main() -> None:
  x = "ACTATTTACGTACT"
  y = "ACGTACGTACGAATACGT"

  score, dp, aligned_x, aligned_y, path = global_alignment(x, y)

  print("Sequence X:", x)
  print("Sequence Y:", y)
  print()
  print("Optimal score:", score)
  print()
  print("Optimal alignment:")
  print(aligned_x)
  print(aligned_y)
  print()
  print("DP matrix:")
  print_dp_matrix(dp, x, y)
  print()
  print("Optimal path:")
  print(path)

  plot_dp_matrix(dp, x, y, path)


if __name__ == "__main__":
  main()