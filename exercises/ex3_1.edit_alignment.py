"""
Global alignment with Edit Distance.

Given two sequences x and y, compute:
  - the full dynamic programming matrix
  - one optimal alignment by backtrace
  - the optimal path in the DP matrix

Cost model:
  - match: 0
  - mismatch: 1
  - insertion: 1
  - deletion: 1

Time complexity:
  - DP construction: O(n * m)
  - Backtrace: O(n + m)

Space complexity:
  - DP matrix: O(n * m)
"""

from typing import List, Tuple
from dp_plot import plot_dp_matrix

def compute_dp_matrix(x: str, y: str) -> List[List[int]]:
  """
  Compute the full edit distance DP matrix.

  dp[i][j] = minimum edit distance between x[:i] and y[:j].

  Time: O(n * m)
  Space: O(n * m)
  """
  n = len(x)
  m = len(y)

  dp = [[0] * (m + 1) for _ in range(n + 1)]

  for i in range(1, n + 1):
    dp[i][0] = i

  for j in range(1, m + 1):
    dp[0][j] = j

  for i in range(1, n + 1):
    for j in range(1, m + 1):
      substitution_cost = 0 if x[i - 1] == y[j - 1] else 1

      dp[i][j] = min(
        dp[i - 1][j] + 1,                    # deletion
        dp[i][j - 1] + 1,                    # insertion
        dp[i - 1][j - 1] + substitution_cost # match / mismatch
      )

  return dp


def backtrace_alignment(
  x: str,
  y: str,
  dp: List[List[int]]
) -> Tuple[str, str, List[Tuple[int, int]]]:
  """
  Reconstruct one optimal alignment from the DP matrix.

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
      substitution_cost = 0 if x[i - 1] == y[j - 1] else 1
      if dp[i][j] == dp[i - 1][j - 1] + substitution_cost:
        aligned_x.append(x[i - 1])
        aligned_y.append(y[j - 1])
        i -= 1
        j -= 1
        path.append((i, j))
        continue

    if i > 0 and dp[i][j] == dp[i - 1][j] + 1:
      aligned_x.append(x[i - 1])
      aligned_y.append("-")
      i -= 1
      path.append((i, j))
      continue

    if j > 0 and dp[i][j] == dp[i][j - 1] + 1:
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
  Compute edit distance and one optimal global alignment.

  Time: O(n * m)
  Space: O(n * m)
  """
  dp = compute_dp_matrix(x, y)
  aligned_x, aligned_y, path = backtrace_alignment(x, y, dp)
  distance = dp[len(x)][len(y)]
  return distance, dp, aligned_x, aligned_y, path

def print_dp_matrix(dp, x, y):
  """
  Pretty print of the DP matrix with fixed-width columns (3 digits).

  Time: O(n * m)
  Space: O(1) extra
  """

  width = 3

  # Header
  print(" " * (width + 2), end="")
  print(f"{'-':>{width}}", end=" ")
  for c in y:
    print(f"{c:>{width}}", end=" ")
  print()

  # Rows
  for i in range(len(dp)):
    row_label = "-" if i == 0 else x[i - 1]
    print(f"{row_label:>2} ", end="")

    for value in dp[i]:
      print(f"{value:>{width}}", end=" ")

    print()


def main() -> None:
  x = "ACTATTTACGTACT"
  y = "ACGTACGTACGAATACGT"

  distance, dp, aligned_x, aligned_y, path = global_alignment(x, y)

  print("Sequence X:", x)
  print("Sequence Y:", y)
  print()
  print("Edit distance:", distance)
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