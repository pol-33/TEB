"""
Semi-global alignment with affine gap penalties.

We align the whole short sequence y against any substring of the long sequence x.
Therefore:
  - prefix gaps in x are free
  - suffix gaps in x are free
  - y must be fully aligned

Matrices:
  - M[i][j]: best score ending with x[i-1] aligned to y[j-1]
  - Ix[i][j]: best score ending with a gap in y
  - Iy[i][j]: best score ending with a gap in x

Scoring:
  - match: +2
  - transition mismatch: -1
  - transversion mismatch: -2
  - gap_open: -3
  - gap_extend: -1

Time complexity:
  - DP construction: O(n * m)
  - Backtrace: O(n + m)

Space complexity:
  - O(n * m)
"""

from typing import Dict, List, Tuple
from dp_plot import plot_dp_matrix

NEG_INF = float("-inf")

GAP_OPEN = -3
GAP_EXTEND = -1

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


def compute_dp_matrices(
  x: str,
  y: str
) -> Tuple[List[List[float]], List[List[float]], List[List[float]]]:
  """
  Compute the three DP matrices for semi-global affine alignment.

  x is the long sequence.
  y is the short sequence that must be fully aligned.

  Time: O(n * m)
  Space: O(n * m)
  """
  n = len(x)
  m = len(y)

  M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
  Ix = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
  Iy = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

  M[0][0] = 0

  # Free prefix in x
  for i in range(1, n + 1):
    M[i][0] = 0

  # Prefix in y is not free
  for j in range(1, m + 1):
    Iy[0][j] = GAP_OPEN + (j - 1) * GAP_EXTEND

  for i in range(1, n + 1):
    for j in range(1, m + 1):
      s = substitution_score(x[i - 1], y[j - 1])

      M[i][j] = max(
        M[i - 1][j - 1],
        Ix[i - 1][j - 1],
        Iy[i - 1][j - 1]
      ) + s

      Ix[i][j] = max(
        M[i - 1][j] + GAP_OPEN,
        Ix[i - 1][j] + GAP_EXTEND
      )

      Iy[i][j] = max(
        M[i][j - 1] + GAP_OPEN,
        Iy[i][j - 1] + GAP_EXTEND
      )

  return M, Ix, Iy


def backtrace_alignment(
  x: str,
  y: str,
  M: List[List[float]],
  Ix: List[List[float]],
  Iy: List[List[float]]
) -> Tuple[float, str, str, List[Tuple[int, int]]]:
  """
  Reconstruct one optimal semi-global alignment.

  The optimal endpoint is the best state in the last column j = len(y).
  Backtrace stops when j = 0, because prefix gaps in x are free.

  Time: O(n + m)
  Space: O(n + m)
  """
  n = len(x)
  m = len(y)

  best_score = NEG_INF
  best_i = 0
  best_state = "M"

  for i in range(n + 1):
    for state_name, value in (("M", M[i][m]), ("Ix", Ix[i][m]), ("Iy", Iy[i][m])):
      if value > best_score:
        best_score = value
        best_i = i
        best_state = state_name

  i = best_i
  j = m
  state = best_state

  aligned_x = []
  aligned_y = []
  path = [(i, j)]

  while j > 0:
    if state == "M":
      s = substitution_score(x[i - 1], y[j - 1])

      if i > 0 and j > 0 and M[i][j] == M[i - 1][j - 1] + s:
        aligned_x.append(x[i - 1])
        aligned_y.append(y[j - 1])
        i -= 1
        j -= 1
        state = "M"
        path.append((i, j))
        continue

      if i > 0 and j > 0 and M[i][j] == Ix[i - 1][j - 1] + s:
        aligned_x.append(x[i - 1])
        aligned_y.append(y[j - 1])
        i -= 1
        j -= 1
        state = "Ix"
        path.append((i, j))
        continue

      if i > 0 and j > 0 and M[i][j] == Iy[i - 1][j - 1] + s:
        aligned_x.append(x[i - 1])
        aligned_y.append(y[j - 1])
        i -= 1
        j -= 1
        state = "Iy"
        path.append((i, j))
        continue

    elif state == "Ix":
      if i > 0 and Ix[i][j] == M[i - 1][j] + GAP_OPEN:
        aligned_x.append(x[i - 1])
        aligned_y.append("-")
        i -= 1
        state = "M"
        path.append((i, j))
        continue

      if i > 0 and Ix[i][j] == Ix[i - 1][j] + GAP_EXTEND:
        aligned_x.append(x[i - 1])
        aligned_y.append("-")
        i -= 1
        state = "Ix"
        path.append((i, j))
        continue

    else:  # state == "Iy"
      if j > 0 and Iy[i][j] == M[i][j - 1] + GAP_OPEN:
        aligned_x.append("-")
        aligned_y.append(y[j - 1])
        j -= 1
        state = "M"
        path.append((i, j))
        continue

      if j > 0 and Iy[i][j] == Iy[i][j - 1] + GAP_EXTEND:
        aligned_x.append("-")
        aligned_y.append(y[j - 1])
        j -= 1
        state = "Iy"
        path.append((i, j))
        continue

    raise RuntimeError("Backtrace failed: inconsistent DP matrices.")

  aligned_x.reverse()
  aligned_y.reverse()
  path.reverse()

  return best_score, "".join(aligned_x), "".join(aligned_y), path


def semiglobal_alignment(
  x: str,
  y: str
) -> Tuple[
  float,
  List[List[float]],
  List[List[float]],
  List[List[float]],
  str,
  str,
  List[Tuple[int, int]]
]:
  """
  Compute semi-global affine alignment.

  Time: O(n * m)
  Space: O(n * m)
  """
  M, Ix, Iy = compute_dp_matrices(x, y)
  score, aligned_x, aligned_y, path = backtrace_alignment(x, y, M, Ix, Iy)
  return score, M, Ix, Iy, aligned_x, aligned_y, path


def print_dp_matrix(dp: List[List[float]], x: str, y: str) -> None:
  """
  Pretty print of a DP matrix with fixed-width columns.

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
      if value == NEG_INF:
        cell = "-inf"
      else:
        cell = str(int(value))
      print(f"{cell:>{width}}", end=" ")
    print()


def matrix_M_for_plot(M: List[List[float]]) -> List[List[float]]:
  """
  Convert matrix M into a plottable matrix.

  NEG_INF entries are replaced by a finite value below the minimum finite score.

  Time: O(n * m)
  Space: O(n * m)
  """
  finite_values = [
    value
    for row in M
    for value in row
    if value != NEG_INF
  ]

  min_finite = min(finite_values)
  replacement = min_finite - 1

  plot_M = []
  for row in M:
    new_row = []
    for value in row:
      new_row.append(replacement if value == NEG_INF else value)
    plot_M.append(new_row)

  return plot_M


def main() -> None:
  # y appears roughly in the middle of x
  y = "TTTACGTCGATTT"
  x = "ACGTCGA"

  score, M, Ix, Iy, aligned_x, aligned_y, path = semiglobal_alignment(x, y)

  print("Long sequence X:", x)
  print("Short sequence Y:", y)
  print()
  print("Optimal semi-global score:", int(score))
  print()
  print("Optimal alignment:")
  print(aligned_x)
  print(aligned_y)
  print()
  print("Matrix M:")
  print_dp_matrix(M, x, y)

  plot_dp_matrix(
    matrix_M_for_plot(M),
    x,
    y,
    path,
    title="Semi-global Affine Alignment: Matrix M"
  )


if __name__ == "__main__":
  main()