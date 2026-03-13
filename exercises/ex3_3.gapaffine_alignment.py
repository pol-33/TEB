"""
Global alignment with affine gap penalties (Gotoh algorithm).

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

Gap model:
  gap(k) = gap_open + (k - 1) * gap_extend

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
  Compute the three DP matrices for affine-gap global alignment.

  Time: O(n * m)
  Space: O(n * m)
  """
  n = len(x)
  m = len(y)

  M = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
  Ix = [[NEG_INF] * (m + 1) for _ in range(n + 1)]
  Iy = [[NEG_INF] * (m + 1) for _ in range(n + 1)]

  M[0][0] = 0

  for i in range(1, n + 1):
    Ix[i][0] = GAP_OPEN + (i - 1) * GAP_EXTEND

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
  Reconstruct one optimal global alignment.

  The returned path is projected only on matrix M for plotting purposes.

  Time: O(n + m)
  Space: O(n + m)
  """
  i = len(x)
  j = len(y)

  states = {
    "M": M[i][j],
    "Ix": Ix[i][j],
    "Iy": Iy[i][j],
  }
  state = max(states, key=states.get)
  score = states[state]

  aligned_x = []
  aligned_y = []
  path = [(i, j)]

  while i > 0 or j > 0:
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

  return score, "".join(aligned_x), "".join(aligned_y), path


def affine_global_alignment(
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
  Compute affine-gap global alignment.

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
  x = "ACTATTTACGTACT"
  y = "ACGTACGTACGAATACGT"

  score, M, Ix, Iy, aligned_x, aligned_y, path = affine_global_alignment(x, y)

  print("Sequence X:", x)
  print("Sequence Y:", y)
  print()
  print("Optimal score:", int(score))
  print()
  print("Optimal alignment:")
  print(aligned_x)
  print(aligned_y)
  print()
  print("Matrix M:")
  print_dp_matrix(M, x, y)

  # Only matrix M is printed and plotted.
  plot_dp_matrix(
    matrix_M_for_plot(M),
    x,
    y,
    path,
    title="Affine Gap Global Alignment: Matrix M"
  )

if __name__ == "__main__":
  main()