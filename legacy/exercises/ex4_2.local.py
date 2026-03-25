"""
Local alignment with affine gap penalties (Smith-Waterman-Gotoh).

We compute:
  - M[i][j]: best local alignment score ending with x[i-1] aligned to y[j-1]
  - Ix[i][j]: best local alignment score ending with a gap in y
  - Iy[i][j]: best local alignment score ending with a gap in x

Scoring:
  - match: +2
  - transition mismatch: -1
  - transversion mismatch: -2
  - gap_open: -3
  - gap_extend: -1

Local alignment rules:
  - alignments may start anywhere
  - alignments may end anywhere
  - negative scores are replaced by 0

Time complexity:
  - DP construction: O(n * m)
  - Backtrace: O(n + m)

Space complexity:
  - O(n * m)
"""

from typing import Dict, List, Tuple

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
  Compute the three DP matrices for local affine alignment.

  Time: O(n * m)
  Space: O(n * m)
  """
  n = len(x)
  m = len(y)

  M = [[0] * (m + 1) for _ in range(n + 1)]
  Ix = [[0] * (m + 1) for _ in range(n + 1)]
  Iy = [[0] * (m + 1) for _ in range(n + 1)]

  for i in range(1, n + 1):
    for j in range(1, m + 1):
      s = substitution_score(x[i - 1], y[j - 1])

      M[i][j] = max(
        0,
        M[i - 1][j - 1] + s,
        Ix[i - 1][j - 1] + s,
        Iy[i - 1][j - 1] + s
      )

      Ix[i][j] = max(
        0,
        M[i - 1][j] + GAP_OPEN,
        Ix[i - 1][j] + GAP_EXTEND
      )

      Iy[i][j] = max(
        0,
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
  Reconstruct one optimal local alignment.

  The alignment starts from the best scoring cell anywhere and stops when
  the current state reaches score 0.

  Time: O(n + m)
  Space: O(n + m)
  """
  n = len(x)
  m = len(y)

  best_score = 0
  best_i = 0
  best_j = 0
  best_state = "M"

  for i in range(n + 1):
    for j in range(m + 1):
      for state_name, value in (
        ("M", M[i][j]),
        ("Ix", Ix[i][j]),
        ("Iy", Iy[i][j]),
      ):
        if value > best_score:
          best_score = value
          best_i = i
          best_j = j
          best_state = state_name

  i = best_i
  j = best_j
  state = best_state

  aligned_x = []
  aligned_y = []
  path = [(i, j)]

  while i > 0 and j > 0:
    if state == "M":
      if M[i][j] == 0:
        break

      s = substitution_score(x[i - 1], y[j - 1])

      if M[i][j] == M[i - 1][j - 1] + s:
        aligned_x.append(x[i - 1])
        aligned_y.append(y[j - 1])
        i -= 1
        j -= 1
        state = "M"
        path.append((i, j))
        continue

      if M[i][j] == Ix[i - 1][j - 1] + s:
        aligned_x.append(x[i - 1])
        aligned_y.append(y[j - 1])
        i -= 1
        j -= 1
        state = "Ix"
        path.append((i, j))
        continue

      if M[i][j] == Iy[i - 1][j - 1] + s:
        aligned_x.append(x[i - 1])
        aligned_y.append(y[j - 1])
        i -= 1
        j -= 1
        state = "Iy"
        path.append((i, j))
        continue

    elif state == "Ix":
      if Ix[i][j] == 0:
        break

      if Ix[i][j] == M[i - 1][j] + GAP_OPEN:
        aligned_x.append(x[i - 1])
        aligned_y.append("-")
        i -= 1
        state = "M"
        path.append((i, j))
        continue

      if Ix[i][j] == Ix[i - 1][j] + GAP_EXTEND:
        aligned_x.append(x[i - 1])
        aligned_y.append("-")
        i -= 1
        state = "Ix"
        path.append((i, j))
        continue

    else:  # state == "Iy"
      if Iy[i][j] == 0:
        break

      if Iy[i][j] == M[i][j - 1] + GAP_OPEN:
        aligned_x.append("-")
        aligned_y.append(y[j - 1])
        j -= 1
        state = "M"
        path.append((i, j))
        continue

      if Iy[i][j] == Iy[i][j - 1] + GAP_EXTEND:
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


def local_alignment(
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
  Compute local affine alignment.

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
      cell = str(int(value))
      print(f"{cell:>{width}}", end=" ")
    print()


def main() -> None:
  # Both sequences contain a strong common infix.
  x = "TTTACGTCGATTT"
  y = "GGGACGTCGACCC"

  score, M, Ix, Iy, aligned_x, aligned_y, path = local_alignment(x, y)

  print("Sequence X:", x)
  print("Sequence Y:", y)
  print()
  print("Optimal local score:", int(score))
  print()
  print("Optimal local alignment:")
  print(aligned_x)
  print(aligned_y)
  print()
  print("Matrix M:")
  print_dp_matrix(M, x, y)

  try:
    from dp_plot import plot_dp_matrix
    plot_dp_matrix(
      M,
      x,
      y,
      path,
      title="Local Affine Alignment: Matrix M"
    )
  except ModuleNotFoundError:
    pass


if __name__ == "__main__":
  main()