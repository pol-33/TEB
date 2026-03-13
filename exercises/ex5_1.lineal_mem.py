"""
Global alignment with Edit Distance in linear space.

Given two sequences x and y, compute only the optimal edit distance score
using dynamic programming with two rows.

Cost model:
  - match: 0
  - mismatch: 1
  - insertion: 1
  - deletion: 1

Time complexity:
  - DP construction: O(n * m)

Space complexity:
  - Two DP rows: O(m)
"""

def global_alignment_linear_space(x: str, y: str) -> int:
  """
  Compute the edit distance using only two DP rows.

  dp[i][j] depends only on:
    - dp[i - 1][j]     (previous row, same column)
    - dp[i][j - 1]     (current row, previous column)
    - dp[i - 1][j - 1] (previous row, previous column)

  Therefore, we only need:
    - prev_row: DP row i - 1
    - curr_row: DP row i

  Returns:
    The optimal global alignment score, i.e. the edit distance.

  Time: O(n * m)
  Space: O(m)
  """
  n = len(x)
  m = len(y)

  prev_row = list(range(m + 1))
  curr_row = [0] * (m + 1)

  for i in range(1, n + 1):
    curr_row[0] = i

    for j in range(1, m + 1):
      substitution_cost = 0 if x[i - 1] == y[j - 1] else 1

      curr_row[j] = min(
        prev_row[j] + 1,                 # deletion
        curr_row[j - 1] + 1,             # insertion
        prev_row[j - 1] + substitution_cost
      )

    prev_row, curr_row = curr_row, prev_row

  return prev_row[m]


def main() -> None:
  x = "ACTATTTACGTACT"
  y = "ACGTACGTACGAATACGT"

  distance = global_alignment_linear_space(x, y)

  print("Sequence X:", x)
  print("Sequence Y:", y)
  print()
  print("Edit distance:", distance)


if __name__ == "__main__":
  main()