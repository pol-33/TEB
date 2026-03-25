"""
Global alignment with Edit Distance using Hirschberg's algorithm.

Given two sequences x and y, compute:
  - the optimal edit distance
  - one optimal global alignment

Cost model:
  - match: 0
  - mismatch: 1
  - insertion: 1
  - deletion: 1

Hirschberg's algorithm uses divide and conquer:
  1. Split x at its middle position.
  2. Compute, in linear space, the cost of aligning the left half of x
     against every prefix of y.
  3. Compute, in linear space and backwards, the cost of aligning the
     right half of x against every suffix of y.
  4. Find the split point of y that minimizes the total cost.
  5. Recurse on both halves and concatenate the resulting alignments.

Time complexity:
  - O(n * m)

Space complexity:
  - O(n + m) extra space, plus the output alignment
"""

from typing import List, Tuple


def compute_last_row(x: str, y: str) -> List[int]:
  """
  Compute the last DP row for aligning x against y using edit distance,
  keeping only two rows.

  Returns:
    A list row such that row[j] is the edit distance between x and y[:j].

  Time: O(len(x) * len(y))
  Space: O(len(y))
  """
  m = len(y)

  prev_row = list(range(m + 1))
  curr_row = [0] * (m + 1)

  for i in range(1, len(x) + 1):
    curr_row[0] = i

    for j in range(1, m + 1):
      substitution_cost = 0 if x[i - 1] == y[j - 1] else 1

      curr_row[j] = min(
        prev_row[j] + 1,                 # deletion
        curr_row[j - 1] + 1,             # insertion
        prev_row[j - 1] + substitution_cost
      )

    prev_row, curr_row = curr_row, prev_row

  return prev_row


def solve_base_case(x: str, y: str) -> Tuple[str, str]:
  """
  Solve a small global alignment problem by building the full DP matrix.

  This is used only for Hirschberg base cases.

  Returns:
    aligned_x, aligned_y

  Time: O(len(x) * len(y))
  Space: O(len(x) * len(y))
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
        dp[i - 1][j] + 1,
        dp[i][j - 1] + 1,
        dp[i - 1][j - 1] + substitution_cost
      )

  i = n
  j = m
  aligned_x = []
  aligned_y = []

  while i > 0 or j > 0:
    if i > 0 and j > 0:
      substitution_cost = 0 if x[i - 1] == y[j - 1] else 1
      if dp[i][j] == dp[i - 1][j - 1] + substitution_cost:
        aligned_x.append(x[i - 1])
        aligned_y.append(y[j - 1])
        i -= 1
        j -= 1
        continue

    if i > 0 and dp[i][j] == dp[i - 1][j] + 1:
      aligned_x.append(x[i - 1])
      aligned_y.append("-")
      i -= 1
      continue

    aligned_x.append("-")
    aligned_y.append(y[j - 1])
    j -= 1

  aligned_x.reverse()
  aligned_y.reverse()

  return "".join(aligned_x), "".join(aligned_y)


def hirschberg_alignment(x: str, y: str) -> Tuple[str, str]:
  """
  Compute one optimal global alignment using Hirschberg's algorithm.

  Returns:
    aligned_x, aligned_y

  Time: O(len(x) * len(y))
  Space: O(len(x) + len(y)) extra, plus output
  """
  n = len(x)
  m = len(y)

  if n == 0:
    return "-" * m, y

  if m == 0:
    return x, "-" * n

  if n == 1 or m == 1:
    return solve_base_case(x, y)

  mid = n // 2

  x_left = x[:mid]
  x_right = x[mid:]

  score_left = compute_last_row(x_left, y)
  score_right = compute_last_row(x_right[::-1], y[::-1])

  split_j = 0
  best_score = float("inf")

  for j in range(m + 1):
    total_score = score_left[j] + score_right[m - j]
    if total_score < best_score:
      best_score = total_score
      split_j = j

  y_left = y[:split_j]
  y_right = y[split_j:]

  aligned_x_left, aligned_y_left = hirschberg_alignment(x_left, y_left)
  aligned_x_right, aligned_y_right = hirschberg_alignment(x_right, y_right)

  return (
    aligned_x_left + aligned_x_right,
    aligned_y_left + aligned_y_right
  )


def alignment_cost(aligned_x: str, aligned_y: str) -> int:
  """
  Compute the edit distance cost of an already constructed alignment.

  Time: O(L), where L is the alignment length
  Space: O(1)
  """
  cost = 0

  for a, b in zip(aligned_x, aligned_y):
    if a == "-" or b == "-":
      cost += 1
    elif a != b:
      cost += 1

  return cost


def global_alignment_hirschberg(x: str, y: str) -> Tuple[int, str, str]:
  """
  Compute edit distance and one optimal global alignment using Hirschberg.

  Returns:
    distance, aligned_x, aligned_y

  Time: O(n * m)
  Space: O(n + m) extra, plus output
  """
  aligned_x, aligned_y = hirschberg_alignment(x, y)
  distance = alignment_cost(aligned_x, aligned_y)
  return distance, aligned_x, aligned_y


def main() -> None:
  x = "ACTATTTACGTACT"
  y = "ACGTACGTACGAATACGT"

  distance, aligned_x, aligned_y = global_alignment_hirschberg(x, y)

  print("Sequence X:", x)
  print("Sequence Y:", y)
  print()
  print("Edit distance:", distance)
  print()
  print("Optimal alignment:")
  print(aligned_x)
  print(aligned_y)


if __name__ == "__main__":
  main()