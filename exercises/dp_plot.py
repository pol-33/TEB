"""
DP matrix visualization for edit distance.

Given a DP matrix, draw:
  - a heatmap of matrix values
  - the numeric value inside each cell
  - optionally, an optimal backtrace path

Time:
  - plotting the matrix values: O(n * m)
  - plotting the path: O(n + m)

Space:
  - O(1) extra, excluding matplotlib internal structures
"""

from typing import List, Tuple, Optional

import matplotlib.pyplot as plt

def plot_dp_matrix(
  dp: List[List[int]],
  x: str,
  y: str,
  path: Optional[List[Tuple[int, int]]] = None,
  title: str = "Edit Distance DP Matrix",
  annotate: bool = True
) -> None:
  """
  Plot a DP matrix as a heatmap and optionally overlay a backtrace path.

  Parameters:
    dp: DP matrix of size (len(x) + 1) x (len(y) + 1)
    x: first sequence
    y: second sequence
    path: optional list of cells (i, j) from (0, 0) to (n, m)
    title: plot title
    annotate: whether to print numeric values inside cells

  Time: O(n * m)
  Space: O(1) extra
  """
  n = len(dp)
  m = len(dp[0])

  # fig, ax = plt.subplots(figsize=(max(6, 0.8 * m), max(5, 0.8 * n)))
  cell_size = 0.6
  fig, ax = plt.subplots(figsize=(max(6, cell_size * m + 2),
                                  max(5, cell_size * n + 2)))

  image = ax.imshow(dp, cmap="RdYlGn", aspect="equal", origin="upper")
  cbar = plt.colorbar(
    image,
    ax=ax,
    shrink=0.7,     # reduce altura
    aspect=25,      # hace la barra más fina
    pad=0.02        # separación del gráfico
  )
  cbar.set_label("DP score", fontsize=10)
  cbar.ax.tick_params(labelsize=8)

  # Axis labels
  x_labels = ["-"] + list(y)
  y_labels = ["-"] + list(x)

  ax.set_xticks(range(m))
  ax.set_yticks(range(n))
  ax.set_xticklabels(x_labels)
  ax.set_yticklabels(y_labels)

  ax.set_xlabel("Sequence X")
  ax.set_ylabel("Sequence Y")
  ax.set_title(title)

  # Grid lines
  ax.set_xticks([j - 0.5 for j in range(1, m)], minor=True)
  ax.set_yticks([i - 0.5 for i in range(1, n)], minor=True)
  ax.grid(which="minor", color="white", linestyle="-", linewidth=1)
  ax.tick_params(which="minor", bottom=False, left=False)

  # Cell annotations
  if annotate:
    max_value = max(max(row) for row in dp)
    threshold = max_value / 2 if max_value > 0 else 0

    for i in range(n):
      for j in range(m):
        text_color = "white" if dp[i][j] > threshold else "black"
        ax.text(
          j,
          i,
          str(dp[i][j]),
          ha="center",
          va="center",
          color=text_color,
          fontsize=10
        )

  # Overlay optimal path
  if path is not None and len(path) > 0:
    xs = [j for i, j in path]
    ys = [i for i, j in path]

    ax.plot(xs, ys, marker="o", linewidth=2)

  plt.tight_layout()
  plt.show()