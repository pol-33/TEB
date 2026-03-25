NUM_RECURSIVE_CALLS = -1

def edit_distance_recursive(x, y, i, j):

    global NUM_RECURSIVE_CALLS
    NUM_RECURSIVE_CALLS += 1
    print("CALL:",x[i:], y[j:])
    if i == 0: return j
    if j == 0: return i
    cost = 0 if x[i-1] == y[j-1] else 1
    return min(edit_distance_recursive(x, y, i-1, j-1) + cost,
                edit_distance_recursive(x, y, i-1, j) + 1,
                edit_distance_recursive(x, y, i, j-1) + 1
                )

def edit_distance_dp(x, y):
    m = len(x)
    n = len(y)
    dp = [[0 for _ in range(n+1)] for _ in range(m+1)]
    for i in range(m+1):
        dp[i][0] = i
    for j in range(n+1):
        dp[0][j] = j
    for i in range(1, m+1):
        for j in range(1, n+1):
            cost = 0 if x[i-1] == y[j-1] else 1
            dp[i][j] = min(dp[i-1][j-1] + cost,
                           dp[i-1][j] + 1,
                           dp[i][j-1] + 1)
    return dp[m][n]

if __name__ == "__main__":
    x = "ACTT"
    y = "AGTGT"
    print(edit_distance_recursive(x, y, len(x), len(y)))
    print("Number of recursive calls: ", NUM_RECURSIVE_CALLS)