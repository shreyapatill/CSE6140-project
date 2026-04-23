def solve_approx(n, adj, edges, cutoff, seed):
    covered = set()  # tracks vertices already in cover
    cover = set()

    for u, v in edges:
        # only act on edges not yet covered
        if u not in covered and v not in covered:
            cover.add(u)
            cover.add(v)
            covered.add(u)
            covered.add(v)

    return cover, None
