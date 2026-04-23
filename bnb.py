import time

def branch_and_bound(adj, n, cutoff):
    start = time.time()

    best_cover = list(adj.keys())
    best_size = n
    trace = []

    def get_uncovered_edges(cover_set):
        uncovered = []
        seen = set()
        for u in adj:
            if u in cover_set:
                continue
            for v in adj[u]:
                if v in cover_set:
                    continue
                edge = (min(u, v), max(u, v))
                if edge not in seen:
                    seen.add(edge)
                    uncovered.append((u, v))
        return uncovered

    def maximal_matching_lb(cover_set):
        """Lower bound: size of maximal matching on uncovered edges."""
        uncovered = get_uncovered_edges(cover_set)
        matched = set()
        count = 0
        for u, v in uncovered:
            if u not in matched and v not in matched:
                matched.add(u)
                matched.add(v)
                count += 1
        return count

    def pick_branch_vertex(cover_set):
        """Pick highest-degree uncovered vertex."""
        uncovered = get_uncovered_edges(cover_set)
        if not uncovered:
            return None
        deg = {}
        for u, v in uncovered:
            deg[u] = deg.get(u, 0) + 1
            deg[v] = deg.get(v, 0) + 1
        return max(deg, key=deg.get)

    # Warm start with greedy cover
    best_cover = list(greedy_upper_bound(adj))
    best_size = len(best_cover)
    trace.append((time.time() - start, best_size))

    # Iterative BnB using an explicit stack
    # Each stack frame: (cover_set, extra_vertices_to_remove_on_pop)
    stack = [(set(), [])]

    while stack:
        if time.time() - start > cutoff:
            break

        cover_set, to_remove = stack.pop()

        # Prune: lower bound check
        lb = maximal_matching_lb(cover_set)
        if len(cover_set) + lb >= best_size:
            continue

        # Check if all edges covered
        uncovered = get_uncovered_edges(cover_set)
        if not uncovered:
            if len(cover_set) < best_size:
                best_size = len(cover_set)
                best_cover = list(cover_set)
                trace.append((time.time() - start, best_size))
            continue

        # Pick vertex to branch on
        v = pick_branch_vertex(cover_set)
        if v is None:
            continue

        # Branch 2: exclude v → must include all its uncovered neighbors
        neighbors = {nb for nb in adj[v] if nb not in cover_set}
        if len(cover_set) + len(neighbors) < best_size:
            new_cover = cover_set | neighbors
            stack.append((new_cover, []))

        # Branch 1: include v (pushed last so explored first)
        new_cover = cover_set | {v}
        stack.append((new_cover, []))

    return best_cover, trace


def greedy_upper_bound(adj):
    """Greedy maximal matching → valid 2-approx vertex cover."""
    cover = set()
    matched = set()
    for u in adj:
        if u in matched:
            continue
        for v in adj[u]:
            if v not in matched:
                cover.add(u)
                cover.add(v)
                matched.add(u)
                matched.add(v)
                break
    return cover