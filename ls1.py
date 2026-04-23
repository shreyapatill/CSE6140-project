import random
import time
import math


def _two_approx_cover(edges, rng):
    # greedy matching for a decent starting point
    cover = set()
    covered = set()
    
    order = list(range(len(edges)))
    rng.shuffle(order)
    
    for i in order:
        if i in covered:
            continue
        u, v = edges[i]
        cover.add(u)
        cover.add(v)
        
        # mark edges touching u or v as covered
        for j, (a, b) in enumerate(edges):
            if a == u or a == v or b == u or b == v:
                covered.add(j)
    
    return cover


def _greedy_shrink(cover, adj):
    # remove nodes we don't actually need
    cover = set(cover)
    changed = True
    while changed:
        changed = False
        
        # try to remove nodes that cover the fewest unique edges first
        candidates = list(cover)
        candidates.sort(key=lambda v: sum(1 for u in adj[v] if u not in cover))
        
        for v in candidates:
            # can we remove v? only if all its neighbors are already in the cover
            needed = False
            for u in adj[v]:
                if u not in cover:
                    needed = True
                    break
            
            if not needed:
                cover.remove(v)
                changed = True
    return cover


def solve_ls1(n, adj, edges, cutoff, seed):
    # setup
    rng = random.Random(seed)
    start_time = time.time()
    m = len(edges)
    
    # get initial cover
    cover = _two_approx_cover(edges, rng)
    cover = _greedy_shrink(cover, adj)
    
    best_cover = set(cover)
    best_size = len(best_cover)
    
    elapsed = time.time() - start_time
    trace = [(elapsed, best_size)]
    
    # simulated annealing params
    T = max(1.0, best_size * 0.05)
    alpha = 0.9999
    T_min = 0.0001
    
    no_improve = 0
    
    while True:
        elapsed = time.time() - start_time
        if elapsed >= cutoff:
            break
            
        cover_list = list(cover)
        if not cover_list:
            # shouldn't happen, but just in case we empty the cover
            cover = set(best_cover)
            cover_list = list(cover)
            
        # pick a few random nodes, and select the one
        # with the fewest neighbors outside the cover.
        candidates = rng.sample(cover_list, min(5, len(cover_list)))
        rem_v = min(candidates, key=lambda v: sum(1 for u in adj[v] if u not in cover))
        
        old_cover = set(cover)
        cover.remove(rem_v)
        
        # fix any edges we just uncovered
        added = set()
        for u in adj[rem_v]:
            if u not in cover:
                # edge (rem_v, u) needs to be covered
                # only u can cover it since we just removed rem_v
                cover.add(u)
                added.add(u)
        
        new_size = len(cover)
        
        # evaluate our new cover
        if new_size < best_size:
            # better! let's try to shrink it more
            cover = _greedy_shrink(cover, adj)
            new_size = len(cover)
            
            best_cover = set(cover)
            best_size = new_size
            elapsed = time.time() - start_time
            trace.append((elapsed, best_size))
            no_improve = 0
            
        elif new_size == best_size:
            # sideways move, accept it to explore
            no_improve += 1
            
        elif rng.random() < math.exp(-(new_size - best_size) / max(T, T_min)):
            # worse move, but SA says accept
            no_improve += 1
            
        else:
            # reject move, go back to old cover
            cover = old_cover
            no_improve += 1
            
        # cool down
        T *= alpha
        
        # restart if we're stuck for a while
        if no_improve > 200:
            no_improve = 0
            T = max(1.0, best_size * 0.05) # heat back up
            
            # restart from our best solution
            cover = set(best_cover)
            
            # randomly drop some nodes to shake things up
            drops = max(1, best_size // 5)
            for v in rng.sample(list(cover), min(drops, len(cover))):
                cover.discard(v)
                
            # add nodes back until we have a valid cover again
            for u, v in edges:
                if u not in cover and v not in cover:
                    # pick the one with more edges
                    add_v = u if len(adj[u]) >= len(adj[v]) else v
                    cover.add(add_v)
                    
            # shrink it
            cover = _greedy_shrink(cover, adj)
            
    return best_cover, trace
