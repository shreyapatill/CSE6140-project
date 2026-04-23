"""
ls2.py -- Local Search 2 for Minimum Vertex Cover.

Algorithm: Memetic Algorithm with NuMVC inner local search and
Adaptive Operator Selection.  We refer to it as MA-NuMVC.

Memetic algorithms hybridise population-based evolutionary search with
single-individual local search [Moscato 1989]. The classical recipe is:
the genetic component generates structurally diverse candidate
solutions via crossover and mutation; the local-search component then
"refines" each offspring to the nearest local optimum before it competes
for a population slot. The local search and the evolutionary operators
play complementary roles: crossover provides global diversification,
local search provides intensification.

For Minimum Vertex Cover and the closely related Maximum Independent
Set, two streams of state-of-the-art results inform our design:

    (a) Memetic / population-based methods [Pullan & Hoos 2006;
        Bouamama, Blum, Boukerram 2012] established the framework we
        adopt at the outer level.

    (b) NuMVC [Cai, Su, Luo, Sattar 2013] is the most-cited modern
        local search for MVC, using edge weighting plus configuration
        checking inside a 2-exchange neighbourhood. We embed a
        budget-bounded variant of NuMVC as the local-search step
        applied to every offspring.

On top of this base recipe we incorporate three additional techniques
drawn from the meta-heuristics literature:

    (i)   Adaptive Operator Selection (AOS) via probability matching
          [Fialho, Da Costa, Schoenauer, Sebag 2010]: the algorithm
          maintains a recent-success score for each crossover operator
          and biases selection toward whichever operator has been more
          productive lately. This removes a hyper-parameter the user
          would otherwise have to tune per instance.

    (ii)  Diversity-preserving steady-state replacement: an offspring
          is only admitted to the population if it is not an exact
          duplicate of an existing member, mitigating premature
          convergence on dense instances.

    (iii) Elite-preserving restart on stagnation: when no improvement
          is observed for a configurable number of generations, all
          but the top-k individuals are regenerated from scratch via
          randomised construction. This combines the breadth of
          multi-start methods with the exploitation of an evolutionary
          core.

This contrasts deliberately with LS1 (simulated annealing) on every
algorithmic axis:

    axis             | LS1 (SA)                    | LS2 (MA-NuMVC)
    -----------------+-----------------------------+----------------------
    state            | single trajectory           | population of P
    move             | 1-flip + cascading repair   | crossover + mutation
                     |                             | + 2-exchange refinement
    selection        | Metropolis (probabilistic)  | tournament + worst-
                     | with cooling temperature    | replacement
    memory           | none                        | population diversity
                     |                             | + AOS reward history
                     |                             | + edge-weight memory
                     |                             | inside NuMVC
    perturbation     | random kick after 200       | mutation + elite-
                     | stagnant iterations         | preserving restart

References (cited in the report):
    P. Moscato (1989). On evolution, search, optimization, genetic
        algorithms and martial arts: Towards memetic algorithms.
        Caltech Concurrent Computation Program, Report 826.
    W. Pullan, H. H. Hoos (2006). Dynamic local search for the maximum
        clique problem. Journal of Artificial Intelligence Research 25,
        159-185.
    S. Bouamama, C. Blum, A. Boukerram (2012). A population-based
        iterated greedy algorithm for the minimum weight vertex cover
        problem. Applied Soft Computing 12(6), 1632-1639.
    A. Fialho, L. Da Costa, M. Schoenauer, M. Sebag (2010). Analyzing
        bandit-based adaptive operator selection mechanisms. Annals of
        Mathematics and Artificial Intelligence 60, 25-64.
    S. Cai, K. Su, C. Luo, A. Sattar (2013). NuMVC: An Efficient Local
        Search Algorithm for Minimum Vertex Cover. JAIR 46, 687-716.
"""

import random
import time


# ----- Hyper-parameters (defaults tuned for the project's dataset) -----
_POP_SIZE         = 12      # population size
_NUMVC_STEPS      = 400     # NuMVC steps applied to each offspring / init
_MUTATION_RATE    = 0.30    # probability of mutating an offspring
_STAGNATION_LIMIT = 200     # generations with no improvement -> restart
_ELITE_KEEP       = 3       # elites preserved across a restart
_AOS_DECAY        = 0.95    # AOS reward decay per acceptance event
_TOURNAMENT_K     = 2       # tournament size for parent selection
_BMS              = 50      # NuMVC: best-from-multiple-selections sample


def solve_ls2(n, adj, edges, cutoff, seed):
    """
    Memetic Algorithm with NuMVC inner local search for MVC.

    Args:
        n      : number of vertices (1-indexed; vertices are 1..n)
        adj    : dict[int -> set[int]] adjacency (built by mvc.parse_graph)
        edges  : list of (u, v) pairs
        cutoff : wall-clock time limit in seconds
        seed   : RNG seed (deterministic per seed)

    Returns:
        (best_cover, trace)
            best_cover : set of vertex indices forming a valid vertex cover
            trace      : list of (timestamp_seconds, cover_size) at every
                         improvement, starting with the initial best.
    """
    rng = random.Random(seed)
    t0 = time.time()
    m = len(edges)

    # ----- Per-vertex incidence: list of (neighbour, edge_id) ----------
    adj_e = [[] for _ in range(n + 1)]
    for ei, (u, v) in enumerate(edges):
        adj_e[u].append((v, ei))
        adj_e[v].append((u, ei))
    deg = [len(adj_e[v]) for v in range(n + 1)]

    # ===================================================================
    # Inner local search: budget-bounded NuMVC
    # ===================================================================
    # Implements the 2-exchange + edge-weighting + configuration-checking
    # algorithm of Cai et al. (2013). Each call refines the input cover
    # for at most `budget` steps and returns the smallest cover seen.
    # State is rebuilt per call so that memetic offspring can be refined
    # independently of one another.
    def numvc_refine(C_init, budget):
        # ---- state init ----
        in_cover = [False] * (n + 1)
        for v in C_init:
            in_cover[v] = True

        weight = [1] * m
        uncov = []
        uncov_pos = [-1] * m
        for ei, (u, v) in enumerate(edges):
            if not in_cover[u] and not in_cover[v]:
                uncov_pos[ei] = len(uncov)
                uncov.append(ei)

        # dscore[v] = signed change in uncovered weight if v's
        # membership flips.  v in C: -sum w over neighbours not in C.
        # v not in C: +sum w over neighbours not in C.
        dscore = [0] * (n + 1)
        for v in range(1, n + 1):
            s = 0
            for (u, ei) in adj_e[v]:
                if not in_cover[u]:
                    s += weight[ei]
            dscore[v] = -s if in_cover[v] else s

        cover_list = [v for v in range(1, n + 1) if in_cover[v]]
        cover_pos = [-1] * (n + 1)
        for i, v in enumerate(cover_list):
            cover_pos[v] = i

        conf_change = [1] * (n + 1)
        age = [0] * (n + 1)

        # ---- O(1) helpers ----
        def add_uncov(ei):
            uncov_pos[ei] = len(uncov)
            uncov.append(ei)

        def del_uncov(ei):
            pos = uncov_pos[ei]
            last = uncov.pop()
            if pos != len(uncov):
                uncov[pos] = last
                uncov_pos[last] = pos
            uncov_pos[ei] = -1

        def add_cover_list(v):
            cover_pos[v] = len(cover_list)
            cover_list.append(v)

        def del_cover_list(v):
            pos = cover_pos[v]
            last = cover_list.pop()
            if pos != len(cover_list):
                cover_list[pos] = last
                cover_pos[last] = pos
            cover_pos[v] = -1

        def numvc_remove(v):
            in_cover[v] = False
            del_cover_list(v)
            for (u, ei) in adj_e[v]:
                if in_cover[u]:
                    dscore[u] -= weight[ei]
                else:
                    dscore[u] += weight[ei]
                    conf_change[u] = 1
                    add_uncov(ei)
            dscore[v] = -dscore[v]
            conf_change[v] = 0

        def numvc_add(v):
            in_cover[v] = True
            add_cover_list(v)
            for (u, ei) in adj_e[v]:
                if in_cover[u]:
                    dscore[u] += weight[ei]
                else:
                    dscore[u] -= weight[ei]
                    conf_change[u] = 1
                    del_uncov(ei)
            dscore[v] = -dscore[v]

        # ---- redundancy removal at start ----
        # Vertices with dscore==0 are redundant (no non-C neighbours)
        # so we can drop them for free.
        redundant = [v for v in cover_list if dscore[v] == 0]
        while redundant:
            v = redundant.pop()
            if in_cover[v] and dscore[v] == 0:
                nbrs_in_C = [u for (u, _) in adj_e[v] if in_cover[u]]
                numvc_remove(v)
                for u in nbrs_in_C:
                    if in_cover[u] and dscore[u] == 0:
                        redundant.append(u)

        # ---- best so far ----
        best_C = list(cover_list)
        best_size = len(best_C)

        # ---- main loop ----
        for step in range(budget):
            # If C is a valid cover, drop a vertex to look for a smaller
            # cover; record current best first.
            if not uncov:
                if len(cover_list) < best_size:
                    best_size = len(cover_list)
                    best_C = list(cover_list)
                k = min(_BMS, len(cover_list))
                if k == 0:
                    break
                sample = rng.sample(cover_list, k)
                rm_v = sample[0]
                best_key = (dscore[rm_v], -age[rm_v])
                for c in sample:
                    key = (dscore[c], -age[c])
                    if key > best_key:
                        rm_v = c
                        best_key = key
                age[rm_v] = step
                numvc_remove(rm_v)
                continue

            # 2-exchange: remove highest-dscore vertex from C, then add
            # an endpoint of a random uncovered edge (CC-preferring).
            k = min(_BMS, len(cover_list))
            sample = rng.sample(cover_list, k)
            rm_v = sample[0]
            best_key = (dscore[rm_v], -age[rm_v])
            for c in sample:
                key = (dscore[c], -age[c])
                if key > best_key:
                    rm_v = c
                    best_key = key
            age[rm_v] = step
            numvc_remove(rm_v)

            ei = uncov[rng.randrange(len(uncov))]
            u, v = edges[ei]
            cu, cv = conf_change[u], conf_change[v]
            if cu and cv:
                add_v = u if (dscore[u], -age[u]) > (dscore[v], -age[v]) else v
            elif cu:
                add_v = u
            elif cv:
                add_v = v
            else:
                add_v = u if (dscore[u], -age[u]) > (dscore[v], -age[v]) else v
            age[add_v] = step
            numvc_add(add_v)

            # Bump weights of uncovered edges; both endpoints are not in
            # C so each gets +1 dscore.
            for euid in uncov:
                weight[euid] += 1
                eu, ev = edges[euid]
                dscore[eu] += 1
                dscore[ev] += 1

        # Capture final state if improved.
        if not uncov and len(cover_list) < best_size:
            best_size = len(cover_list)
            best_C = list(cover_list)

        return set(best_C)

    # ===================================================================
    # Genetic operators
    # ===================================================================
    def repair(C):
        """Add high-degree endpoint of any uncovered edge."""
        C = set(C)
        for u, v in edges:
            if u not in C and v not in C:
                C.add(u if deg[u] >= deg[v] else v)
        return C

    def construct():
        """Randomised greedy 2-approximation; refined by NuMVC."""
        order = list(range(m))
        rng.shuffle(order)
        C = set()
        covered = [False] * m
        for ei in order:
            if covered[ei]:
                continue
            u, v = edges[ei]
            if deg[u] > deg[v]:
                pick = u
            elif deg[v] > deg[u]:
                pick = v
            else:
                pick = u if rng.random() < 0.5 else v
            C.add(pick)
            for (_, ej) in adj_e[pick]:
                covered[ej] = True
        return numvc_refine(C, _NUMVC_STEPS)

    def crossover_union(p1, p2):
        """Union the parent covers, then refine.

        Always feasible (a superset of either parent is a valid cover).
        Inherits anchor vertices from both parents but starts oversized;
        NuMVC trims it.
        """
        return numvc_refine(p1 | p2, _NUMVC_STEPS)

    def crossover_common(p1, p2):
        """Keep p1 ∩ p2; randomly inherit half of (p1 △ p2); repair; refine.

        Yields more diverse children than crossover_union because the
        symmetric-difference filter discards redundancy. May produce
        infeasible intermediates, which `repair` then fixes.
        """
        common = p1 & p2
        diff = sorted(p1 ^ p2)        # sort for reproducibility
        child = set(common)
        for v in diff:
            if rng.random() < 0.5:
                child.add(v)
        return numvc_refine(repair(child), _NUMVC_STEPS)

    def mutate(C):
        """Drop k random vertices, then repair + refine.

        k ∼ N(2, 1) clipped to [1, |C|-1]. Different in spirit from LS1's
        "drop ~20% then repair" kick: this is a local perturbation.
        """
        C = set(C)
        if len(C) > 1:
            k = max(1, min(len(C) - 1, int(round(rng.gauss(2.0, 1.0)))))
            for v in rng.sample(sorted(C), k):
                C.discard(v)
        return numvc_refine(repair(C), _NUMVC_STEPS)

    def tournament(pop):
        """Smallest-cover member of a random K-subset."""
        cands = rng.sample(range(len(pop)), _TOURNAMENT_K)
        best = cands[0]
        for c in cands[1:]:
            if len(pop[c]) < len(pop[best]):
                best = c
        return pop[best]

    # ===================================================================
    # Initialisation
    # ===================================================================
    population = [construct() for _ in range(_POP_SIZE)]
    population.sort(key=len)
    best_cover = set(population[0])
    best_size = len(best_cover)
    trace = [(time.time() - t0, best_size)]

    # ----- AOS bookkeeping --------------------------------------------
    op_reward = {'union': 1.0, 'common': 1.0}

    no_improve = 0

    # ===================================================================
    # Main evolutionary loop
    # ===================================================================
    while time.time() - t0 < cutoff:
        # 1. Parent selection (binary tournament)
        p1 = tournament(population)
        p2 = tournament(population)

        # 2. Adaptive crossover-operator selection (probability matching)
        total = op_reward['union'] + op_reward['common']
        if rng.random() < op_reward['union'] / total:
            op_used = 'union'
            child = crossover_union(p1, p2)
        else:
            op_used = 'common'
            child = crossover_common(p1, p2)

        # 3. Optional mutation
        if rng.random() < _MUTATION_RATE:
            child = mutate(child)

        # Track global best as we go, even if not accepted into population
        if len(child) < best_size:
            best_size = len(child)
            best_cover = set(child)
            trace.append((time.time() - t0, best_size))
            no_improve = 0

        # 4. Diversity-aware steady-state replacement
        is_dup = any(child == ind for ind in population)
        worst_idx = max(range(_POP_SIZE), key=lambda k: len(population[k]))

        if (not is_dup) and len(child) < len(population[worst_idx]):
            population[worst_idx] = child
            op_reward[op_used] += 1.0
            op_reward['union'] *= _AOS_DECAY
            op_reward['common'] *= _AOS_DECAY
            if len(child) >= best_size:
                no_improve += 1
        else:
            no_improve += 1

        # 5. Elite-preserving restart on stagnation
        if no_improve >= _STAGNATION_LIMIT:
            no_improve = 0
            population.sort(key=len)
            elite = population[:_ELITE_KEEP]
            population = list(elite)
            while len(population) < _POP_SIZE:
                population.append(construct())
            op_reward = {'union': 1.0, 'common': 1.0}

    return best_cover, trace
