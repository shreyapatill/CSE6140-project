#!/usr/bin/env python3
import sys
import argparse
import os

def parse_graph(instance_path):
    # try to find the file
    filepath = instance_path + ".in"
    if not os.path.exists(filepath):
        for folder in ["data/test/", "data/small/", "data/large/"]:
            candidate = folder + os.path.basename(instance_path) + ".in"
            if os.path.exists(candidate):
                filepath = candidate
                break
    
    # read the graph data
    with open(filepath, 'r') as f:
        n, m = map(int, f.readline().split())
        
        adj = {i: set() for i in range(1, n + 1)}
        edges = []
        
        for _ in range(m):
            line = f.readline().strip()
            if not line: continue
            u, v = map(int, line.split())
            adj[u].add(v)
            adj[v].add(u)
            edges.append((u, v))
            
    return n, adj, edges

def write_solution(instance_name, algorithm, cutoff, seed, cover, trace):
    # get just the dataset name, like "small1"
    base_name = os.path.basename(instance_name)
    
    # figure out what to name our files
    if algorithm in ("LS1", "LS2"):
        sol_filename = f"{base_name}_{algorithm}_{cutoff}_{seed}.sol"
        trace_filename = f"{base_name}_{algorithm}_{cutoff}_{seed}.trace"
    elif algorithm == "BnB":
        sol_filename = f"{base_name}_{algorithm}_{cutoff}.sol"
        trace_filename = f"{base_name}_{algorithm}_{cutoff}.trace"
    else:  # Approx
        sol_filename = f"{base_name}_{algorithm}_{cutoff}.sol"
        trace_filename = None
        
    # save the solution
    sorted_cover = sorted(cover)
    with open(sol_filename, 'w') as f:
        f.write(f"{len(sorted_cover)}\n")
        f.write(" ".join(str(v) for v in sorted_cover) + "\n")
        
    # save the trace if we have one
    if trace_filename and trace:
        with open(trace_filename, 'w') as f:
            for timestamp, quality in trace:
                f.write(f"{timestamp:.4f} {quality}\n")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-inst", required=True)
    parser.add_argument("-alg", required=True, choices=["BnB", "Approx", "LS1", "LS2"])
    parser.add_argument("-time", type=int, required=True)
    parser.add_argument("-seed", type=int, required=True)
    args = parser.parse_args()
    
    # load up the graph
    n, adj, edges = parse_graph(args.inst)
    
    # run the chosen algorithm
    if args.alg == "LS1":
        from ls1 import solve_ls1
        cover, trace = solve_ls1(n, adj, edges, args.time, args.seed)
    elif args.alg == "LS2":
        from ls2 import solve_ls2
        cover, trace = solve_ls2(n, adj, edges, args.time, args.seed)
    elif args.alg == "BnB":
        from bnb import solve_bnb
        cover, trace = solve_bnb(n, adj, edges, args.time, args.seed)
    elif args.alg == "Approx":
        from approx import solve_approx
        cover, trace = solve_approx(n, adj, edges, args.time, args.seed)
        
    # write out the results
    write_solution(args.inst, args.alg, args.time, args.seed, cover, trace)

if __name__ == "__main__":
    main()
