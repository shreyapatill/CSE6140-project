# CSE6140 Project

Implementation of four algorithms for the Minimum Vertex Cover (MVC) problem.

---

## Team Members

- Shreya Patil
- Siddharth Sawhney
- Isidora Gajic
- Yash Bhalchandra

---

## Project Structure

```
CSE6140-project/
├── mvc.py              # Main entry point — parses arguments and dispatches to algorithms
├── bnb.py              # Exact Branch and Bound algorithm
├── approx.py           # 2-approximation algorithm via maximal matching
├── ls1.py              # Local Search 1
├── ls2.py              # Local Search 2
├── solutions_bnb/      # Output files from Branch and Bound
├── solutions_approx/   # Output files from Approximation
├── solutions_ls1/      # Output files from Local Search 1
├── solutions_ls2/      # Output files from Local Search 2
└── README.md
```

---

## Requirements

- Python 3.8 or higher
- No external libraries required (standard library only)

---

## How to Run

All algorithms are invoked through `mvc.py` using the following command format:

```bash
python mvc.py -inst <instance> -alg [BnB|Approx|LS1|LS2] -time <cutoff> -seed <random seed>
```

### Arguments

| Argument | Description |
|----------|-------------|
| `-inst`  | Path to the dataset instance, without the `.in` extension (e.g. `data/small/small1`) |
| `-alg`   | Algorithm to run: `BnB`, `Approx`, `LS1`, or `LS2` |
| `-time`  | Time limit in seconds |
| `-seed`  | Random seed (used by LS1 and LS2; ignored for BnB and Approx) |


## Output Files

Each run writes its output files to the current working directory (i.e.,
the directory from which `python3 mvc.py ...` is executed). For BnB,
LS1, and LS2 this includes both a `.sol` and a `.trace` file; for
Approx, only a `.sol` file is produced.

### Solution file (`.sol`)

Named `<instance>_<alg>_<cutoff>.sol` for deterministic algorithms, or
`<instance>_<alg>_<cutoff>_<seed>.sol` for randomized algorithms.

```
12
1 4 7 9 10 11 15 18 19 20 22 23
```

Line 1 is the size of the vertex cover. Line 2 is the space-separated list of vertex indices.

### Trace file (`.trace`)

*Produced for `BnB`, `LS1`, and `LS2`. Not produced for `Approx`, since
the approximation algorithm has no anytime improvements to record (per
handout §7.3).*

Named `<instance>_<alg>_<cutoff>.trace` or `<instance>_<alg>_<cutoff>_<seed>.trace`.

```
0.000123 21
0.001456 19
0.012310 18
```

Each line contains a timestamp (seconds elapsed) and the best solution size found at that point.

---

## Algorithm Overview

| Algorithm | Type | Exact? | Seed used? |
|-----------|------|--------|------------|
| `BnB`     | Branch and Bound | Yes | No |
| `Approx`  | Approximation Algorithm | No | No |
| `LS1`     | Local Search 1 | No | Yes |
| `LS2`     | Local Search 2 | No | Yes |

---

## Input Format

Input files follow the `.in` extension. The first line contains two integers `n` and `m`
(number of vertices and edges). Each of the following `m` lines contains two integers
`u` and `v` representing an undirected edge. Vertices are indexed from 1 to `n`.

```
5 6
1 2
1 3
2 3
2 4
3 5
4 5
```
