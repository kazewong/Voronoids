# Voronoids
A parallel code to compute the Voronoi diagram in Rust.

This code implements parallel insertion algorithms in 2D and 3D to construct the Delaunay graph, which is the dual graph of the Voronoi diagram.

## Installation

### Rust

The rust crate is available on crates.io, and it can be installed by the command

```
cargo add voronoids
```

### Python wrapper

We provide a python wrapper to the rust library that can be installed by the command

```
pip install voronoids
```

Note that the Python wrapper does not expose all the functionality in Rust, and is slightly less performant compared to the Rust version of the code because we need to wrap the output of the Rust library into data type in Python.
To construct a Delaunay graph from a number of points, follow the snippet below:

```
import numpy as np
import voronoids

pts = np.random.uniform(size=(10000,3))
delaunay_graph = voronoids.delaunay(pts)
```

## Attribution

The algorithm implemented in this code follows these two following papers:

```
@article{BOISSONNAT1993339,
title = {On the randomized construction of the Delaunay tree},
journal = {Theoretical Computer Science},
volume = {112},
number = {2},
pages = {339-354},
year = {1993},
issn = {0304-3975},
doi = {https://doi.org/10.1016/0304-3975(93)90024-N},
url = {https://www.sciencedirect.com/science/article/pii/030439759390024N},
author = {Jean-Daniel Boissonnat and Monique Teillaud},
abstract = {The Delaunay tree is a hierarchical data structure which is defined from the Delaunay triangulation and, roughly speaking, represents a triangulation as a hierarchy of balls. It allows a smidynamic construction of the Delaunay triangulation of a finite set of n points in any dimension. In this paper, we prove that a randomized construction of the Delaunay tree (and, thus, of the Delaunay triangulation) can be done in O(n log n) expected time in the plane and in O(n⌈d2⌉) expected time in d-dimensional space. These results are optimal for fixed d. The algorithm is extremely simple and experimental results are given.}
}

@article{10.1145/3402819,
author = {Blelloch, Guy E. and Gu, Yan and Shun, Julian and Sun, Yihan},
title = {Parallelism in Randomized Incremental Algorithms},
year = {2020},
issue_date = {October 2020},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
volume = {67},
number = {5},
issn = {0004-5411},
url = {https://doi.org/10.1145/3402819},
doi = {10.1145/3402819},
abstract = {In this article, we show that many sequential randomized incremental algorithms are in fact parallel. We consider algorithms for several problems, including Delaunay triangulation, linear programming, closest pair, smallest enclosing disk, least-element lists, and strongly connected components.We analyze the dependencies between iterations in an algorithm and show that the dependence structure is shallow with high probability or that, by violating some dependencies, the structure is shallow and the work is not increased significantly. We identify three types of algorithms based on their dependencies and present a framework for analyzing each type. Using the framework gives work-efficient polylogarithmic-depth parallel algorithms for most of the problems that we study.This article shows the first incremental Delaunay triangulation algorithm with optimal work and polylogarithmic depth. This result is important, since most implementations of parallel Delaunay triangulation use the incremental approach. Our results also improve bounds on strongly connected components and least-element lists and significantly simplify parallel algorithms for several problems.},
journal = {J. ACM},
month = {sep},
articleno = {27},
numpages = {27},
keywords = {strongly connected components, smallest enclosing disk, linear programming, least-element lists, closest pair, Randomized incremental algorithms, Delaunay triangulation}
}
```
