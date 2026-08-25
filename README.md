# HP-MC

## Description

HP-MC is a C++/OpenMP code for quickly computing the exact maximum clique size of large sparse graphs. It uses multiple techniques such as branching on pivot vertices, pruning induced connected components, iterative coloring for tighter bounds on induced subgraphs and work-stealing parallelization. For more information on the implementation, see the [Publication](#publication) section below. 

## Requirements
- GCC 11 or later (tested with GCC 13.3)
- C++17 or later
- OpenMP support

## Installation
Follow these instructions to install, compile, and run the code.

1. Clone the code using "git clone https://github.com/burtscher/HP-MC.git"
2. Navigate to the HP-MC directory with "cd HP-MC"
3. Execute "make"
4. Run "./hpmc graphs/internet.egr"

To convert your own graphs into our format, follow the instructions at https://userweb.cs.txstate.edu/~burtscher/research/ECLgraph/.

## Publication

If you use HP-MC in your work, please cite the following publication:
Cameron Bradley and Martin Burtscher. HP-MC: Quickly Computing Maximum Clique Sizes of Large Sparse Graphs. Proceedings of the IEEE International Conference on Cluster Computing. September 2026. <br> 

[[paper]](https://userweb.cs.txstate.edu/~burtscher/papers/cluster26b.pdf)

This work has been supported by the U.S. Department of Education under a Graduate Assistance in Areas of National Need (GAANN) award.
