# An Adaptive Sampling Algorithm for the Top-K Group Betweenness Centrality
 This repository implements the Adaptive Sampling Algorithm in our paper "An Adaptive Sampling Algorithm for the Top-K Group Betweenness Centrality". \
 We proposed a novel Adaptive Sampling Algorithm to solve the Top-K Group Betweenness Centrality problem, which is to find a group of K nodes from a network so that the total fraction of shortest paths that pass through the K nodes is maximized. 
# Installation
 The algorithm is completed by c++, and the software requires *Microsoft Visual Studio 2023*.
# Running
 The folder `include` contains the header files supporting our algorithm and you can use `main.cpp` to run our algorithm.
# Datasets
 We adopt eight real-world networks and two synthetic networks, which can be found in `datasets`.
# Results
 Experimental results with real-world networks demonstrate that the numer of samples used by the proposed algorithm is from 2 to 18 times smaller than the state-of-the-art, while the centrality of the group found by the algorithm is comparable with the baseline.
