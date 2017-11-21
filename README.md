# Variable density 3-d nodes for RBF computations
<p align="center">
<img src="https://raw.githubusercontent.com/OVlasiuk/3dRBFnodes/master/nodes.png" width="300">
</p>

---
This is a collection of Matlab routines for distributing nodes (discrete sets) with variable density. The nodes are drawn from a sequence of either periodic Riesz minimizers or irrational lattices, with variable number of elements. The principal application for such nodes are RBF-based PDE solvers. The algorithm processes the working domain by dividing it into a number of equal voxels and putting nodes with constant density in each of them, then applying a repel procedure using the weighted hypersingular Riesz kernel. Currently only the 3-dimensional case is supported; a higher-dimensional implementation can be produced by trivial adjustments to the present code. For theoretic details see the [companion paper][1].

# Usage

The **RUNME** routine is intended to be the starting point for exploring the code. Upon execution, it will offer to run one of the three main examples: **node_cloud**, **node_earth**, and **node_shell**. The subfolder **./Output** will contain all the output saved to disk, as well as scripts to plot the results with gnuplot, should the user wish to do so (most of the plotting is done within Matlab itself). 
To see details about a specific function/script, type **help _name_of_the_function_** in Matlab prompt while in the 3dRBFnodes folder (it may be necessary to execute *RUNME.m* first to adjust the Matlab PATH).
A number of parameters as well as the distribution density can be set in the preambles of **node_{example}.m** scripts. Executing either one of them adds the *helpers* folder to the Matlab path.
All the helper functions are collected in the **./helpers** subfolder. The ones that are perhaps the easiest to use for generic purposes are:
- repel
- saturate

The routines that have names starting with *paper_* will reproduce the figures in the corresponding section of the [paper][1]. It should be noted that the GPU helpers in **./helpers_gpu** are unstable (and almost certainly unusable) at the time of writing.

# Contributors

Based on the ideas of N. Flyer, B. Fornberg, T. Michaels, and O. Vlasiuk.
See [here](https://github.com/OVlasiuk/3dRBFnodes/graphs/contributors) for the list of individual contributors to the source code.   

# References

The accompanying paper can be found on the [arXiv][1].

For a different approach to 2-dimensional distributions, see [Fast generation of 2-D node distributions for mesh-free PDE discretizations][2].

For how to solve PDEs with RBFs, see [Solving PDEs with radial basis functions][3].

[1]: https://arxiv.org/abs/1710.05011
[2]: https://doi.org/10.1016/j.camwa.2015.01.009
[3]: https://doi.org/10.1017/S0962492914000130
