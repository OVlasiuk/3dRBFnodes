# Boxed lattices for RBF computations
<p align="center">
<img src="https://raw.githubusercontent.com/OVlasiuk/3dRBFnodes/master/nodes.png" width="400">
</p>

---
This is a collection of Matlab routines for distributing nodes drawn from an irrational lattice with variable density. The principal application for such nodes are RBF-based PDE solvers. The algorithm processes the working domain by dividing it into a number of equal boxes and putting nodes with constant density in each of them, then applying a repel procedure using the hypersingular Riesz kernel. Currently only the 3-dimensional case is supported; it isn't difficult to make a higher-dimensional implementation along the same lines. Theoretic details will be described in a companion paper.

To see details about a specific function/script, type **help _name_of_the_function_** in Matlab prompt while in the 3dRBFnodes folder (it may be necessary to execute *RUNME.m* first to adjust the Matlab PATH).
A number of parameters as well as the distribution density can be set in the preamble of **node_dis.m** and **node_earth.m**. Executing either one of them adds the *helpers* folder to the Matlab path.

# Usage

Use **RUNME** to run **node_dis.m** or **node_earth.m**. The subfolder **./Output** will contain all the output saved to disk, as well as scripts to plot the results with gnuplot. Number of bins in the histogram is adjusted in **repel.m**.

# Contributors

Based on the ideas of N. Flyer, B. Fornberg, T. Michaels and O. Vlasiuk.

For a different approach to 2-dimensional distributions, see
[Fast generation of 2-D node distributions for mesh-free PDE discretizations.](https://doi.org/10.1016/j.camwa.2015.01.009)

For how to solve PDEs with RBFs, see [Solving PDEs with radial basis functions.](https://doi.org/10.1017/S0962492914000130)
