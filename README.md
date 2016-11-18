# Boxed lattices for RBF computations.
**3dRBFnodes**

This is a Matlab code for placing scaled rational nodes in the unit cube with
variable density.  The theoretic details *will be* described in the
documentation. The algorithm processes the domain by dividing it into a number
of boxes and putting nodes with constant density in each of them. 

A number of parameters can be set in the preamble of **par_node_dis.m**, which
contains the main routine.

# Usage

Run *par_node_dis.m*. The subfolder **./Output** contains all the output:
the coordinates of the generated nodes (*cnf.txt*), an illustration (*nodes.fig*), and a histogram containing the distribution of distances to the nearest neighbor(*histogram.png*). Number of bins in the histogram is adjusted in **repel.m**.

# Contributors

Based on the ideas of N. Flyer, B. Fornberg, T. Michaels and O. Vlasiuk.

For a related, although completely different results, see
[Fast generation of 2-D node distributions for mesh-free PDE discretizations.](https://doi.org/10.1016/j.camwa.2015.01.009)
