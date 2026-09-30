# Author: Simon Telen
# Date: September 15, 2026
#
# This file contains code accompanying Chapter 3 in the book
# "Toric Geometry: Theory and Practice".

using Oscar

# Example 3.1.12:
# Construct the lattice polygon from Figure 3.3 and compute basic
# combinatorial data and its Ehrhart polynomial.

P = convex_hull([0 0;
                 1 0;
                 0 1;
                 2 1;
                 1 2])

dim(P)
length(vertices(P))
length(facets(P))

E_P = ehrhart_polynomial(P)

# Exercise 3.3.11:
# Compute the degree of the projective toric variety X_A in two ways:
# from its toric ideal, and from the normalized lattice volume of Conv(A).

A = [1 0 0 1 2 2 1;
     0 1 2 2 1 0 1]

Ahat = [1 0 0 1 2 2 1;
        0 1 2 2 1 0 1;
        1 1 1 1 1 1 1]

# The homogeneous toric ideal of the affine cone over X_A.
IAhat = toric_ideal(transpose(Ahat))

# Compute the lattice index of the affine lattice Z'_A using the Smith
# normal form of the matrix obtained by translating the first column to zero.

S, P_snf, Q = snf_with_transform(
    matrix_space(ZZ, size(A)...)(A .- A[:, 1])
)

lattice_index = prod(diagonal(S))

# Kushnirenko's theorem: deg(X_A) equals the normalized volume of Conv(A),
# divided by the affine-lattice index.

degree(IAhat)
lattice_volume(convex_hull(transpose(A))) // lattice_index

# Example 3.5.21: 
# Checking very ampleness and normality in Oscar: 

P = convex_hull([0 0; 1 0; 0 1; 2 1; 1 2]);
is_very_ample(P) 
is_normal(P) 