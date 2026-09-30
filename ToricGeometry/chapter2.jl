# Author: Simon Telen
# Date: September 15, 2026
#
# This file contains code accompanying Chapter 2 in the book
# "Toric Geometry: Theory and Practice".

using Oscar

# Example 2.2.12:
# Construct the cone from Example 2.2.3 and compute its dual.

σ = positive_hull([0 1;
                   1 2;
                   2 1]);

σ_dual = polarize(σ);

# Compute the dimension and facet inequalities of σ.

dim_σ = dim(σ)
facets_σ = facets(σ)

# Check basic properties of the cone.
# The expected output is (true, false, true, true).

is_simplicial(σ), is_smooth(σ), is_fulldimensional(σ), is_pointed(σ)

# Example 2.3.9:
# Compute the Hilbert basis of the dual cone and the toric ideal of the
# corresponding minimal affine embedding of the normal toric variety Y_σ.

H = hilbert_basis(σ_dual)
I = toric_ideal(H)

# Example 2.4.11:
# The cone is the cone over the two-dimensional permutohedron.

σ_perm = positive_hull([1 2 3;
                        2 1 3;
                        1 3 2;
                        3 1 2;
                        2 3 1;
                        3 2 1]);

# Construct the associated normal affine toric variety.
# It has dimension three.

Y_σ_perm = affine_normal_toric_variety(σ_perm)
dim(Y_σ_perm)

# The Hilbert basis of the dual cone has 15 elements. The associated
# affine embedding has a toric ideal with 77 binomial generators.

H_perm = hilbert_basis(polarize(σ_perm))
I_perm = toric_ideal(Y_σ_perm)