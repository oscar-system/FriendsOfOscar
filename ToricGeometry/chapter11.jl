# Author: Simon Telen
# Date: September 15, 2026
#
# This file contains code accompanying Chapter 11 in the book
# "Toric Geometry: Theory and Practice".

using Oscar

# Example 11.3.10:
# Compute the normal fan of the pentagon from Example 3.1.12.

P = convex_hull([0 0;
                 1 0;
                 0 1;
                 2 1;
                 1 2])

Σ = normal_fan(P)

# Inspect the maximal cones and rays of the normal fan.
cones(Σ, 2)
rays(Σ)

# Example 11.4.15:
# Construct the normal toric variety associated with the normal fan of
# the polygon with vertices (0,15), (0,1), (2,0), and (10,0).

P = convex_hull([0 15;
                 0 1;
                 2 0;
                 10 0])

Σ = normal_fan(P)
X_Σ = normal_toric_variety(Σ)

# The variety is covered by four normal affine toric surfaces.
# Only one of these affine charts is smooth.

cover = affine_open_covering(X_Σ)
[issmooth(Y) for Y in cover]

# The affine charts correspond to the maximal cones of Σ. For each such
# cone, compute the Hilbert basis of its dual cone. The Hilbert basis
# gives a minimal monomial embedding of the corresponding normal affine
# toric surface.

maximal_cones = cones(Σ, 2)

dual_cones = [polarize(σ) for σ in maximal_cones]
hilbert_bases = [hilbert_basis(σ_dual) for σ_dual in dual_cones]

# The first two maximal cones are σ_12 and σ_23 in Example 11.4.15.
# Their Hilbert bases yield the embeddings
#
# Y_{σ_12} ≃ {xz - y^2 = 0} ⊂ C^3,
# Y_{σ_23} ≃ {uw - v^2 = 0} ⊂ C^3.

A_σ12 = hilbert_bases[1]
A_σ23 = hilbert_bases[2]

I_σ12 = toric_ideal(A_σ12)
I_σ23 = toric_ideal(A_σ23)