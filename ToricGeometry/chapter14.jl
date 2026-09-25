# Author: Simon Telen
# Date: September 15, 2026
#
# This file contains code accompanying Chapter 14 in the book
# "Toric Geometry: Theory and Practice".

using Oscar

# Example 14.2.13:
# Compute the Cox ring and irrelevant ideal of the normal toric surface
# from Example 11.4.15.

P = convex_hull([0 15;
                 0 1;
                 2 0;
                 10 0])

Σ = normal_fan(P)
X_Σ = normal_toric_variety(Σ)

set_coordinate_names(X_Σ, ["y$i" for i = 1:length(rays(Σ))])

S = cox_ring(X_Σ)
B = irrelevant_ideal(X_Σ)

# Compute the matrix F of primitive ray generators used for the grading
# and irrelevant ideal.

F = hcat(rays(Σ)...)
for i = 1:size(F, 2)
    F[:, i] = F[:, i] * lcm(denominator.(F[:, i]))
end
F

# Example 14.4.3:
# Detect boundary solutions of two Laurent equations using homogenization
# in the Cox ring. The parameter z=1 gives a boundary solution, whereas
# z=2 does not.

S, y = polynomial_ring(QQ, :y => 1:7)

z = 1

f1h = y[1] * y[2] * y[3] * y[4] * y[5] * y[6] * y[7] +
      y[3] * y[4]^2 * y[5]^2 * y[6] +
      y[1]^2 * y[2]^2 * y[3] * y[6] * y[7]^2 +
      y[1] * y[5]^2 * y[6]^2 * y[7]^2 +
      y[1] * y[2]^2 * y[3]^2 * y[4]^2

f2h = y[4]^2 * y[5]^2 * y[6] * y[7] +
      y[1] * y[2] * y[4] * y[5] * y[6] * y[7]^2 +
      y[2] * y[3] * y[4]^3 * y[5] +
      y[1] * y[2]^2 * y[3] * y[4]^2 * y[7] +
      z * y[1]^2 * y[2]^2 * y[6] * y[7]^3

Ih = ideal([f1h, f2h])

# The irrelevant ideal of the Cox ring for the normal fan of P1 + P2.

B = ideal([
    y[3] * y[4] * y[5] * y[6] * y[7],
    y[4] * y[5] * y[6] * y[7] * y[1],
    y[5] * y[6] * y[7] * y[1] * y[2],
    y[6] * y[7] * y[1] * y[2] * y[3],
    y[7] * y[1] * y[2] * y[3] * y[4],
    y[1] * y[2] * y[3] * y[4] * y[5],
    y[2] * y[3] * y[4] * y[5] * y[6]
])

# The saturation criterion in Theorem 14.4.2 tests whether all solutions
# lie in the dense torus, equivalently whether there are no boundary solutions.

saturation(Ih + ideal(prod(y)), B) == ideal(S(1)) # false if z = 1