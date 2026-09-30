# Author: Simon Telen
# Date: September 15, 2026
#
# This file contains code accompanying Chapter 6 in the book
# "Toric Geometry: Theory and Practice".

using Oscar

# Example 6.2.8:
# Compute the Minkowski sum of the polygons Conv(A1) and Conv(A2).

A1 = [0 -1 1 0 0;
      0  0 0 -1 1];

A2 = [0 1 0 1 2;
      0 0 1 1 0];

P = minkowski_sum(
    convex_hull(transpose(A1)),
    convex_hull(transpose(A2))
)

length(vertices(P)) # returns 7

# Example 6.3.3:
# Compute the multidegree coefficients of the surface X_A in P^4 × P^4
# from Example 6.1.4.

# Cayley configuration and toric ideal.
Cay = [A1 A2;
       ones(Int, 1, 5)  zeros(Int, 1, 5);
       zeros(Int, 1, 5) ones(Int, 1, 5)];

I = toric_ideal(transpose(Cay))

# Transfer the ideal to a polynomial ring with variables x and y,
# corresponding to the two projective factors.

R, x, y = polynomial_ring(QQ, :x => 1:5, :y => 1:5)
phi = hom(base_ring(I), R, [x; y])
I = phi(I)

# For i = 2,1,0, impose i generic linear equations in x and 2-i
# generic linear equations in y. Saturation by the irrelevant ideal
# removes solutions with all x- or all y-coordinates equal to zero.

for i = 2:-1:0
    J = [(rand(-100:100, 1, 5) * x)[1] for j = 1:i]
    J = [J; [(rand(-100:100, 1, 5) * y)[1] for j = 1:2-i]]
    K = saturation(I + ideal(J), intersect(ideal(x), ideal(y)))
    println("delta_($i,$(2-i)): $(degree(K))")
end