# Author: Simon Telen
# Date: September 15, 2026
#
# This file contains code accompanying Chapter 9 in the book
# "Toric Geometry: Theory and Practice".

using Oscar
using HomotopyContinuation

# Example 9.1.13:
# Solve a generic unmixed system of Laurent polynomials with support A.
# The number of torus solutions equals the normalized volume of Conv(A).

A = [-3 -3 -2 -2 2 3;
      0 -1  0 -1 0 0];

d, n = size(A)

# Declare the torus variables t[1], ..., t[d].
@var t[1:d]

# Construct the Laurent monomials specified by the columns of A.
mons = [prod(t .^ A[:, i]) for i = 1:n]

# Choose a generic real coefficient matrix and form the d Laurent polynomials.
z = randn(2, 6)
f = z * mons

# HomotopyContinuation.jl accepts ordinary polynomials, so clear the
# negative exponents in the Laurent system.
f = prod([t[i]^(-minimum(A[i, :])) for i = 1:d]) * f

# Solve in the torus. The expected number of solutions is 7.
R = HomotopyContinuation.solve(f, only_non_zero = true)
length(solutions(R)) # output: 7


# Example 9.2.8:
# Compute the intersection of the multiprojective toric surface X_A with
# two hyperplanes in P^4 × P^4. The degree is the mixed volume MV(P1,P2)=5.

A1 = [0 -1 1 0 0;
      0  0 0 -1 1];

A2 = [0 1 0 1 2;
      0 0 1 1 0];

# Form the Cayley configuration and its toric ideal.
Cay = [A1 A2;
       ones(Int, 1, 5)  zeros(Int, 1, 5);
       zeros(Int, 1, 5) ones(Int, 1, 5)];

I = toric_ideal(transpose(Cay))

# Move the ideal to a polynomial ring whose variables are grouped into
# coordinates x and y for the two factors P^4 × P^4.
R, x, y = polynomial_ring(QQ, :x => 1:5, :y => 1:5)
phi = hom(base_ring(I), R, [x; y])
I = phi(I)

# The choice z_ij = 1 gives the two hyperplane equations sum(x)=0 and sum(y)=0.
J = ideal([sum(x); sum(y)])

# Saturation by the irrelevant ideal removes components on which all
# x-coordinates or all y-coordinates vanish.
K = saturation(I + J, intersect(ideal(x), ideal(y)))

# The output is (true, 2, 5): the ideal is radical, its affine multicone
# has dimension 2, and the multiprojective intersection consists of
# five points.
is_radical(K), dim(K), degree(K)