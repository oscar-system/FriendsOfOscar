# Author: Simon Telen
# Date: September 15, 2026
#
# This file contains code accompanying Chapter 15 in the book
# "Toric Geometry: Theory and Practice".

using LinearAlgebra
using Oscar

# Example 15.1.4:
# Compute a feasible integer transport plan for the optimal transport problem
# from Exercise 15.0.2, using normal forms modulo the ideal from
# Exercise 1.3.7.

A = [1 1 1 0 0 0;
     0 0 0 1 1 1;
     1 0 0 1 0 0;
     0 1 0 0 1 0;
     0 0 1 0 0 1];

w = [3; 5; 7; 11; 13; 17];
b = [66; 51; 36; 26; 55];

# The variables x encode the transport plan; the variables t encode the
# row and column constraints, together with the localization variable.

R, vrs = polynomial_ring(
    QQ,
    [["x_$i" for i = 1:6]; ["t_$j" for j = 1:6]]
);

x = vrs[1:6];
t = vrs[7:end];

J = ideal([
    prod(t) - 1;
    [x[i] - prod(t[1:end-1] .^ (A[:, i])) for i = 1:6]
]);

NF = normal_form(prod(t[1:end-1] .^ b), J)

# The normal form is the monomial x^u, from which we read the feasible plan.
u = [36; 26; 4; 0; 0; 51];


# Continue Example 15.1.4:
# Compute the cost-minimizing transport plan by taking a normal form
# modulo the toric ideal with respect to the weight order induced by w.

IA = toric_ideal(transpose(A));
S = base_ring(IA);
x = gens(S);

# The weight ordering is refined by lexicographic order.
wo = weight_ordering(w, lex(S));

with_ordering(S, wo) do
    normal_form(prod(x .^ u), IA)
end


# Example 15.3.5:
# Iterative proportional scaling for the entropically regularized optimal
# transport problem from Exercise 15.0.2.

# We remove one redundant row from the transport matrix and its right-hand
# side. The preprocessing below restores a row so that A has constant
# column sum, as required for the iteration in (15.3.2).

A = [1 1 1 0 0 0;
     0 0 0 1 1 1;
     1 0 0 1 0 0;
     0 1 0 0 1 0];

d, n = size(A);

w = [3, 5, 7, 11, 13, 17];
b = [66, 51, 36, 26];

ε = 1; # regularization parameter

# Preprocessing from Exercise 15.3.1: append a row if necessary so that
# all columns of A have the same sum.

vt = ones(Int, 1, d) * A;
newrow = maximum(vt) * ones(Int, 1, n) - vt;

if norm(newrow) != 0
    Anew = [A; newrow];

    A_oscar = matrix_space(QQ, d, n)(A);
    newrow_oscar = matrix_space(QQ, 1, n)(newrow);

    bnew = [b; Rational{Int64}((solve(A_oscar, newrow_oscar) * b)[1])];

    A = Anew;
    b = bnew;
end

c = sum(A[:, 1]);      # constant column sum
x = exp.(-w / ε);      # initial point x^(0)
tol = 1e-8;            # stopping tolerance
err = Inf;
maxiter = 1000;
iter = 0;

# Iteration (15.3.2).

while err > tol && iter < maxiter
    xold = x;

    x = x .* (
        [prod((b ./ (A * x)) .^ (A[:, i])) for i = 1:n]
    ) .^ (1 / c);

    err = norm(x - xold) / norm(x);
    iter += 1;

    if mod(iter, 10) == 5
        println("iter = $iter, error = $err, coordinate sum = $(sum(x))")
    end
end