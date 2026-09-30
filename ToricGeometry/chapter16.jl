# Author: Simon Telen
# Date: September 15, 2026
#
# This file contains code accompanying Chapter 16 in the book
# "Toric Geometry: Theory and Practice".

using LinearAlgebra
using Oscar

# Example 16.2.5:
# Set up the iterative proportional scaling computation for the toric
# statistical model from Examples 16.1.2 and 16.2.2.

A = [1 1 1 0 0 0;
     1 0 0 1 0 0;
     0 1 0 0 1 0;
     1 1 1 1 1 1];

d, n = size(A)

# The weights of the model are all equal to one. In the notation of
# Section 15.3, exp(-w/ε) equals the scaling vector of the model.

w = -log.([1, 1, 1, 1, 1, 1])

# The right-hand side consists of the sufficient statistics A*u/N from
# Example 16.2.2, together with the appended coordinate 1.

b = [59//180, 1//6, 1//3, 1]

ε = 1 # regularization parameter

# The matrix A already has constant column sum c=2, so no preprocessing
# step as in Exercise 15.3.1 is needed.

c = sum(A[:, 1])
x = exp.(-w / ε)

tol = 1e-8
err = Inf
maxiter = 1000
iter = 0

# Apply the iteration (15.3.2). Its limit is the unique positive solution
# \hat p from Proposition 16.2.4.

while err > tol && iter < maxiter
    xold = x

    x = x .* (
        [prod((b ./ (A * x)) .^ (A[:, i])) for i = 1:n]
    ) .^ (1 / c)

    err = norm(x - xold) / norm(x)
    iter += 1
end

x