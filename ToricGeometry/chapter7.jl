# Author: Simon Telen
# Date: September 15, 2026
#
# This file contains code accompanying Chapter 7 in the book
# "Toric Geometry: Theory and Practice".

using Oscar

# Section 7.3:
# Auxiliary data for elimination and Horn-uniformization computations.

function get_aux_variables(A)
    d, n = size(A)
    Ahat = [A; ones(Int, 1, n)]
    Ahat = matrix_space(ZZ, size(Ahat)...)(Ahat)

    a = [A[:, i] for i = 1:n]
    aplus = [[maximum([aa, 0]) for aa in a[j]] for j = 1:n]
    aminus = [[minimum([aa, 0]) for aa in a[j]] for j = 1:n]

    return d, n, Ahat, aplus, aminus
end

# Proposition 7.3.4:
# Compute the elimination ideal defining the A-discriminant variety.

function get_A_discriminant(A)
    d, n, Ahat, aplus, aminus = get_aux_variables(A)

    R, t, v, z = polynomial_ring(
        QQ,
        :t => 1:d,
        :v => 1:d,
        :z => 1:n
    )

    D = diagonal_matrix([
        prod(t .^ aplus[i]) * prod(v .^ (-aminus[i]))
        for i = 1:n
    ])

    eqs1 = Ahat * D * z
    eqs2 = [v[i] * t[i] - 1 for i = 1:d]

    E = eliminate(ideal([eqs1; eqs2]), [t; v])
end

# Example 7.3.9:
# The discriminant variety of the Segre threefold P^1 x P^2 is the
# rank-one locus of a 2-by-3 matrix.

A = [1 1 1 0 0 0;
     0 0 0 1 1 1;
     1 0 0 1 0 0;
     0 1 0 0 1 0;
     0 0 1 0 0 1];

get_A_discriminant(A)

# Theorem 7.3.12:
# Compute the A-discriminant using the Horn uniformization.

function get_A_discriminant_via_Horn(A)
    d, n, Ahat, aplus, aminus = get_aux_variables(A)

    B = nullspace(Ahat)[2]
    m = size(B, 2)

    R, t, v, u, z = polynomial_ring(QQ,:t => 1:d,:v => 1:d,:u => 1:m,:z => 1:n)

    Bu = B * u

    eqs1 = [
        z[i] * prod(t .^ aplus[i]) -
        prod(t .^ (-aminus[i])) * Bu[i]
        for i = 1:n
    ]

    eqs2 = [v[i] * t[i] - 1 for i = 1:d]

    E = eliminate(ideal([eqs1; eqs2]), [t; v; u])
end

# Exercise 7.3.10:
# Compute the discriminant and its Newton polytope.

A = [4 2 4 -3 -2;
     2 -1 -1 0 -2];

ΔA = get_A_discriminant(A)
P = newton_polytope(gens(ΔA)[1])

dim(P), length(vertices(P))

# Exercise 7.3.13: 
# Compute the same A-discriminant via the Horn uniformization: 
ΔA = get_A_discriminant_via_Horn(A)
