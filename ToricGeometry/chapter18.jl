# Author: Simon Telen
# Date: September 15, 2026
#
# This file contains code accompanying Chapter 18 in the book
# "Toric Geometry: Theory and Practice".

using Oscar

# Example 18.2.8:
# Implement Algorithm 3, the subduction algorithm. Given a polynomial q
# and generators g, the function either expresses q as a polynomial in g
# or returns "FAIL".

function subduce(q, g, o)
    # The columns of A are the leading exponent vectors of the generators.
    LE = [
        collect(exponents(leading_term(gg; ordering = o)))[1]
        for gg in g
    ]

    A = hcat(LE...)
    d, n = size(A)

    S, x = polynomial_ring(QQ, :x => 1:n)
    f = S(0)

    while degree(ideal(q)) > 0 && q != 0
        LT = leading_term(q; ordering = o)
        rhs = collect(exponents(LT))[1]
        c = collect(coefficients(LT))[1]

        # Lattice points u of the polyhedron below encode monomials
        # in the generators whose leading term can cancel LT.
        P = polyhedron(
            (-identity_matrix(ZZ, n), zeros(Int, n)),
            (A, rhs)
        )

        if dim(P) >= 0
            V = Int.(vertices(P)[1])

            f = f + c * prod(x .^ V)
            q = q - c * prod(g .^ V)
        else
            return "FAIL"
        end

        println(q)
    end

    if q != 0
        return f + collect(coefficients(q))[1]
    else
        return f
    end
end

# Example 18.2.8:
# Subduce g1*g3 - g2^3 for the parametrization from Example 18.1.2.

R, t = polynomial_ring(QQ, :t => 1:2)

g = [
    t[1]^2 * t[2] + 1,
    t[1] * t[2],
    t[1] * t[2]^2 + 1
]

q = g[1] * g[3] - g[2]^3

w = [-1, -1]
o = lex(R)
oW = weight_ordering(-w, o)

f = subduce(q, g, oW)
evaluate(f, g) == q


# Example 18.3.3:
# Verify the Khovanskii basis property for the parametrization of a
# degree-five surface in P^4.

R, t, s = polynomial_ring(QQ, :t => 1:2, :s => 1:1)

g = s[1] .* [
    t[1] * (t[1]^2 + t[2]^2),
    t[2] * (t[1]^2 + t[2]^2),
    t[1],
    t[2],
    1
]

w = -[2; 1; 0]
o = lex(R)
oW = weight_ordering(-w, o)

# The leading exponents define the toric special fiber.
Ahat = hcat([
    leading_exponent(gg; ordering = oW)
    for gg in g
]...)

IA = toric_ideal(transpose(Ahat))

# Subduce the images of binomial generators of the toric ideal.
f = [subduce(evaluate(b, g), g, oW) for b in gens(IA)]

# Verify that the subduction output represents the original polynomials.
[
    evaluate(f[i], g) == evaluate(gens(IA)[i], g)
    for i = 1:length(f)
]


# Example 18.3.4:
# Compute the Ehrhart polynomial of the polytope of a toric degeneration
# of Gr(2,6) in its Plücker embedding.

n = 6

R, t, s = polynomial_ring(QQ, :t => (1:2, 1:n), :s => 1:1)
tmat = matrix_space(R, 2, n)(t)

g = s .* [
    det(tmat[:, [I[1]; I[2]]])
    for I in subsets(n, 2)
]

w = -[1:n; 2:2:2*n; 0]

o = lex(R)
oW = weight_ordering(-w, o)

Ahat = hcat([
    leading_exponent(gg; ordering = oW)
    for gg in g
]...)

P = convex_hull(transpose(Ahat[1:end-1, :]))
ehrhart_polynomial(P)