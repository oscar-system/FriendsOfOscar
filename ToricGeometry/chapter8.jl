# Author: Simon Telen
# Date: September 15, 2026
#
# This file contains code accompanying Chapter 8 in the book
# "Toric Geometry: Theory and Practice".

using Oscar

# The following functions are from Chapter 4. They compute the multiplicity
# of an affine or projective toric variety along a torus orbit.

function get_lattice_basis(A)
    S, P, Q = snf_with_transform(A)
    return inv(P) * S[:, 1:rank(A)]
end

function sublattice_in_linspace(A, L)
    N = nullspace(transpose(L))[2]
    V = nullspace(transpose(N) * A)[2]
    return get_lattice_basis(A * V)
end

function get_lattice_index(A)
    return prod(diagonal(snf(A))[1:size(A, 1)])
end

function get_face_lattice_index(τ, A)
    inds = findall(i -> A[:, i] in τ, 1:size(A, 2))
    L = A[:, inds]
    Rτ = sublattice_in_linspace(A, L)
    newA = solve(Rτ, A[:, inds]; side = :right)
    return get_lattice_index(newA)
end

function SDV(A)
    A = A[:,findall(i->sum(abs.(A[:,i]))!=0, 1:size(A,2))]
    d, n = size(A)
    Awith0 = zero(matrix_space(ZZ, d, n + 1))
    Awith0[:, 2:end] = A
    V1 = volume(convex_hull(transpose(Awith0)))
    V2 = volume(convex_hull(transpose(A)))
    vol = V1 - V2
    return factorial(d) * vol / get_lattice_index(A)
end

function get_multiplicity(τ, A)
    A = matrix_space(ZZ, size(A)...)(A)
    A = hcat(matrix_space(ZZ, size(A, 1), 1)(zeros(Int, size(A, 1))), A)
    inds = findall(i -> A[:, i] in τ, 1:size(A, 2))
    B = nullspace(transpose(A[:, inds]))[2]
    Amodτ = transpose(B) * A
    return get_face_lattice_index(τ, A) * SDV(Amodτ)
end

function get_multiplicity_proj(A, Q)
    if dim(Q) == 0
        At = matrix_space(ZZ, size(A)...)(A .- Array(lattice_points(Q)[1]))
        return get_multiplicity([zeros(Int, size(A, 1))], At)
    else
        Qinds = findall(p -> A[:, p] in Q, 1:size(A, 2))
        vtcs = vertices(convex_hull(transpose(A)))

        # Choose a vertex of Conv(A) contained in Q and work in its affine chart.
        vtx = vtcs[findfirst(v -> v in Q, vtcs)]
        At = matrix_space(ZZ, size(A)...)(A .- Array(vtx))
        Qt = A[:, Qinds] .- Array(vtx)
        τQ = positive_hull([Qt[:, j] for j = 1:size(Qt, 2)])
        return get_multiplicity(τQ, At)
    end
end

# Auxiliary function used for resultant and discriminant computations.
function get_aux_variables(A)
    d, n = size(A)
    Ahat = [A; ones(Int, 1, n)]
    Ahat = matrix_space(ZZ, size(Ahat)...)(Ahat)

    a = [A[:, i] for i = 1:n]
    aplus = [[maximum([aa, 0]) for aa in a[j]] for j = 1:n]
    aminus = [[minimum([aa, 0]) for aa in a[j]] for j = 1:n]

    return d, n, Ahat, aplus, aminus
end

# Proposition 8.1.4:
# Eliminate the torus variables from d+1 Laurent polynomials with common
# support A to compute the defining ideal of the A-resultant variety.

function get_A_resultant(A)
    d, n, _, aplus, aminus = get_aux_variables(A)

    R, t, v, z = polynomial_ring(
        QQ,
        :t => 1:d,
        :v => 1:d,
        :z => (0:d, 1:n)
    )

    eqs1 = matrix(z) * [
        prod(t .^ aplus[i]) * prod(v .^ (-aminus[i]))
        for i = 1:n
    ]

    eqs2 = [v[i] * t[i] - 1 for i = 1:d]

    E = eliminate(ideal([eqs1; eqs2]), [t; v])
end

# Section 8.2:
# Compute an A-discriminant, optionally expressing its defining ideal in
# a prescribed polynomial ring. This is used for discriminants of faces.

function get_A_discriminant(A; vrs = [])
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

    if !isempty(vrs)
        S = parent(vrs[1])
        phi = hom(R, S, [S.(ones(Int, 2*d)); vrs])
        return phi(E)
    else
        return E
    end
end

# Theorem 8.2.3:
# Return the face discriminants and their multiplicities in the principal
# A-determinant. A face with non-hypersurface discriminant contributes 1.

function get_principal_A_det(A)
    P = convex_hull(transpose(A))
    d, n = size(A)

    factorlist = []
    R, z = graded_polynomial_ring(QQ, :z => 1:n)

    if dim(P) != d
        println("P is not full dimensional")
        return []
    end

    for i = 0:dim(P)
        for Q in faces(P, i)
            mult = get_multiplicity_proj(A, Q)

            Qinds = findall(i -> A[:, i] in Q, 1:n)
            AQ = A[:, Qinds]

            Adisc = get_A_discriminant(Array(AQ); vrs = z[Qinds])

            if length(gens(Adisc)) > 1
                Adisc = ideal(R(1))
            end

            push!(factorlist, (Adisc, mult))
        end
    end

    return factorlist
end

# Example usage: 

A = [0 1 0 2 1 0; 0 0 1 0 1 2]
get_A_resultant(A)
get_principal_A_det(A)