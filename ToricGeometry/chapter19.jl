# Author: Simon Telen
# Date: September 15, 2026
#
# This file contains code accompanying Chapter 19 in the book
# "Toric Geometry: Theory and Practice".

using Oscar

# Example 19.1.3:
# Compute the Bézoutian matrix and Chow form of the rational normal
# curve of degree e.

e = 5
pluckinds = subsets(e + 1, 2)

R, p = polynomial_ring(QQ, :p => pluckinds)

function get_entry(l, m)
    inds = findall(
        s -> s[1] <= l && s[2] >= l + 1 && s[1] + s[2] == l + m + 1,
        pluckinds
    )

    return sum(p[inds])
end

M = reshape([get_entry(i, j) for i = 1:e for j = 1:e], e, e)
M = matrix_space(R, size(M)...)(M)

# The determinant is the Chow form in primal Plücker coordinates.
Chow = det(M)


# Section 19.4:
# Compute the tact invariant of two plane conics by eliminating the
# projective coordinates of a tangency point.

R, a, b, x = polynomial_ring(QQ, :a => 1:6, :b => 1:6, :x => 1:3)

f = a[1] * x[1]^2 +
    a[2] * x[1] * x[2] +
    a[3] * x[1] * x[3] +
    a[4] * x[2]^2 +
    a[5] * x[2] * x[3] +
    a[6] * x[3]^2

g = b[1] * x[1]^2 +
    b[2] * x[1] * x[2] +
    b[3] * x[1] * x[3] +
    b[4] * x[2]^2 +
    b[5] * x[2] * x[3] +
    b[6] * x[3]^2

# The ideal imposes f=g=0 and rank at most one for the 2-by-3
# Jacobian matrix of the two conics.

eqs = [f; g; minors(jacobi_matrix([f; g])[end-2:end, :], 2)]

# Saturation by the irrelevant ideal removes the origin in affine
# coordinates before eliminating the point variables.

I = saturation(ideal(eqs), ideal(x));
tact_invariant = eliminate(I, x)