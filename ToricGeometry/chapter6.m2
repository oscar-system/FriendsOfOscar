-- Author: Simon Telen
-- Date: September 15, 2026
--
-- This file contains code accompanying Chapter 6 in the book
-- "Toric Geometry: Theory and Practice".

-- Example 6.3.16:
-- Compute the Cayley configuration, its toric ideal, the multigraded
-- Hilbert polynomial, and mixed volumes for the running example.

A1 = matrix {{0,-1,1,0,0},{0,0,0,-1,1}};
A2 = matrix {{0,1,0,1,2},{0,0,1,1,0}};

n = {rank(source A1), rank(source A2)}

-- Form the Cayley configuration Cay(A1,A2).
toprows = transpose(transpose(A1) || transpose(A2));
r1 = join(apply(n_0, i -> 1), apply(n_1, i -> 0));
r2 = join(apply(n_0, i -> 0), apply(n_1, i -> 1));
bottomrows = matrix {r1,r2};
Cay = toprows || bottomrows;

-- Work in the Z^2-graded coordinate ring of P^4 × P^4.
needsPackage "QuasiDegrees";

degs = join(apply(n_0, i -> {1,0}), apply(n_1, i -> {0,1}));
S = QQ[X_1..X_(n_0), Y_1..Y_(n_1), Degrees => degs];

IA = toricIdeal(Cay, S);

-- Compute the multigraded Hilbert polynomial.
needsPackage "CorrespondenceScrolls";

M = cokernel matrix gens IA;
multiHilbertPolynomial(M)

-- Compute the mixed volumes giving the leading coefficients of the
-- multigraded Hilbert polynomial.

needsPackage "MixedMultiplicity";

listA1 = apply(n_0, i -> entries A1_i);
listA2 = apply(n_1, i -> entries A2_i);

(mMixedVolume {listA1, listA1})/(2!*0!)
(mMixedVolume {listA1, listA2})/(1!*1!)
(mMixedVolume {listA2, listA2})/(0!*2!)

-- Compute the Minkowski sum of the two Newton polytopes.

needsPackage "Polyhedra";

P1 = convexHull A1;
P2 = convexHull A2;
minksum = P1 + P2;