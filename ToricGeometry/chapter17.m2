-- Author: Simon Telen
-- Date: September 15, 2026
--
-- This file contains code accompanying Chapter 17 in the book
-- "Toric Geometry: Theory and Practice".

-- Example 17.1.12:
-- Construct the GKZ system for A = (0 1 2) in the Weyl algebra and
-- compute its holonomic rank.

needsPackage "Dmodules";

Ahat = matrix {{0,1,2},{1,1,1}};

b = {-4,-3};

n = rank source Ahat;

-- Create the Weyl algebra in variables z_1,...,z_n and differential
-- operators dz_1,...,dz_n.
D = makeWA(QQ[z_1..z_n]);

theta = apply(n, i -> z_(i+1) * dz_(i+1));

-- The binomial operator corresponding to the toric ideal I_Ahat.
I1 = ideal(dz_2^2 - dz_1 * dz_3);

-- The Euler and homogeneity operators.
I2 = ideal(Ahat * (transpose matrix {theta}) - (transpose matrix {b}));

-- The left ideal encoding the GKZ system.
I = I1 + I2;

-- The holonomic rank equals 2 = 1! Vol(A).
holonomicRank I