# Author: Simon Telen
# Date: September 15, 2026
#
# This file contains code accompanying Chapter 13 in the book
# "Toric Geometry: Theory and Practice".

using Oscar

# Example 13.2.10:
# Compute the divisor class group of the normal toric surface from
# Example 11.4.15.

P = convex_hull([0 15;
                 0 1;
                 2 0;
                 10 0])

Σ = normal_fan(P)
X_Σ = normal_toric_variety(Σ)

# The class group is computed from the ray matrix, as described in
# Theorem 13.2.4 and Proposition 13.2.5.

class_group(X_Σ)
