# Rational plane tropical curves

This code is based on the algorithm described in Section 7.2. of the paper "Enumerative tropical Algebraic Geometry in R2" by Mikhalkin.

It generates all rational plane tropical curves of a fixed degree d with 3d-1 points, computes their multiplicities, and produces the corresponding Newton subdivisions. Note that this code counts and constructs reducible tropical curves as well. 

For a detailed example with explanations, see Example_Jupyter_Notebook.

The source code is contained in RationalPlaneTropicalCurves.jl
The file Example.jl contains the computation for d=3.

This code was written by Alejandro Ovalle, supervised by Nathan Pflueger and Alheydis Geiger during an internship at Max-Planck-Institute for Mathematics in the Sciences.
