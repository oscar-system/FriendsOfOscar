# Toric Geometry: Theory and Practice

## Authors / Creators

- [Simon Telen](https://sites.google.com/view/simontelen/home-page) 

## Description

This repository contains the Julia and Macaulay2 code accompanying the book:

> **Toric Geometry: Theory and Practice**  
> Simon Telen, September 15, 2026.

A draft of the book is available [here](https://sites.google.com/view/simontelen/teaching/toric-book).

The code is organized by chapter:

- `chapterX.jl` contains the Julia/Oscar code appearing in Chapter `X`.
- `chapterX.m2` contains the Macaulay2 code appearing in Chapter `X`, when applicable.

The code implements or illustrates selected examples, exercises, and algorithms from the text. Comments point to the relevant statements in the book while keeping the original computations unchanged.

## Software

Most Julia files use:

- [Oscar.jl](https://www.oscar-system.org/) (v1.8.2) for commutative algebra, polyhedral geometry, toric varieties, and Gröbner-basis computations.

Chapter9.jl also uses: 

- [HomotopyContinuation.jl](https://www.juliahomotopycontinuation.org/) (v2.22.4) for numerical solution of sparse polynomial systems in Chapter 9.

The Macaulay2 files use several standard packages, including:

- `QuasiDegrees`
- `CorrespondenceScrolls`
- `MixedMultiplicity`
- `Polyhedra`
- `DModules`

## Contents

The repository includes computations involving:

- Smith normal forms, monomial maps, and toric ideals;
- cones, Hilbert bases, normal affine toric varieties, and local multiplicities;
- Ehrhart polynomials, lattice volumes, and degrees of projective toric varieties;
- multiprojective toric varieties, Cayley configurations, multidegrees, and mixed volumes;
- A-discriminants, A-resultants, principal A-determinants, and the Horn uniformization;
- sparse polynomial solving and mixed root counts;
- normal fans, abstract normal toric varieties, divisor class groups, and Cox rings;
- integer programming, entropic regularization, iterative proportional scaling, and toric statistical models;
- GKZ systems, Khovanskii bases, toric degenerations, Bézoutian determinants, and tact invariants.

Some files reproduce code snippets directly from the text; others include supporting helper functions needed to make the chapter-level computations self-contained.
