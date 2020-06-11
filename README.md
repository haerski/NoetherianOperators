# NoetherianOperators
This repository contains Macaulay2 code for computing Noetherian operators symbolically and numerically.
Algorithms are implemented in the package `DualSpaces`, and examples can be found in the `Examples` directory.

## Quick start guide
Open a Macaulay2 instance in this directory. Load the package by running `needsPackage "DualSpaces"`.

The function `noetherianOperators(I,P)` computes a set of Noetherian operators for the (isolated) P-primary component of I. The function `noetherianOperators I` corresponds to `noetherianOperators(I, radical I)`.

To compute specialized Noetherian operators at a point p numerically, use `numNoethOpsAtPoint(I,p)`. The output is a set of specialized Noetherian operators for the primary ideal corresponding to the irreducible component containing p.

Rational function interpolation is implemented in `rationalInterpolation(pts,vals, numBasis, denBasis)`, where `pts` is a list of points, `vals` is a list of values of the rational function at each point, and `numBasis`,`denBasis` are the monomial ansatzes for the numerators and denominator respectively.
