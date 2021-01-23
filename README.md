# NoetherianOperators
This repository contains Macaulay2 code accompanying the paper "Noetherian operators and primary decomposision" by J. Chen, M. Härkönen, R. Krone, and A. Leykin.
Algorithms are implemented in the package `NoetherianOperators` (distributed with Macaulay2 since version 1.17), and examples can be found in the `examples` directory. In particular, the file `examples/paper_examples.m2` will run through every example in the aforementioned paper.

## Quick start guide
Open a Macaulay2 instance in this directory. Load the package by running `needsPackage "NoetherianOperators"`.

The function `noetherianOperators(I,P)` computes a set of Noetherian operators for the (isolated) P-primary component of I. The function `noetherianOperators I` corresponds to `noetherianOperators(I, radical I)`.

The function `noetherianOperators(I,P, Strategy => "Hybrid")` symbolically computes a set of Noetherian operators for the (isolated) P-primary component of I by using numerical data as a "template". We often see performance improvements compared to the purely symbolic `noetherianOperators`.

To compute specialized Noetherian operators at a point p numerically, use `specializedNoetherianOperators(I,p)`. The output is a set of specialized Noetherian operators for the primary ideal corresponding to the irreducible component containing p.

The function `numericalNoetherianOperators(I,pts,DependentSet => {...})` numerically computes Noetherian operators for the primary component of `I` containing the points in `pts`. At the moment, the option `DependentSet` is mandatory. To find a suitable dependent set numerically, one can use the package `NumericalImplicitization`. See the end of the file `examples/carpet.m2` for details.

Rational function interpolation is implemented in `rationalInterpolation(pts,vals, numBasis, denBasis)`, where `pts` is a list of points, `vals` is a list of values of the rational function at each point, and `numBasis`,`denBasis` are the monomial ansatzes for the numerators and denominator respectively.

## Documentation
The full documentation of the package can be viewed by running
```macaulay2
needsPackage "NoetherianOperators"
viewHelp NoetherianOperators
```