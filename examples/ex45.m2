restart
-- Symbolic primary decomposition via Noetherian Operators
needsPackage "NumericalAlgebraicGeometry"
needsPackage "Bertini"
needsPackage "K3Carpets"
needsPackage "MinimalPrimes"
needsPackage "DualSpaces"

installMinprimes()

setRandomSeed 1

I = carpet(3,3, Characteristic => 0)
#I_*
codim I
R =ring I

J = ideal((gens I) * random(R^10,R^5));
codim J

-- These will not finish
--elapsedTime minimalPrimes J
--elapsedTime primaryDecomposition J

-- Numerical version
S = CC monoid R
J' = sub(J,S)
elapsedTime nid = bertiniPosDimSolve(J', BertiniInputConfiguration => {RandomSeed => 1})




ws = (components nid)#0
pts = bertiniSample(100, ws, BertiniInputConfiguration => {RandomSeed => 1});
elapsedTime numericalNoetherianOperators(J', pts, DependentSet => {1, 2, 4, 5, 6} / (i -> S_i), InterpolationTolerance => 1e-12, NoetherianDegreeLimit => 2)

ws = (components nid)#1
pts = bertiniSample(20, ws, BertiniInputConfiguration => {RandomSeed => 1});
elapsedTime numericalNoetherianOperators(J', pts, DependentSet => {0, 1, 2, 3, 4} / (i -> S_i), InterpolationTolerance => 1e-6, NoetherianDegreeLimit => 3)