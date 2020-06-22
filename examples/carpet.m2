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

J = (0,3,5,6,9) / (i -> I_i) // ideal
codim J

elapsedTime (minimalPrimes J / (P -> sort noetherianOperators(J,P)))

-- Numerical version
S = CC monoid R
J' = sub(J,S)
elapsedTime nid = bertiniPosDimSolve(J', BertiniInputConfiguration => {RandomSeed => 1})

ws = (components nid)#0
pts = bertiniSample(50, ws, BertiniInputConfiguration => {RandomSeed => 1});
numericalNoetherianOperators(J', pts, DependentSet => {S_"x_1", S_"x_2", S_"x_3", S_"y_2", S_"y_3"}, InterpolationTolerance => 1e-6, Saturate => false, NoetherianDegreeLimit => 2)

ws = (components nid)#1
pts = bertiniSample(350, ws, BertiniInputConfiguration => {RandomSeed => 1});
numericalNoetherianOperators(J', pts, DependentSet => {S_"x_0", S_"x_1", S_"x_2", S_"y_1", S_"y_2"}, InterpolationTolerance => 1e-6, Saturate => false, NoetherianDegreeLimit => 3)

ws = (components nid)#2
pts = bertiniSample(100, ws, BertiniInputConfiguration => {RandomSeed => 1});
numericalNoetherianOperators(J', pts, DependentSet => {S_"x_0", S_"x_1", S_"y_0", S_"y_1", S_"y_2"}, InterpolationTolerance => 1e-6, Saturate => false, NoetherianDegreeLimit => 2)

ws = (components nid)#3
pts = bertiniSample(350, ws, BertiniInputConfiguration => {RandomSeed => 1});
numericalNoetherianOperators(J', pts, DependentSet => {S_"x_1", S_"x_2", S_"y_1", S_"y_2", S_"y_3"}, InterpolationTolerance => 1e-6, Saturate => false, NoetherianDegreeLimit => 3)

ws = (components nid)#4
pts = bertiniSample(120, ws, BertiniInputConfiguration => {RandomSeed => 1});
numericalNoetherianOperators(J', pts, DependentSet => {S_"x_0", S_"x_1", S_"x_2", S_"x_3", S_"y_1"}, InterpolationTolerance => 1e-6, Saturate => false, NoetherianDegreeLimit => 2)

ws = (components nid)#5
pts = bertiniSample(120, ws, BertiniInputConfiguration => {RandomSeed => 1});
numericalNoetherianOperators(J', pts, DependentSet => {S_"x_0", S_"y_0", S_"y_1", S_"y_2", S_"y_3"}, InterpolationTolerance => 1e-6, Saturate => false, NoetherianDegreeLimit => 2)

ws = (components nid)#6
pts = bertiniSample(120, ws, BertiniInputConfiguration => {RandomSeed => 1});
numericalNoetherianOperators(J', pts, DependentSet => {S_"x_0", S_"x_1", S_"x_2", S_"y_0", S_"y_1"}, InterpolationTolerance => 1e-6, Saturate => false, NoetherianDegreeLimit => 2)




-- heuristic to check indep set
foo = subsets(gens S, 5)
apply(foo, j -> #numNoethOpsAtPoint(J', pts#0, DependentSet => j, DegreeLimit => 1))
m = min oo
positions(ooo, i -> i == m)
foo = foo_oo

apply(foo, j -> #numNoethOpsAtPoint(J', pts#0, DependentSet => j, DegreeLimit => 2))
m = min oo
positions(ooo, i -> i == m)
foo = foo_oo

apply(foo, j -> #numNoethOpsAtPoint(J', pts#0, DependentSet => j, DegreeLimit => 3))
m = min oo
positions(ooo, i -> i == m)
foo = foo_oo

apply(foo, j -> #numNoethOpsAtPoint(J', pts#0, DependentSet => j, DegreeLimit => 4))
m = min oo
positions(ooo, i -> i == m)
foo = foo_oo
