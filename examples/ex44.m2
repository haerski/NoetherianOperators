restart
needsPackage "NumericalAlgebraicGeometry"
needsPackage "Bertini"
loadPackage ("DualSpaces", Reload => true)
R = CC[x_0..x_5]
P = minors(2,matrix{{x_0,x_1,x_3,x_4},{x_1,x_2,x_4,x_5}});
f1 = x_1^4 - 2*x_0*x_1^2*x_2 + x_0^2*x_2^2 + x_1*x_2*x_3*x_4 - x_0*x_2*x_4^2 - x_1^2*x_3*x_5 + x_0*x_1*x_4*x_5
f2 = x_1^4 - 2*x_0*x_1^2*x_2 + x_0^2*x_2^2 + x_1*x_2*x_3*x_4 - x_1^2*x_4^2 - x_0*x_2*x_3*x_5 + x_0*x_1*x_4*x_5
f3 = x_2^2*x_3*x_4 - x_1*x_2*x_4^2 + x_4^4 - x_1*x_2*x_3*x_5 + x_1^2*x_4*x_5 - 2*x_3*x_4^2*x_5 + x_3^2*x_5^2
I = ideal(f1,f2,f3)
nvb = elapsedTime numericalIrreducibleDecomposition (I, Software => BERTINI)
-- choose the component that is P primary
membershipTest = components nvb / (w->clean(1e-6, evaluate(gens P, w#Points#0))) / norm
ws = (components nvb)#(minPosition membershipTest)
-- sample 200 points on this component
pts = bertiniSample(500, ws);
-- Some coefficients will not be properly interpolated, they will output a ?
-- This means that more points are necessary
numericalNoetherianOperators(I,take(pts,2), DependentSet => {x_1,x_3,x_4}, InterpolationTolerance => 1e-10, Saturate => false)

idx := 0;
noethOpsAtPoints = pts / (p -> (<<(idx=idx+1)<<"/"<<#pts<<endl; numNoethOpsAtPoint(I, p, DependentSet => {x_1, x_3, x_4}, DegreeLimit => 5)));
-- remove bad points, i.e. points where the noetherian operators look different than the majority
monLists = noethOpsAtPoints / (i -> i/monomials);
most = commonest tally monLists;
goodIdx = positions(noethOpsAtPoints, i -> (i / monomials) == most#0);
<<"Num good points: " << #goodIdx << " / " << #noethOpsAtPoints << endl;
noethOpsAtPoints = noethOpsAtPoints_goodIdx;
pts = pts_goodIdx;
vals = noethOpsAtPoints /
	(i -> i#2) / coefficients / last / (j -> j#1) @@ flatten @@ entries / (i -> sub(i,CC));
rationalInterpolation(pts, vals, basis(0,4,R), basis(0,4,R), Tolerance => 1e-4)
rationalInterpolation(pts, vals, R, Tolerance => 1e-4)



restart
needsPackage "NumericalAlgebraicGeometry"
needsPackage "Bertini"
needsPackage "SymbolicPowers"
loadPackage ("DualSpaces", Reload => true)
S = QQ[x_0..x_4]
P = minors(2,matrix{{x_0,x_2,x_3},{x_1,x_3,x_4}});
--ones = unique flatten apply(2, i->apply(3,j->(random(0,9),i)))
--M = matrix table(10,2,(i,j) -> if member((i,j), ones) then 1_R else 0)
--I = ideal(gens(P^6) * M)
I = ideal apply(2, i-> sum (P^3_*)_(first random subsets(10,3)))
codim I
IS = ideal(-x_1^3*x_2^3+3*x_0*x_1^2*x_2^2*x_3-3*x_0^2*x_1*x_2*x_3^2+x_0^3*x_3^3-x_1^3*x_3^3-x_1*x_2*x_3^4+x_0*x_3^5+3*x_0*x_1^2*x_3^2*x_4+2*x_1*x_2^2*x_3^2*x_4-2*x_0*x_2*x_3^3*x_4-x_1*x_2^3*x_4^2-3*x_0^2*x_1*x_3*x_4^2+x_0*x_2^2*x_3*x_4^2+x_0^3*x_4^3,-x_1^2*x_2*x_3^3+x_0*x_1*x_3^4-x_1*x_3^5-x_3^6+x_1^2*x_2^2*x_3*x_4-x_0^2*x_3^3*x_4+2*x_1*x_2*x_3^3*x_4+x_0*x_3^4*x_4+3*x_2*x_3^4*x_4-x_0*x_1*x_2^2*x_4^2+x_0^2*x_2*x_3*x_4^2-x_1*x_2^2*x_3*x_4^2-2*x_0*x_2*x_3^2*x_4^2-3*x_2^2*x_3^2*x_4^2+x_0*x_2^2*x_4^3+x_2^3*x_4^3)

--elapsedTime primaryDecomposition I


R = CC monoid S
I = ideal(-x_1^3*x_2^3+3*x_0*x_1^2*x_2^2*x_3-3*x_0^2*x_1*x_2*x_3^2+x_0^3*x_3^3-x_1^3*x_3^3-x_1*x_2*x_3^4+x_0*x_3^5+3*x_0*x_1^2*x_3^2*x_4+2*x_1*x_2^2*x_3^2*x_4-2*x_0*x_2*x_3^3*x_4-x_1*x_2^3*x_4^2-3*x_0^2*x_1*x_3*x_4^2+x_0*x_2^2*x_3*x_4^2+x_0^3*x_4^3,-x_1^2*x_2*x_3^3+x_0*x_1*x_3^4-x_1*x_3^5-x_3^6+x_1^2*x_2^2*x_3*x_4-x_0^2*x_3^3*x_4+2*x_1*x_2*x_3^3*x_4+x_0*x_3^4*x_4+3*x_2*x_3^4*x_4-x_0*x_1*x_2^2*x_4^2+x_0^2*x_2*x_3*x_4^2-x_1*x_2^2*x_3*x_4^2-2*x_0*x_2*x_3^2*x_4^2-3*x_2^2*x_3^2*x_4^2+x_0*x_2^2*x_4^3+x_2^3*x_4^3)


nvb = elapsedTime numericalIrreducibleDecomposition (I, Software => BERTINI)
components nvb
-- choose the component that is P primary
membershipTest = components nvb / (w->clean(1e-6, evaluate(gens P, w#Points#0))) / norm
ws = (components nvb)#(minPosition membershipTest)
-- sample 200 points on this component
pts = bertiniSample(1000, ws);
-- Some coefficients will not be properly interpolated, they will output a ?
-- This means that more points are necessary
elapsedTime numericalNoetherianOperators(I,pts, DependentSet => {x_2,x_3}, InterpolationTolerance => 1e-2)





restart
needsPackage "NumericalAlgebraicGeometry"
needsPackage "Bertini"
needsPackage "SymbolicPowers"
loadPackage ("DualSpaces", Reload => true)
S = QQ[x_0..x_5]
P = minors(2,matrix{{x_0,x_1,x_3,x_4},{x_1,x_2,x_4,x_5}});
--ones = unique flatten apply(2, i->apply(3,j->(random(0,9),i)))
--M = matrix table(10,2,(i,j) -> if member((i,j), ones) then 1_R else 0)
--I = ideal(gens(P^6) * M)
I = ideal apply(3, i-> sum ((P^2_*)_(first random subsets(21,3)) / (i -> i*random(QQ, Height => 1000000))))
codim I
I = P^2 * (ideal apply(3, i->random(1,S)))


--elapsedTime primaryDecomposition I


R = CC monoid S
I = sub(I, R)


nvb = elapsedTime numericalIrreducibleDecomposition (I, Software => BERTINI)
components nvb
-- choose the component that is P primary
membershipTest = components nvb / (w->clean(1e-6, evaluate(gens P, w#Points#0))) / norm
ws = (components nvb)#(minPosition membershipTest)
-- sample 200 points on this component
pts = bertiniSample(20, ws);
-- Some coefficients will not be properly interpolated, they will output a ?
-- This means that more points are necessary
elapsedTime numericalNoetherianOperators(I,pts, DependentSet => {x_1,x_3,x_4}, InterpolationTolerance => 1e-2)