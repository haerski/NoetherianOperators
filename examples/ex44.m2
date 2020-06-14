needsPackage "NumericalAlgebraicGeometry"
needsPackage "Bertini"
needsPackage "DualSpaces"
x = symbol x;
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
-- sample points on this component
pts = bertiniSample(60, ws);
numericalNoetherianOperators(I,pts, DependentSet => {x_1,x_3,x_4})