restart
setRandomSeed 1
needsPackage "NoetherianOperators"
R = QQ[x_0..x_5]
P = minors(2,matrix{{x_0,x_1,x_3,x_4},{x_1,x_2,x_4,x_5}});
f1 = x_1^4 - 2*x_0*x_1^2*x_2 + x_0^2*x_2^2 + x_1*x_2*x_3*x_4 - x_0*x_2*x_4^2 - x_1^2*x_3*x_5 + x_0*x_1*x_4*x_5
f2 = x_1^4 - 2*x_0*x_1^2*x_2 + x_0^2*x_2^2 + x_1*x_2*x_3*x_4 - x_1^2*x_4^2 - x_0*x_2*x_3*x_5 + x_0*x_1*x_4*x_5
f3 = x_2^2*x_3*x_4 - x_1*x_2*x_4^2 + x_4^4 - x_1*x_2*x_3*x_5 + x_1^2*x_4*x_5 - 2*x_3*x_4^2*x_5 + x_3^2*x_5^2
I = ideal(f1,f2,f3)
primes = minimalPrimes I

-- In this particular example, computing the kernel via Gaussian reduction is faster
-- than using Grobner bases
-- This will take about two minutes
elapsedTime nops = primes / (P -> elapsedTime noetherianOperators(I,P, KernelStrategy => "Gaussian"))

netList nops#0
netList nops#1
netList nops#2
netList nops#3
netList nops#4

-- Example 4.4
witness = first components bertiniPosDimSolve(P, BertiniInputConfiguration => {RandomSeed => 1})
p = first bertiniSample(1, witness, BertiniInputConfiguration => {RandomSeed => 1})
-- Computing specialized Noetherian operators is very fast (p must lie on the component of interest)
elapsedTime specializedNoetherianOperators(I, p, DependentSet => gens R - set support first independentSets P)
-- Using the hybrid method, we reduce the computation time to around four seconds.
elapsedTime noetherianOperators(I,P, Strategy => "Hybrid", KernelStrategy => "Gaussian")
-- This computes a primary decomposition using the Hybrid strategy
elapsedTime nops = primes / (P -> elapsedTime noetherianOperators(I,P, Strategy => "Hybrid", KernelStrategy => "Gaussian"))



-- Example 4.8
witness = first components bertiniPosDimSolve(P, BertiniInputConfiguration => {RandomSeed => 2})
-- The coefficients that couldn't be interpolated within the interpolation degree limit
-- and tolerance will appear as question marks
netList numericalNoetherianOperators(I,
	Sampler => (n, I) -> bertiniSample(n,witness),
	InterpolationDegreeLimit => 1,
	DependentSet => gens R - set support first independentSets P)