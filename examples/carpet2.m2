restart
-- Symbolic primary decomposition via Noetherian Operators
needsPackage "K3Carpets"
needsPackage "NoetherianOperators"

setRandomSeed 1

J = carpet(3,3, Characteristic => 0)
#J_*
codim J
R =ring J

I = ideal((gens J) * random(R^10,R^5));
codim I

-- These will not finish
--elapsedTime minimalPrimes I
--elapsedTime primaryDecomposition I

-- Numerical version
S = CC monoid R
I' = sub(I,S)
elapsedTime nid = bertiniPosDimSolve(I', BertiniInputConfiguration => {RandomSeed => 1})

-- We want a sampler that can sample points on a witness set called ws
sampler = (n,I) -> bertiniSample(n, ws, BertiniInputConfiguration => {RandomSeed => 1});

ws = (components nid)#0
elapsedTime numericalNoetherianOperators(I', Sampler => sampler, DependentSet => {1, 2, 4, 5, 6} / (i -> S_i), InterpolationTolerance => 1e-12, Saturate => false, NoetherianDegreeLimit => 2)

ws = (components nid)#1
elapsedTime numericalNoetherianOperators(I', Sampler => sampler, DependentSet => {0, 1, 2, 3, 4} / (i -> S_i), InterpolationTolerance => 1e-6, Saturate => false, NoetherianDegreeLimit => 3)


-- Noetherian operators of J
noetherianOperators(J)