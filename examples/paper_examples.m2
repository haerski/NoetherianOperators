restart
needsPackage "NoetherianOperators"

-- Example 2.2
R = QQ[x,y]
n = 4
I = ideal((x+y+1)^n)
noetherianOperators I

-- Example 3.17, through our algorithm
x = symbol x;
R = QQ[x_1..x_3]
Q = ideal((x_1^2 - x_3)^2, x_2 - x_3 * (x_1^2 - x_3))
P = radical Q
noetherianOperators Q

-- Example 3.17, through punctual Hilbert scheme
y = symbol y;
F = frac(R/P)
S = F[y_1,y_2]
gamma = map(S,R, {y_1 + x_1, y_2 + x_2, x_3})
I = gamma Q + (ideal(y_1,y_2))^2
zeroDimensionalDual(origin S, I)

-- Example 3.18, see the file ex318.m2

-- Example 4.3
R = QQ[t,x,y]
I = ideal(x^2, y^2 - x*t)
specializedNoetherianOperators(I, point{{1,0,0}}, DependentSet => {x,y})
-- Symbolic version
noetherianOperators(I, DependentSet => {x,y}) // netList


-- Example 4.4, see the file ex318.m2


-- Example 4.5
-- numerical version
R = QQ[t,x,y]
I = ideal(x^2, y^2 - x*t)
pts = toList(1..6) / (i -> matrix{{i_QQ, 0, 0}})
nops = pts / (p -> specializedNoetherianOperators(I, p, DependentSet => {x,y}));
netList nops
-- interpolating e.g. last coefficient
vals = {1, 1/2, 1/3, 1/4, 1/5, 1/6} -- coefficient in _dx*dy
ts = {matrix{{1}},matrix{{2}},matrix{{3}},matrix{{4}},matrix{{5}},matrix{{6}}} -- corresponding value of t
S = RR[t]
numBasis = denBasis = basis(0,1,S)
rationalInterpolation(ts, vals, numBasis, denBasis) -- output is in the form (numerator, denominator)
-- hence we conclude that the coefficient of dx*dy is 1/t.
--Symbolically
noetherianOperators(I, DependentSet => {x,y}) // netList


-- The interpoltation procedure above is automated in numericalNoetherianOperators
-- Due to limitations of Macaulay2, this function does not clear denominators
numericalNoetherianOperators(I, DependentSet => {x,y})

-- Example 4.7 is in the file carpet2.m2
-- Example 4.8 is in the file ex318.m2





-- Example 4.4 is in the file ex44.m2
-- Example 4.5 is in the file ex45.m2