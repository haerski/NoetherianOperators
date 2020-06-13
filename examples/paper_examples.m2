restart
-- Example 2.2
needsPackage "DualSpaces"
R = QQ[x,y]
n = 7
I = ideal((x+y+1)^n)
noetherianOperators I
noetherianOperators(I, DependentSet => {R_1})