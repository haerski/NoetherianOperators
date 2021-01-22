-- Example computation from the paper “Primary Ideals and Their
-- Differential Equations” by Yairon Cid-Ruiz, Roser Homs and Bernd Sturmfels
----------------------------------------------------
restart
needsPackage "NoetherianOperators"
----------------------------------------------------
-- Example 1
---------------------------------------------------
U= QQ[x_1,x_2,x_3,x_4,u_1,u_2,u_3,u_4,y_1,y_2]
A = matrix {{u_3,u_1,u_2},{u_1,u_2,u_4}}
PP = minors(2,A)
JJ=ideal{PP,x_1-u_1-y_1,x_2-u_2-y_2,x_3-u_3,x_4-u_4,y_1^3,y_2+x_2*y_1^2}
J=ideal{eliminate(JJ,{u_1,u_2,u_3,u_4,y_1,y_2})}

R=QQ[x_1,x_2,x_3,x_4]
F=map(R,U)
Q=F(J)
noetherianOperators(Q, Strategy => "MacaulayMatrix")

----------------------------------------------------
----------------------------------------------------
-- Example 2
---------------------------------------------------

R=QQ[x_1,x_2,x_3,x_4]
Q=ideal{x_1^2,x_1*x_2,x_1*x_3,x_1*x_4-x_3^2+x_1,x_3^2*x_4-x_2^2,x_3^2*x_4-x_3^2-x_2*x_3+2*x_1}
noetherianOperators(Q, Strategy => "MacaulayMatrix")

------------------------------------------

R=QQ[x_1,x_2,x_3,x_4]
MM = matrix {{x_3,x_1,x_2},{x_1,x_2,x_4}}
P = minors(2,MM)
M=ideal{x_1^2,x_2^2,x_3^2,x_4^2}
Q=joinIdeals(P,M)
noetherianOperators(Q, Strategy => "MacaulayMatrix")

----------------------------------------------------
----------------------------------------------------
-- Example 3: Palamodov's example
---------------------------------------------------
R = QQ[x_1, x_2, x_3]
Q = ideal(x_1^2, x_2^2, x_1-x_2*x_3)
noetherianOperators(Q, Strategy => "MacaulayMatrix")
----------------------------------------------------
----------------------------------------------------

----------------------------------------------------
----------------------------------------------------
-- Example 4: taken from page 143 of "Solving Systems of Polynomial Equations"
---------------------------------------------------
R = QQ[x_1, x_2, x_3, x_4]
Q = ideal(x_1^3*x_4^2-x_2^5, x_1^2*x_4^3-x_3^5, x_1*x_3^2-x_2^3, x_2^2*x_4 - x_3^3)
Q1 = ideal(x_1*x_4-x_2*x_3, x_1*x_3^2-x_2^3, x_2^2*x_4-x_3^3)
Q2 = ideal(x_1^2, x_2^2, x_3^2)
Q3 = ideal(x_2^2, x_3^2, x_4^2)
Q4 = ideal(x_1^3, x_2^3, x_3^3, x_4^3, x_1*x_3^2, x_2^2*x_4)
assert(Q == intersect(Q1, Q2, Q3, Q4)) -- check that we copied correctly

---- the Noetherian operators of Q1
isPrime Q1
-- since it is prime we can choose 1 as the Noetherian operator

---- the Noetherian operators of Q2
isPrime Q2
P2 = radical Q2 -- it is equal to (x_1, x_2, x_3)
noetherianOperators(Q2, Strategy => "MacaulayMatrix")

---- the Noetherian operators of Q3
isPrime Q3
P3 = radical Q3 -- it is equal to (x2, x3, x4)
F = map(R, R, {x_4, x_1, x_2, x_3}) -- change of variables to get Noether normal position
Q3' = F Q3
noetherianOperators(Q3', Strategy => "MacaulayMatrix")

---- the Noetherian operators of Q4
isPrime Q4
P4 = radical Q4 -- it is equal to (x1, x2, x3, x4)
-- no need to "take-out" anyone
noetherianOperators(Q4, Strategy => "MacaulayMatrix")
----------------------------------------------------
----------------------------------------------------

----------------------------------------------------
----------------------------------------------------
-- Example 5: some random primary ideal
---------------------------------------------------
R = QQ[x_1,x_2,x_3]
Q = ideal(random(3, R), random(2, R), random(2, R), random(4, R))
assert(dim Q == 0)
noetherianOperators(Q, Strategy => "MacaulayMatrix")
----------------------------------------------------
----------------------------------------------------

----------------------------------------------------
----------------------------------------------------
-- Example 6 : a small example
---------------------------------------------------
R = QQ[x_1,x_2,x_3]
Q = ideal(x_1^2, x_2^2, x_3^2, x_1*x_2 + x_1*x_3 +x_2*x_3)
noetherianOperators(Q, Strategy => "MacaulayMatrix")
----------------------------------------------------
----------------------------------------------------

----------------------------------------------------
----------------------------------------------------
-- Example 7:
---------------------------------------------------
R = QQ[x_1,x_2,x_3,x_4]
J = ideal(x_1^4 + x_2*x_3*x_4, x_2^4 + x_1*x_3*x_4, x_3^4 + x_1*x_2*x_4)
dim J
primDec = primaryDecomposition J
-- here we will only take care of the first primary component...
Q = primDec_0
noetherianOperators(Q, Strategy => "MacaulayMatrix")
-- Alternatively, we can compute these without doing a primary decomposition
noetherianOperators(J, (minimalPrimes J)#0, Strategy => "MacaulayMatrix")
----------------------------------------------------
----------------------------------------------------

----------------------------------------------------
----------------------------------------------------
-- Example 8: powers of the maximal irrelevant ideal
---------------------------------------------------
R = QQ[x_1,x_2,x_3]
mm= ideal vars R
n=4
Q=mm^n
noetherianOperators(Q, Strategy => "MacaulayMatrix")
----------------------------------------------------
----------------------------------------------------

----------------------------------------------------
----------------------------------------------------
-- Example 9:
---------------------------------------------------
R = QQ[x_1,x_2,x_3]
Q = ideal(x_1^2,x_2^2,x_3^2)
noetherianOperators(Q, Strategy => "MacaulayMatrix")
----------------------------------------------------
----------------------------------------------------

