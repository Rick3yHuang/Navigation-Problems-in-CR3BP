-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
A1 = (transpose dirA1) * sA1
A2 = (transpose dirA2) * sA2
B1 = (transpose dirB1) * sB1
B2 = (transpose dirB2) * sB2
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
(CPolynomialCoeffMatrixA,CPolynomialCoeffMatrixB) = toSequence CPolynomialCoeffMatrixList;
orbitContraintTuples = {{A1,CPolynomialCoeffMatrixA,cA},{A2,CPolynomialCoeffMatrixA,cA},
			{B1,CPolynomialCoeffMatrixB,cB},{B2,CPolynomialCoeffMatrixB,cB}};
netList orbitContraintTuples
distanceConstraintTuples = {{A1,B1,distAB1},{A2,B2,distAB2}};
netList distanceConstraintTuples
end

-----------------------------------------------------------------------------
---------- start M2 in "../../NavigationProblemsInCR3BP.m2" -------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Lyapunov";		       -- Type of orbit to fit
jacobiConstantDegree = 3;	       -- Degree in Jacobi constant of the model polynomials
modelDegree = 6;		       -- Degree of the model polynomials

X = ZZ/7772777			       -- Finite field for symbolic computation
Y = X[cA,cB,sA1,sA2,sB1,sB2]	       -- Polynomial ring over finite field

-- Random measurements
(dirA1,dirA2,dirB1,dirB2,distAB1,distAB2) = toSequence (apply(4, i -> random(X^1, X^2))|{random(X),random(X)});
CPolynomialCoeffMatrixList = apply(2, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));		 -- Random model coeffs

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Four-Spacecraft-IOD.m2")
elapsedTime findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree})
