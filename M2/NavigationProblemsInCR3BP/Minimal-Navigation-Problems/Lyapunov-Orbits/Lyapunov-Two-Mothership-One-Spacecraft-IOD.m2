-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
A = matrix{{xA},{yA}};
M1 = transpose M1;
M2 = transpose M2;
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
CPolynomialCoeffMatrix = first CPolynomialCoeffMatrixList;
orbitContraintTuples = {{A,CPolynomialCoeffMatrix,c}}
netList orbitContraintTuples
distanceConstraintTuples = {{A,M1,distAM1},{A,M2,distAM2}};
netList distanceConstraintTuples
end

-----------------------------------------------------------------------------
------------ start M2 in "../../NavigationProblemsInCR3BP.m2" ---------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Lyapunov";		       -- Type of orbit to fit
jacobiConstantDegree = 3;	       -- Degree of the Jacobi constant polynomial
modelDegree = 6;		       -- Degree of the orbit model polynomials

X = ZZ/7772777
Y = X[c,xA,yA]

-- Random measurements
(M1,M2,distAM1,distAM2) = (random(X^1, X^2),random(X^1,X^2),random(X),random(X));
CPolynomialCoeffMatrixList = apply(1, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Two-Mothership-One-Spacecraft-IOD.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree}) -- found 84 complex solutions
