-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
A = matrix{{xA},{yA}};
B = matrix{{xB},{yB}};
M = transpose M;
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
CPolynomialCoeffMatrix = first CPolynomialCoeffMatrixList;
orbitContraintTuples = {{A,CPolynomialCoeffMatrix,c},
                        {B,CPolynomialCoeffMatrix,c}};
netList orbitContraintTuples
distanceConstraintTuples = {{A,M,distAM},{B,M,distBM},{A,B,distAB}};
netList distanceConstraintTuples
end

-----------------------------------------------------------------------------
---------- start M2 in "../../NavigationProblemsInCR3BP.m2" -------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Lyapunov";		       -- Type of orbit to fit
jacobiConstantDegree = 3;	       -- Degree of the Jacobi constant polynomial
modelDegree = 4;		       -- Degree of the orbit model polynomials

X = ZZ/7772777
Y = X[c,xA,yA,xB,yB]

-- Random measurements
(M,distAM,distBM,distAB) = (random(X^1, X^2),random(X),random(X),random(X));
CPolynomialCoeffMatrixList = apply(1, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-One-Mothership-Two-Spacecraft-IOD.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree}) -- found 84 complex solutions
