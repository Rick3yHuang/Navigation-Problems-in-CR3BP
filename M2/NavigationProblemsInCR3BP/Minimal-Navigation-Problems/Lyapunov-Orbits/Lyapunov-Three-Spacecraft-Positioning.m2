-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
A = matrix{{xA},{yA}};
B = matrix{{xB},{yB}};
C = matrix{{xC},{yC}};
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
(CPolynomialCoeffMatrixA,CPolynomialCoeffMatrixB,CPolynomialCoeffMatrixC) = toSequence CPolynomialCoeffMatrixList;
orbitContraintTuples = {{A,CPolynomialCoeffMatrixA,scaledCATrue},
                        {B,CPolynomialCoeffMatrixB,scaledCBTrue},
			{C,CPolynomialCoeffMatrixC,scaledCCTrue}};
netList orbitContraintTuples
distanceConstraintTuples = {{A,B,distAB},{A,C,distAC},{B,C,distBC}};
netList distanceConstraintTuples
end

-----------------------------------------------------------------------------
------------- start M2 in "../../NavigationProblemsInCR3BP.m2" --------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Lyapunov";					  -- Type of orbit to fit
jacobiConstantDegree = 3;	       -- Degree in Jacobi constant of the model polynomials
modelDegree = 4;		       -- Degree of the model polynomials

X = ZZ/7772777
Y = X[xA,yA,xB,yB,xC,yC]

-- Random measurements
(distAB,distAC,distBC,scaledCATrue,scaledCBTrue,scaledCCTrue) = toSequence flatten entries random(X^1,X^6);
CPolynomialCoeffMatrixList = apply(3, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Three-Spacecraft-Positioning.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree}) -- found 256 complex solutions
