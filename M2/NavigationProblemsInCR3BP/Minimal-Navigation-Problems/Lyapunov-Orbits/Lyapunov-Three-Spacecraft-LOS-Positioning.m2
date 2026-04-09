-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
A = matrix{{xA},{yA}};
B = A - distAB*matrix{{cosThetaAB},{sinThetaAB}};
C = A - distAC*matrix{{cosThetaAC},{sinThetaAC}};
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
(CPolynomialCoeffMatrixA,CPolynomialCoeffMatrixB) = toSequence CPolynomialCoeffMatrixList;
orbitContraintTuples = {{A,CPolynomialCoeffMatrixA,cA},{B,CPolynomialCoeffMatrixB,cB},{C,CPolynomialCoeffMatrixB,cB}};
netList orbitContraintTuples
distanceConstraintTuples = {{A,B,distAB,cosThetaAB,sinThetaAB},{A,C,distAC,cosThetaAC,sinThetaAC}};
netList distanceConstraintTuples
end

-----------------------------------------------------------------------------
------------ start M2 in "../../NavigationProblemsInCR3BP.m2" ---------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Lyapunov";		       -- Type of orbit to fit
jacobiConstantDegree = 3;	       -- Degree in Jacobi constant of the model polynomials
modelDegree = 4;		       -- Degree of the model polynomials

X = ZZ/7772777			       -- Finite field for symbolic computation
Y = X[xA,yA,cB,sinThetaAB,sinThetaAC]	       -- Polynomial ring over finite field

-- Random measurements
(distAB,distAC,distBC,cA) = toSequence flatten entries random(X^1,X^4);
(cosThetaAB, cosThetaAC) = toSequence flatten entries random(X^1,X^2);
CPolynomialCoeffMatrixList = apply(2, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));		 -- Random model coeffs

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Three-Spacecraft-LOS-Positioning.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree})
