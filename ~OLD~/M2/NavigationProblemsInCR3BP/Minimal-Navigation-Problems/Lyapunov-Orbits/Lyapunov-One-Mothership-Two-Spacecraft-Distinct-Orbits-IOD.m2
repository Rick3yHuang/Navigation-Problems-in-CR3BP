-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
massRatioEarthMoon = 1.215058560962404e-2;
A = matrix{{xA},{yA}};
B = matrix{{xB},{yB}};
M = transpose M;
vA = matrix{{vxA},{vyA}};
vB = matrix{{vxB},{vyB}};
vM = matrix{{0},{0}};
sqrtDistA = matrix{{r1A},{r2A}};
sqrtDistB = matrix{{r1B},{r2B}};
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
(CPolynomialCoeffMatrixA,CPolynomialCoeffMatrixB) = toSequence CPolynomialCoeffMatrixList;
orbitContraintTuples = {{A,vA,CPolynomialCoeffMatrixA,cA,massRatioEarthMoon,sqrtDistA},
                        {B,vB,CPolynomialCoeffMatrixB,cB,massRatioEarthMoon,sqrtDistB}};
netList orbitContraintTuples
distanceConstraintTuples = {{A,M,vA,vM,distAM,distDerivativeAM},
                            {B,M,vB,vM,distBM,distDerivativeBM},
			    {A,B,vA,vB,distAB,distDerivativeAB}};
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
modelDegree = 4;		       -- Degree of the model polynomials

X = ZZ/7772777			       -- Finite field for symbolic computation
Y = X[cA,cB,xA,yA,vxA,vyA,xB,yB,vxB,vyB,r1A,r2A,r1B,r2B]	       -- Polynomial ring over finite field

-- Random measurements
(M,distA1M,distB1M,distA1B1,distA2M,distB2M,distA2B2) = toSequence ({random(X^1, X^2)}|flatten entries random(X^1,X^6));
CPolynomialCoeffMatrixList = apply(2, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));		 -- Random model coeffs

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-One-Mothership-Two-Spacecraft-Distinct-Orbits-IOD.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree},DerivativesOfDistances => true, OrbitScenario => orbitType)
