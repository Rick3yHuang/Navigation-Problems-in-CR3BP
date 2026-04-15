-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
powerBasisList = apply({cA,cB}, x -> transpose matrix{apply(jacobiConstantDegree+1, i -> x^i)})
heightOrbitCoefficientsList = apply(#heightCPolynomialCoeffMatrixList, i -> (heightCPolynomialCoeffMatrixList_i)*(powerBasisList_i))
uvList = {{uA1,vA1},{uA2,vA2},{uB1,vB1},{uB2,vB2}};
wList = apply(#uvList, i -> (
	(u,v) = toSequence uvList_i;
	modelBasis = returnModelBasis({u,v},modelDegree);
	first flatten entries((matrix{{1}}|matrix{modelBasis})*heightOrbitCoefficientsList_(floor(i/2)))
	))
A1 = matrix{{uA1},{vA1},{wList_0}};
A2 = matrix{{uA2},{vA2},{wList_1}};
B1 = matrix{{uB1},{vB1},{wList_2}};
B2 = matrix{{uB2},{vB2},{wList_3}};
M = transpose M;
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
(CPolynomialCoeffMatrixA,CPolynomialCoeffMatrixB) = toSequence CPolynomialCoeffMatrixList;
(heightCPolynomialCoeffMatrixA,heightCPolynomialCoeffMatrixB) = toSequence heightCPolynomialCoeffMatrixList;
orbitContraintTuples = {{A1,CPolynomialCoeffMatrixA,cA},{A2,CPolynomialCoeffMatrixA,cA},
                        {B1,CPolynomialCoeffMatrixB,cB},{B2,CPolynomialCoeffMatrixB,cB}};
netList orbitContraintTuples
distanceConstraintTuples = {{A1,M,distA1M},{B1,M,distB1M},{A1,B1,distA1B1},
                            {A2,M,distA2M},{B2,M,distB2M},{A2,B2,distA2B2}};
netList distanceConstraintTuples
end

-----------------------------------------------------------------------------
---------- start M2 in "../../FittingAlgebraicOrbitEquation.m2" -------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Halo";		       -- Type of orbit to fit
jacobiConstantDegree = 3;	       -- Degree in Jacobi constant of the model polynomials
modelDegree = 6;		       -- Degree of the model polynomials

X = ZZ/7772777			       -- Finite field for symbolic computation
Y = X[cA,cB,uA1,vA1,uB1,vB1,vA2,uA2,uB2,vB2]	       -- Polynomial ring over finite field

-- Random measurements
(M,distA1M,distB1M,distA1B1,distA2M,distB2M,distA2B2) = toSequence ({random(X^1, X^3)}|flatten entries random(X^1,X^6));
CPolynomialCoeffMatrixList = apply(2, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));		 -- Random model coeffs
heightCPolynomialCoeffMatrixList = apply(2, i -> random(X^(sub((modelDegree^2/4+1)+modelDegree,ZZ)), X^(jacobiConstantDegree+1))); -- Random height model coeffs

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-One-Mothership-Four-Spacecraft-IOD.m2")
elapsedTime findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree},OrbitScenario => orbitType)
