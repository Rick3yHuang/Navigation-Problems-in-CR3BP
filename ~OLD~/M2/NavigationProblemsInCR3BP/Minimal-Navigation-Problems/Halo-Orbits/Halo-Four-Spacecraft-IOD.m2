-- This is underdetermined case, need more modifications
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
(heightCPolynomialCoeffMatrixA,heightCPolynomialCoeffMatrixB) = toSequence heightCPolynomialCoeffMatrixList;
orbitContraintTuples = {{A1,CPolynomialCoeffMatrixA,heightCPolynomialCoeffMatrixA,cA},{A2,CPolynomialCoeffMatrixA,heightCPolynomialCoeffMatrixA,cA},
			{B1,CPolynomialCoeffMatrixB,heightCPolynomialCoeffMatrixB,cB},{B2,CPolynomialCoeffMatrixB,heightCPolynomialCoeffMatrixB,cB}};
netList orbitContraintTuples
distanceConstraintTuples = {{A1,B2,distAB1},{A2,B2,distAB2}};
netList distanceConstraintTuples
end

-----------------------------------------------------------------------------
---------- start M2 in "../../FittingAlgebraicOrbitEquation.m2" -------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0
needsPackage "FittingAlgebraicOrbitEquations"
orbitType = "Halo";		       -- Type of orbit to fit
jacobiConstantDegree = 3;	       -- Degree in Jacobi constant of the model polynomials
modelDegree = 6;		       -- Degree of the model polynomials

X = ZZ/7772777			       -- Finite field for symbolic computation
Y = X[cA,cB,sA1,sA2,sB1,sB2]	       -- Polynomial ring over finite field

-- Random measurements
(dirA1,dirA2,dirB1,dirB2,distAB1,distAB2) = toSequence (apply(4, i -> random(X^1, X^2))|{random(X),random(X)});
CPolynomialCoeffMatrixList = apply(2, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));		 -- Random model coeffs
heightCPolynomialCoeffMatrixList = apply(2, i -> random(X^(sub((modelDegree^2/4+1)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));	 -- Random height model coeffs

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Four-Spacecraft-IOD.m2")
elapsedTime findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree}, OrbitScenario => orbitType)
