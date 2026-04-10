-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
powerBasisList = apply({cA,cB}, x -> transpose matrix{apply(jacobiConstantDegree+1, i -> x^i)})
heightOrbitCoefficientsList = apply(#heightCPolynomialCoeffMatrixList, i -> (heightCPolynomialCoeffMatrixList_i)*(powerBasisList_i))
uvList = {{uA,vA}};
wList = apply(#uvList, i -> (
	(u,v) = toSequence uvList_i;
	modelBasis = returnModelBasis({u,v},modelDegree);
	first flatten entries((matrix{{1}}|matrix{modelBasis})*heightOrbitCoefficientsList_0)
	))
A = matrix{{uA},{uA},{wList_0}};
B = A - distAB*matrix{{uDir},{vDir},{wDir}};
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
(CPolynomialCoeffMatrixA,CPolynomialCoeffMatrixB) = toSequence CPolynomialCoeffMatrixList;
(heightCPolynomialCoeffMatrixA,heightCPolynomialCoeffMatrixB) = toSequence heightCPolynomialCoeffMatrixList;
orbitContraintTuples = {{A,CPolynomialCoeffMatrixA,cA},{B,CPolynomialCoeffMatrixB,cB,heightCPolynomialCoeffMatrixB}};
netList orbitContraintTuples
distanceConstraintTuples = {{}};
netList distanceConstraintTuples
end

-----------------------------------------------------------------------------
------------ start M2 in "../../NavigationProblemsInCR3BP.m2" ---------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Halo";		       -- Type of orbit to fit
jacobiConstantDegree = 3;	       -- Degree in Jacobi constant of the model polynomials
modelDegree = 6;		       -- Degree of the model polynomials

X = ZZ/7772777			       -- Finite field for symbolic computation
Y = X[cA,uA,vA]	       -- Polynomial ring over finite field

-- Random measurements
(distAB,cB,uDir,vDir,wDir) = toSequence flatten entries random(X^1,X^5);
CPolynomialCoeffMatrixList = apply(2, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));		 -- Random model coeffs
heightCPolynomialCoeffMatrixList = apply(2, i -> random(X^(sub((modelDegree^2/4+1)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Two-Spacecraft-Positioning.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree})
