-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
powerBasisList = transpose matrix{apply(jacobiConstantDegree+1, i -> cA^i)}
heightOrbitCoefficientA = (heightCPolynomialCoeffMatrixList_0)*(powerBasisList_0)
uvList = {{uA,vA}};
wList = apply(#uvList, i -> (
	(u,v) = toSequence uvList_i;
	modelBasis = returnModelBasis({u,v},modelDegree);
	first flatten entries((matrix{{1}}|matrix{modelBasis})*heightOrbitCoefficientA)
	))
A = matrix{{uA},{vA},{wList_0}};
B = A - distAB*matrix{{uDirAB},{vDirAB},{wDirAB}};
C = A - distAC*matrix{{uDirAC},{vDirAC},{wDirAC}};
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
(CPolynomialCoeffMatrixA,CPolynomialCoeffMatrixB,CPolynomialCoeffMatrixC) = toSequence CPolynomialCoeffMatrixList;
(heightCPolynomialCoeffMatrixA,heightCPolynomialCoeffMatrixB,heightCPolynomialCoeffMatrixC) = toSequence heightCPolynomialCoeffMatrixList;
orbitContraintTuples = {{A,CPolynomialCoeffMatrixA,cA},{B,CPolynomialCoeffMatrixB,cB,heightCPolynomialCoeffMatrixB},{C,CPolynomialCoeffMatrixC,cC,heightCPolynomialCoeffMatrixC}};
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
Y = X[uA,vA,cA,cB,cC]	       -- Polynomial ring over finite field

-- Random measurements
(distAB,distAC,distBC,uDirAB,vDirAB,wDirAB,uDirAC,vDirAC,wDirAC) = toSequence flatten entries random(X^1,X^9);
(cosThetaAB, cosThetaAC, cosPhiAB, cosPhiAC) = toSequence flatten entries random(X^1,X^4);
CPolynomialCoeffMatrixList = apply(3, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));		 -- Random model coeffs
heightCPolynomialCoeffMatrixList = apply(3, i -> random(X^(sub((modelDegree^2/4+1)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Three-Spacecraft-LOS-Positioning.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree})
