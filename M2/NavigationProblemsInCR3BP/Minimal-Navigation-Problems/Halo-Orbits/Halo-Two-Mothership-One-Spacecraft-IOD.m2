-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
powerBasisList = apply({c}, x -> transpose matrix{apply(jacobiConstantDegree+1, i -> x^i)})
heightOrbitCoefficientsList = apply(#heightCPolynomialCoeffMatrixList, i -> (heightCPolynomialCoeffMatrixList_i)*(powerBasisList_i))
uvList = {{uA,vA}};
wList = apply(#uvList, i -> (
	(u,v) = toSequence uvList_i;
	modelBasis = returnModelBasis({u,v},modelDegree);
	first flatten entries((matrix{{1}}|matrix{modelBasis})*heightOrbitCoefficientsList_i)
	))
A = matrix{{uA},{vA},{wList_0}};
M1 = transpose M1;
M2 = transpose M2;
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
CPolynomialCoeffMatrix = first CPolynomialCoeffMatrixList;
heightCPolynomialCoeffMatrix = first heightCPolynomialCoeffMatrixList;
orbitContraintTuples = {{A,CPolynomialCoeffMatrix,c}}
netList orbitContraintTuples
distanceConstraintTuples = {{A,M1,distAM1},{A,M2,distAM2}};
netList distanceConstraintTuples
end

-----------------------------------------------------------------------------
---------- start M2 in "../../FittingAlgebraicOrbitEquation.m2" -------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Halo";		       -- Type of orbit to fit
jacobiConstantDegree = 3;	       -- Degree of the Jacobi constant polynomial
modelDegree = 6;		       -- Degree of the orbit model polynomials

X = ZZ/7772777
Y = X[c,uA,vA]

-- Random measurements
(M1,M2,distAM1,distAM2) = (random(X^1, X^3),random(X^1,X^3),random(X),random(X));
CPolynomialCoeffMatrixList = apply(1, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));
heightCPolynomialCoeffMatrixList = apply(1, i -> random(X^(sub((modelDegree^2/4+1)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Two-Mothership-One-Spacecraft-IOD.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree},OrbitScenario => "Halo")
