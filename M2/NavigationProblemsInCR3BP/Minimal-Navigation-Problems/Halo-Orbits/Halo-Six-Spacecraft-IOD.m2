-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
powerBasisList = apply({JA,JB,JC,JD,JE,JF}, x -> transpose matrix{apply(jacobiConstantDegree+1, i -> x^i)})
heightOrbitCoefficientsList = apply(#heightCPolynomialCoeffMatrixList, i -> (heightCPolynomialCoeffMatrixList_i)*(powerBasisList_i))
uvList = {{uA,vA},{uB,vB},{uC,vC},{uD,vD},{uE,vE},{uF,vF}};
wList = apply(#uvList, i -> (
	(u,v) = toSequence uvList_i;
	modelBasis = returnModelBasis({u,v},modelDegree);
	first flatten entries((matrix{{1}}|matrix{modelBasis})*heightOrbitCoefficientsList_i)
	))
A = matrix{{uA},{vA},{wList_0}};
B = matrix{{uB},{vB},{wList_1}};
C = matrix{{uC},{vC},{wList_2}};
D = matrix{{uD},{vD},{wList_3}};
E = matrix{{uE},{vE},{wList_4}};
F = matrix{{uF},{vF},{wList_5}};
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
orbitContraintTuples = apply(6,i -> {{A,B,C,D,E,F}_i, CPolynomialCoeffMatrixList_i, (JA,JB,JC,JD,JE,JF)_i});
netList orbitContraintTuples
distanceConstraintTuples = {{A,B,distAB},{A,C,distAC},{A,D,distAD},{B,C,distBC},{B,D,distBD},{C,D,distCD},{A,E,distAE},{B,E,distBE},{C,E,distCE},{A,F,distAF},{B,F,distBF},{C,F,distCF}};
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
modelDegree = 4;		       -- Degree of the model polynomials

X = ZZ/7772777			       -- Finite field for symbolic computation
-- Polynomial ring over finite field
Y = X[JA,JB,JC,JD,JE,JF,uA,vA,uB,vB,uC,vC,uD,vD,uE,vE,uF,vF]

-- Random measurements
(distAB,distAC,distAD,distBC,distBD,distCD,distAE,distBE,distCE,distAF,distBF,distCF) = toSequence apply(12,i -> random(X));
CPolynomialCoeffMatrixList = apply(6, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));		 -- Random model coeffs
heightCPolynomialCoeffMatrixList = apply(6, i -> random(X^(sub((modelDegree^2/4+1)+modelDegree,ZZ)), X^(jacobiConstantDegree+1))); -- Random height model coeffs

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Six-Spacecraft-IOD.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree})
