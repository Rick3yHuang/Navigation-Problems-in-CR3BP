-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
A = matrix{{uA},{vA},{wA}};
B = matrix{{uB},{vB},{wB}};
C = matrix{{uC},{vC},{wC}};
D = matrix{{uD},{vD},{wD}};
E = matrix{{uE},{vE},{wE}};
F = matrix{{uF},{vF},{wF}};
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
orbitContraintTuples = apply(6,i -> {{A,B,C,D,E,F}_i, CPolynomialCoeffMatrixList_i, heightCPolynomialCoeffMatrixList_i, (JA,JB,JC,JD,JE,JF)_i});
netList orbitContraintTuples
distanceConstraintTuples = {{A,B,distAB},{A,C,distAC},{A,D,distAD},{B,C,distBC},{B,D,distBD},{C,D,distCD},{A,E,distAE},{B,E,distBE},{C,E,distCE},{A,F,distAF},{B,F,distBF},{C,F,distCF}};
netList distanceConstraintTuples
end

-----------------------------------------------------------------------------
---------- start M2 in "../../FittingAlgebraicOrbitEquation.m2" -------------
-----------------------------------------------------------------------------

-----------------------------------------------------------------------------
---------------- Experiment with fabricated data  ---------------------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0;
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Halo";					     -- Type of orbit to fit
-- ID Intervals for saved coeffs
IDIntervals = {{200,400},{200,400},{400,600},{400,600},{600,800},{600,800}};
(orbitIDA,orbitIDB,orbitIDC,orbitIDD,orbitIDE,orbitIDF) = (250,350,450,500,650,700);
(pointIndexA,pointIndexB,pointIndexC, pointIndexD,pointIndexE,pointIndexF) = (10,100,50,300,200,100);
jacobiConstantDegree = 3;				     -- degree of Jacobi constant polynomial model
modelDegree = 6;					     -- degree of orbit polynomial model
fabricateFromPointsOnModels = true;			     -- Whether to fabricate observed points from the fitted models

R = CC[JA,JB,JC,JD,JE,JF,uA,vA,wA,uB,vB,wB,uC,vC,wC,uD,vD,wD,uE,vE,wE,uF,vF,wF];
-- Data fabrication
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Six-Spacecraft-IOD-Fabrication.m2")

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Six-Spacecraft-IOD.m2")
system = constructMinimalProblems(orbitContraintTuples,distanceConstraintTuples,R,{jacobiConstantDegree,modelDegree},OrbitScenario => orbitType);
netList system
sub(matrix{system}, matrix(CC,{trueSols}))

-- Solve minimal problem
needsPackage "NumericalAlgebraicGeometry"
sols = solveSystem system;						 -- All complex solutions
#sols									 -- Number of complex solutions
realSols = realPoints sols;						 -- Real solutions
s = select(sols, sol -> (norm_2(coordinates sol - trueSols) < 5e-1))	 -- Solutions close to true solution
coordinates first s - trueSols						 -- Error in the found solution

-----------------------------------------------------------------------------
---------------------- Symbolic Computation  --------------------------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Halo";		       -- Type of orbit to fit
jacobiConstantDegree = 3;	       -- Degree in Jacobi constant of the model polynomials
modelDegree = 6;		       -- Degree of the model polynomials

X = ZZ/7772777			       -- Finite field for symbolic computation
-- Polynomial ring over finite field
Y = X[JA,JB,JC,JD,JE,JF,uA,vA,wA,uB,vB,wB,uC,vC,wC,uD,vD,wD,uE,vE,wE,uF,vF,wF]

-- Random measurements
(distAB,distAC,distAD,distBC,distBD,distCD,distAE,distBE,distCE,distAF,distBF,distCF) = toSequence apply(12,i -> random(X));
CPolynomialCoeffMatrixList = apply(6, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));		 -- Random model coeffs
heightCPolynomialCoeffMatrixList = apply(6, i -> random(X^(sub((modelDegree^2/4+1)+modelDegree,ZZ)), X^(jacobiConstantDegree+1))); -- Random height model coeffs

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Six-Spacecraft-IOD.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree}, OrbitScenario => orbitType)
