-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
A = matrix{{xA},{yA},{zA}};
B = matrix{{xB},{yB},{zB}};
C = matrix{{xC},{yC},{zC}};
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
(CPolynomialCoeffMatrixA,CPolynomialCoeffMatrixB,CPolynomialCoeffMatrixC) = toSequence CPolynomialCoeffMatrixList;
(heightCPolynomialCoeffMatrixA,heightCPolynomialCoeffMatrixB,heightCPolynomialCoeffMatrixC) = toSequence heightCPolynomialCoeffMatrixList;
orbitContraintTuples = {{A,CPolynomialCoeffMatrixA,heightCPolynomialCoeffMatrixA,scaledCATrue},
                        {B,CPolynomialCoeffMatrixB,heightCPolynomialCoeffMatrixB,scaledCBTrue},
			{C,CPolynomialCoeffMatrixC,heightCPolynomialCoeffMatrixC,scaledCCTrue}};
netList orbitContraintTuples
distanceConstraintTuples = {{A,B,distAB},{A,C,distAC},{B,C,distBC}};
netList distanceConstraintTuples
end

-----------------------------------------------------------------------------
---------- start M2 in "../../FittingAlgebraicOrbitEquation.m2" -------------
-----------------------------------------------------------------------------

-----------------------------------------------------------------------------
---------------- Experiment with fabricated data  ---------------------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Halo";					  -- Type of orbit to fit
IDIntervals = {{200,400},{400,600},{400,600}};		  -- ID Intervals for saved coeffs
(orbitIDA,orbitIDB,orbitIDC) = (350,450,500);		  -- Orbit IDs for two spacecraft
(pointIndexA,pointIndexB,pointIndexC) = (10,100,300);	  -- Indices of observed points for spacecraft A, B, C
modelDegree = 4;		       -- Degree of the model polynomials
jacobiConstantDegree = 3;	       -- Degree in Jacobi constant of the model polynomials
fabricateFromPointsOnModels = true;			  -- Whether to fabricate observed points from the fitted models

R = CC[xA,yA,zA,xB,yB,zB,xC,yC,zC];
-- Data fabrication
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Three-Spacecraft-Positioning-Fabrication.m2")

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Three-Spacecraft-Positioning.m2")
system = constructMinimalProblems(orbitContraintTuples,distanceConstraintTuples,R,{jacobiConstantDegree,modelDegree},OrbitScenario => orbitType);
netList system
sub(matrix{system}, matrix(CC,{trueSols}))

-- Solve minimal problem
needsPackage "NumericalAlgebraicGeometry"
sols = solveSystem system;	  					 -- All complex solutions
#sols									 -- Number of complex solutions (6710 - 02/26/2026)
realSols = realPoints sols;						 -- Real solutions
#realSols								 -- Number of real solutions
s = select(realSols, sol -> (norm_2(coordinates sol - trueSols) < 1e-2)) -- Solutions close to true solution
coordinates first s - trueSols						 -- Error in the found solution

-----------------------------------------------------------------------------
---------------------- Symbolic Computation  --------------------------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Halo";					  -- Type of orbit to fit
jacobiConstantDegree = 3;	       -- Degree in Jacobi constant of the model polynomials
modelDegree = 4;		       -- Degree of the model polynomials

X = ZZ/7772777
Y = X[xA,yA,zA,xB,yB,zB,xC,yC,zC]

-- Random measurements
(distAB,distAC,distBC,scaledCATrue,scaledCBTrue,scaledCCTrue) = toSequence flatten entries random(X^1,X^6);
CPolynomialCoeffMatrixList = apply(3, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));
heightCPolynomialCoeffMatrixList = apply(3, i -> random(X^(sub((modelDegree^2/4+1)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Three-Spacecraft-Positioning.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree},OrbitScenario => orbitType) -- found 256 complex solutions





















-*-------------------------------------------------------------------------------------
----------------- Verify the model polynomials directly--------------------------------
---------------------------------------------------------------------------------------
allCoeffMatrixA = readAllCoefficients (savedModelsDirectory | "Jacobi-Constant-Degree-" | jacobiConstantDegree | "-Model-Degree-" | modelDegree | "/" | orbitType | "-All-Coefficients.txt");
S = CC[xA,yA,c];
fA = modelPolynomial(S,CPolynomialCoeffMatrixList_0,jacobiConstantDegree,modelDegree)
sub(fA, sub(observedPointA|matrix{{scaledCList#(orbitIDA-1)}},CC))

T = CC[xA,yA];
gA = baseModelPolynomial(T,allCoeffMatrixA,orbitIDA,modelDegree)
sub(gA, sub(observedPointA,CC))
*-------------------------------------------------------------------------------------
