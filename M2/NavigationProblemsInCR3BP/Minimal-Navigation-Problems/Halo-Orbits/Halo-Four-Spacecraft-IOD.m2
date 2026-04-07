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

-----------------------------------------------------------------------------
---------------- Experiment with fabricated data  ---------------------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Halo";		       -- Type of orbit to fit
IDIntervals = {{200,400},{400,600}};   -- ID Intervals for saved coeffs
(orbitIDA,orbitIDB) = (350,450);       -- Orbit IDs for two spacecraft
pointIndicesA = {1,100};	       -- Indices of observed points for spacecraft A
pointIndicesB = {1,100};	       -- Indices of observed points for spacecraft B
jacobiConstantDegree = 3;	       -- Degree in Jacobi constant of the model polynomials
modelDegree = 4;		       -- Degree of the model polynomials
fabricateFromPointsOnModels = true;    -- Whether to fabricate observed points from the fitted models

R = CC[cA,cB,sA1,sA2,sB1,sB2]
-- Data fabrication
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Four-Spacecraft-IOD-Fabrication.m2")

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Four-Spacecraft-IOD.m2")
system = constructMinimalProblems(orbitContraintTuples,distanceConstraintTuples,R,{jacobiConstantDegree,modelDegree},OrbitScenario => orbitType);
netList system
sub(matrix{system}, matrix(CC,{trueSols}))

-- Solve minimal problem
needsPackage "NumericalAlgebraicGeometry"
sols = solveSystem system;						 -- All complex solutions
#sols									 -- Number of complex solutions
realSols = realPoints sols;						 -- Real solutions
s = select(sols, sol -> (norm_2(coordinates sol - trueSols) < 9e-2))	 -- Solutions close to true solution
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
Y = X[cA,cB,sA1,sA2,sB1,sB2]	       -- Polynomial ring over finite field

-- Random measurements
(dirA1,dirA2,dirB1,dirB2,distAB1,distAB2) = toSequence (apply(4, i -> random(X^1, X^2))|{random(X),random(X)});
CPolynomialCoeffMatrixList = apply(2, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));		 -- Random model coeffs
heightCPolynomialCoeffMatrixList = apply(2, i -> random(X^(sub((modelDegree^2/4+1)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));	 -- Random height model coeffs

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Four-Spacecraft-IOD.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree}, OrbitScenario => orbitType)






















-*-------------------------------------------------------------------------------------
----------------- Verify the model polynomials directly--------------------------------
---------------------------------------------------------------------------------------
allCoeffMatrixA = readAllCoefficients (savedModelsDirectory | "Jacobi-Constant-Degree-" | jacobiConstantDegree | "-Model-Degree-" | modelDegree | "/" | orbitType | "-All-Coefficients.txt");
S = CC[xA,yA,c];
fA = modelPolynomial(S,CPolynomialCoeffMatrixList_0,3,4)
sub(fA, sub(A1|matrix{{scaledCList#(orbitIDA-1)}},CC))
sub(fA, sub(A2|matrix{{scaledCList#(orbitIDA-1)}},CC))

T = CC[xA,yA];
gA = baseModelPolynomial(T,allCoeffMatrixA,orbitIDA,4)
sub(gA, sub(A1,CC))
sub(gA, sub(A2,CC))
*-------------------------------------------------------------------------------------
