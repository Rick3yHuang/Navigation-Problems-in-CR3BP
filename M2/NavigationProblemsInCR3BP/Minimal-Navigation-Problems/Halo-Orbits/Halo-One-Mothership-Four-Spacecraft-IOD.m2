-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
A1 = matrix{{xA1},{yA1},{zA1}};
A2 = matrix{{xA2},{yA2},{zA2}};
B1 = matrix{{xB1},{yB1},{zB1}};
B2 = matrix{{xB2},{yB2},{zB2}};
M = transpose M;
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
(CPolynomialCoeffMatrixA,CPolynomialCoeffMatrixB) = toSequence CPolynomialCoeffMatrixList;
(heightCPolynomialCoeffMatrixA,heightCPolynomialCoeffMatrixB) = toSequence heightCPolynomialCoeffMatrixList;
orbitContraintTuples = {{A1,CPolynomialCoeffMatrixA,heightCPolynomialCoeffMatrixA,cA},{A2,CPolynomialCoeffMatrixA,heightCPolynomialCoeffMatrixA,cA},
                        {B1,CPolynomialCoeffMatrixB,heightCPolynomialCoeffMatrixB,cB},{B2,CPolynomialCoeffMatrixB,heightCPolynomialCoeffMatrixB,cB}};
netList orbitContraintTuples
distanceConstraintTuples = {{A1,M,distA1M},{B1,M,distB1M},{A1,B1,distA1B1},
                            {A2,M,distA2M},{B2,M,distB2M},{A2,B2,distA2B2}};
netList distanceConstraintTuples
end

-----------------------------------------------------------------------------
---------- start M2 in "../../FittingAlgebraicOrbitEquation.m2" -------------
-----------------------------------------------------------------------------

-----------------------------------------------------------------------------
---------------- Experiment with fabricated data ----------------------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Halo";		       -- Type of orbit to fit
IDIntervals = {{200,400},{400,600}};   -- ID Intervals for saved coeffs
(orbitIDA,orbitIDB) = (350,450);       -- Orbit IDs for both spacecraft
pointIndicesA = {50,150};	       -- Indices of observed points for spacecraft A
pointIndicesB = {50,150};	       -- Indices of observed points for spacecraft B
jacobiConstantDegree = 3;	       -- Degree in Jacobi constant of the model pol
modelDegree = 6;		       -- Degree of the model polynomials
fabricateFromPointsOnModels = true;    -- Whether to fabricate observed points from the fitted models

R = CC[cA,cB,xA1,yA1,zA1,xB1,yB1,zB1,xA2,yA2,zA2,xB2,yB2,zB2];
-- Data fabrication
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-One-Mothership-Four-Spacecraft-IOD-Fabrication.m2")

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-One-Mothership-Four-Spacecraft-IOD.m2")
system = constructMinimalProblems(orbitContraintTuples,distanceConstraintTuples,R,{jacobiConstantDegree,modelDegree},OrbitScenario => orbitType);
netList system
sub(matrix{system}, matrix(CC,{trueSols}))

-- Solve minimal problem
needsPackage "NumericalAlgebraicGeometry"
elapsedTime sols = solveSystem system;				      	 -- All complex solutions
#sols									 -- Number of complex solutions
realSols = realPoints sols;						 -- Real solutions
s = select(sols, sol -> (norm_2(coordinates sol - trueSols) < 5e-1))	 -- Solutions close to true solution
coordinates first s - trueSols

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
Y = X[cA,cB,xA1,yA1,zA1,xB1,yB1,zB1,xA2,yA2,zA2,xB2,yB2,zB2]	       -- Polynomial ring over finite field

-- Random measurements
(M,distA1M,distB1M,distA1B1,distA2M,distB2M,distA2B2) = toSequence ({random(X^1, X^3)}|flatten entries random(X^1,X^6));
CPolynomialCoeffMatrixList = apply(2, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));		 -- Random model coeffs
heightCPolynomialCoeffMatrixList = apply(2, i -> random(X^(sub((modelDegree^2/4+1)+modelDegree,ZZ)), X^(jacobiConstantDegree+1))); -- Random height model coeffs

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-One-Mothership-Four-Spacecraft-IOD.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree},OrbitScenario => orbitType)






























-*-------------------------------------------------------------------------------------
----------------- Verify the model polynomials directly--------------------------------
---------------------------------------------------------------------------------------
S = CC[x,y,c];
fA = modelPolynomial(S,CPolynomialCoeffMatrixList_0)
sub(fA, sub(observedPointA|matrix{{scaledCList#(orbitID-1)}},CC))
sub(fA, sub(observedPointB|matrix{{scaledCList#(orbitID-1)}},CC))

T = CC[x,y];
gA = baseModelPolynomial(T,allCoeffMatrixA,orbitID)
sub(gA, sub(observedPointA,CC))
sub(gA, sub(observedPointB,CC))
*-------------------------------------------------------------------------------------
