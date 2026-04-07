-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
A = matrix{{uA},{vA},{wA}};
B = matrix{{uB},{vB},{wB}};
M = transpose M
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
orbitContraintTuples = apply(2,i -> {{A,B}_i, CPolynomialCoeffMatrixList_0, heightCPolynomialCoeffMatrixList_0, c});
netList orbitContraintTuples
distanceConstraintTuples = {{A,B,distAB},{A,M,distAM},{B,M,distBM}};
netList distanceConstraintTuples
end

-----------------------------------------------------------------------------
---------- start M2 in "../../FittingAlgebraicOrbitEquation.m2" -------------
-----------------------------------------------------------------------------

-----------------------------------------------------------------------------
---------------- Experiment with fabricated data  ---------------------------
-----------------------------------------------------------------------------
restart
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Halo";					  -- Type of orbit to fit
IDIntervals = {{400,600}};				  -- ID Intervals for saved coeffs
orbitID = 450;						  -- Orbit IDs for two spacecraft
pointIndexA = 50;					  -- Index of observed point for the spacecraft A
pointIndexB = 135;					  -- Index of observed point for the spacecraft B
pointIndexList = {pointIndexA,pointIndexB};    		  -- List of point indices for the two spacecraft
jacobiConstantDegree = 3;				  -- degree of Jacobi constant polynomial model
modelDegree = 4;					  -- degree of orbit polynomial model
fabricateFromPointsOnModels = true;			  -- Whether to fabricate observed points from the fitted models

R = CC[c,uA,vA,wA,uB,vB,wB];
-- Data fabrication
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-One-Mothership-Two-Spacecraft-IOD-Fabrication.m2")

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/"  | orbitType | "-One-Mothership-Two-Spacecraft-IOD.m2")
system = constructMinimalProblems(orbitContraintTuples,distanceConstraintTuples,R,{jacobiConstantDegree,modelDegree},OrbitScenario => orbitType);
netList system
sub(matrix{system}, matrix(CC,{trueSols}))

-- Solve minimal problem
needsPackage "NumericalAlgebraicGeometry"
sols = solveSystem system;						-- All complex solutions
#sols									-- Number of complex solutions
realSols = realPoints sols    						-- Real solutions
#realSols;								-- Number of real solutions
s = select(sols, sol -> (norm_2(coordinates sol - trueSols) < 1e-2))	-- Solutions close to true solution
coordinates first s - trueSols	      					-- Error in the found solution

-----------------------------------------------------------------------------
---------------------- Symbolic Computation  --------------------------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Halo";		       -- Type of orbit to fit
jacobiConstantDegree = 3;	       -- Degree of the Jacobi constant polynomial
modelDegree = 4;		       -- Degree of the orbit model polynomials
fabricateFromPointsOnModels = true;    -- Whether to fabricate observed points from the fitted models

X = ZZ/7772777
Y = X[c,uA,vA,wA,uB,vB,wB]

-- Random measurements
(M,distAM,distBM,distAB) = (random(X^1, X^3),random(X),random(X),random(X));
CPolynomialCoeffMatrixList = apply(1, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));
heightCPolynomialCoeffMatrixList = apply(1, i -> random(X^(sub((modelDegree^2/4+1)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-One-Mothership-Two-Spacecraft-IOD.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree},OrbitScenario => orbitType)
