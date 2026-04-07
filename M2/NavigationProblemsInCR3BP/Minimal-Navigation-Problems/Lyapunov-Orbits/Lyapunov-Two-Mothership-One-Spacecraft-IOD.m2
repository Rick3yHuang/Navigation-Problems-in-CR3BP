-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
A = matrix{{xA},{yA}};
M1 = transpose M1;
M2 = transpose M2;
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
CPolynomialCoeffMatrix = first CPolynomialCoeffMatrixList;
orbitContraintTuples = {{A,CPolynomialCoeffMatrix,c}}
netList orbitContraintTuples
distanceConstraintTuples = {{A,M1,distAM1},{A,M2,distAM2}};
netList distanceConstraintTuples
end

-----------------------------------------------------------------------------
---------- start M2 in "../../NavigationProblemsInCR3BP.m2" -------------
-----------------------------------------------------------------------------

-----------------------------------------------------------------------------
---------------- Experiment with fabricated data ----------------------------
-----------------------------------------------------------------------------
restart
setRandomSeed 0
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Lyapunov";		       -- Type of orbit to fit
IDIntervals = {{400,600}};	       -- ID Intervals for saved coeffs
orbitID = 420;			       -- Orbit IDs for the spacecraft
pointIndices = {50};		       -- Indices of observed points for the spacecraft
jacobiConstantDegree = 3;	       -- Degree of the Jacobi constant polynomial
modelDegree = 6;		       -- Degree of the orbit model polynomials
fabricateFromPointsOnModels = true;    -- Whether to fabricate observed points from the fitted models

R = CC[c,xA,yA]
-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Two-Mothership-One-Spacecraft-IOD-Fabrication.m2")

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Two-Mothership-One-Spacecraft-IOD.m2")
system = constructMinimalProblems(orbitContraintTuples,distanceConstraintTuples,R,{jacobiConstantDegree,modelDegree});
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
orbitType = "Lyapunov";		       -- Type of orbit to fit
jacobiConstantDegree = 3;	       -- Degree of the Jacobi constant polynomial
modelDegree = 6;		       -- Degree of the orbit model polynomials
fabricateFromPointsOnModels = true;    -- Whether to fabricate observed points from the fitted models

X = ZZ/7772777
Y = X[c,xA,yA]

-- Random measurements
(M1,M2,distAM1,distAM2) = (random(X^1, X^2),random(X^1,X^2),random(X),random(X));
CPolynomialCoeffMatrixList = apply(1, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-Two-Mothership-One-Spacecraft-IOD.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree}) -- found 84 complex solutions




















-*-------------------------------------------------------------------------------------
----------------- Verify the model polynomials directly--------------------------------
---------------------------------------------------------------------------------------
allCoeffMatrixA = readAllCoefficients (savedModelsDirectory | "Jacobi-Constant-Degree-" | jacobiConstantDegree | "-Model-Degree-" | modelDegree | "/" | orbitType | "-All-Coefficients.txt");
S = CC[xA,yA,c];
fA = modelPolynomial(S,CPolynomialCoeffMatrixList_0,jacobiConstantDegree,modelDegree)
sub(fA, sub(A|matrix{{scaledCList#(orbitID-1)}},CC))

T = CC[xA,yA];
gA = baseModelPolynomial(T,allCoeffMatrixA,orbitID,modelDegree)
sub(gA, sub(A,CC))
*-------------------------------------------------------------------------------------
