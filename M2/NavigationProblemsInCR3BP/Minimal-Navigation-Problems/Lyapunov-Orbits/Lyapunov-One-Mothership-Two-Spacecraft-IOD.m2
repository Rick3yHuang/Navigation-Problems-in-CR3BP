-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
A = matrix{{xA},{yA}};
B = matrix{{xB},{yB}};
M = transpose M;
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
CPolynomialCoeffMatrix = first CPolynomialCoeffMatrixList;
orbitContraintTuples = {{A,CPolynomialCoeffMatrix,c},
                        {B,CPolynomialCoeffMatrix,c}};
netList orbitContraintTuples
distanceConstraintTuples = {{A,M,distAM},{B,M,distBM},{A,B,distAB}};
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
orbitID = 450;			       -- Orbit IDs for both spacecraft
pointIndexA = 50;	               -- Index of observed point for the spacecraft A
pointIndexB = 135;	               -- Index of observed point for the spacecraft B
jacobiConstantDegree = 3;	       -- Degree of the Jacobi constant polynomial
modelDegree = 4;		       -- Degree of the orbit model polynomials
fabricateFromPointsOnModels = true;    -- Whether to fabricate observed points from the fitted models
fabricateOneExtraPoint = true;	       

R = CC[c,xA,yA,xB,yB]
-- Data fabrication
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-One-Mothership-Two-Spacecraft-IOD-Fabrication.m2")

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-One-Mothership-Two-Spacecraft-IOD.m2")
system = constructMinimalProblems(orbitContraintTuples,distanceConstraintTuples,R,{jacobiConstantDegree,modelDegree});
netList system
sub(matrix{system}, matrix(CC,{trueSols}))

-- Solve minimal problem
needsPackage "NumericalAlgebraicGeometry"
sols = solveSystem system;						-- All complex solutions
#sols									-- Number of complex solutions
realSols = realPoints sols    						-- Real solutions
#realSols								-- Number of real solutions
S = CC[x,y,c];
fA = modelPolynomial(S,CPolynomialCoeffMatrixList_0,jacobiConstantDegree,modelDegree)
p = select(realSols, sol -> (abs sub(fA, sub(P|matrix{{(coordinates sol)_0}},CC)) < 1e-10));
apply(p, sol -> sub(fA, sub(P|matrix{{(coordinates sol)_0}},CC)))
#p
s = select(sols, sol -> (norm_2(coordinates sol - trueSols) < 1e-2))	-- Solutions close to true solution
coordinates first s - trueSols	      					-- Error in the found solution
s === p

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
Y = X[c,xA,yA,xB,yB]

-- Random measurements
(M,distAM,distBM,distAB) = (random(X^1, X^2),random(X),random(X),random(X));
CPolynomialCoeffMatrixList = apply(1, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-One-Mothership-Two-Spacecraft-IOD.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree}) -- found 84 complex solutions




















-*-------------------------------------------------------------------------------------
----------------- Verify the model polynomials directly--------------------------------
---------------------------------------------------------------------------------------
allCoeffMatrixA = readAllCoefficients (savedModelsDirectory | "Jacobi-Constant-Degree-" | jacobiConstantDegree | "-Model-Degree-" | modelDegree | "/" | orbitType | "-All-Coefficients.txt");
S = CC[xA,yA,c];
fA = modelPolynomial(S,CPolynomialCoeffMatrixList_0,jacobiConstantDegree,modelDegree)
sub(fA, sub(observedPointA|matrix{{scaledCList#(orbitID-1)}},CC))
sub(fA, sub(observedPointB|matrix{{scaledCList#(orbitID-1)}},CC))

T = CC[xA,yA];
gA = baseModelPolynomial(T,allCoeffMatrixA,orbitID,modelDegree)
sub(gA, sub(observedPointA,CC))
sub(gA, sub(observedPointB,CC))
*-------------------------------------------------------------------------------------
