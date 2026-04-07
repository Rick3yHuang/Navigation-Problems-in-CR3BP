-----------------------------------------------------------------------------
--------------- Describe Spacecraft Positions  ------------------------------
-----------------------------------------------------------------------------
A = matrix{{xA},{yA}};
B = matrix{{xB},{yB}};
M = transpose M;
vA = matrix{{vxA},{vyA}};
vB = matrix{{vxB},{vyB}};
vM = matrix{{0},{0}};
sqrtDistA = matrix{{r1A},{r2A}};
sqrtDistB = matrix{{r1B},{r2B}};
-----------------------------------------------------------------------------
--------------------- Build Minimal Problem  --------------------------------
-----------------------------------------------------------------------------
(CPolynomialCoeffMatrixA,CPolynomialCoeffMatrixB) = toSequence CPolynomialCoeffMatrixList;
orbitContraintTuples = {{A,vA,CPolynomialCoeffMatrixA,cA,massRatioEarthMoon,sqrtDistA},
                        {B,vB,CPolynomialCoeffMatrixB,cB,massRatioEarthMoon,sqrtDistB}};
netList orbitContraintTuples
distanceConstraintTuples = {{A,M,vA,vM,distAM,distDerivativeAM},
                            {B,M,vB,vM,distBM,distDerivativeBM},
			    {A,B,vA,vB,distAB,distDerivativeAB}};
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
IDIntervals = {{200,400},{400,600}};   -- ID Intervals for saved coeffs
(orbitIDA,orbitIDB) = (350,450);       -- Orbit IDs for both spacecraft
pointIndexA = 50;		       -- Indices of observed points for spacecraft A
pointIndexB = 135;		       -- Indices of observed points for spacecraft B
jacobiConstantDegree = 3;	       -- Degree in Jacobi constant of the model pol
modelDegree = 4;		       -- Degree of the model polynomials
fabricateFromPointsOnModels = true;    -- Whether to fabricate observed points from the fitted models

R = CC[cA,cB,xA,yA,vxA,vyA,xB,yB,vxB,vyB,r1A,r2A,r1B,r2B]	       -- Polynomial ring for minimal problem;
-- Data fabrication
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-One-Mothership-Two-Spacecraft-Distinct-Orbits-IOD-Fabrication.m2")

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-One-Mothership-Two-Spacecraft-Distinct-Orbits-IOD.m2")
system = constructMinimalProblems(orbitContraintTuples,distanceConstraintTuples,R,{jacobiConstantDegree,modelDegree},DerivativesOfDistances => true);
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
orbitType = "Lyapunov";		       -- Type of orbit to fit
jacobiConstantDegree = 3;	       -- Degree in Jacobi constant of the model polynomials
modelDegree = 4;		       -- Degree of the model polynomials

X = ZZ/7772777			       -- Finite field for symbolic computation
Y = X[cA,cB,xA,yA,vxA,vyA,xB,yB,vxB,vyB,r1A,r2A,r1B,r2B]	       -- Polynomial ring over finite field

-- Random measurements
(M,distA1M,distB1M,distA1B1,distA2M,distB2M,distA2B2) = toSequence ({random(X^1, X^2)}|flatten entries random(X^1,X^6));
CPolynomialCoeffMatrixList = apply(2, i -> random(X^(sub((modelDegree^2/4)+modelDegree,ZZ)), X^(jacobiConstantDegree+1)));		 -- Random model coeffs

-- Build minimal problem
needs (minimalProblemDirectory | orbitType | "-Orbits/" | orbitType | "-One-Mothership-Four-Spacecraft-IOD.m2")
findDegree(orbitContraintTuples,distanceConstraintTuples,Y,{jacobiConstantDegree,modelDegree})






























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
