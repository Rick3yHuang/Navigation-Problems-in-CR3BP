-----------------------------------------------------------------------------
-------- start M2 in "../../NavigationProblemsInCR3BP.m2" ---------------
-----------------------------------------------------------------------------
restart
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Lyapunov";					  -- Type of orbit to fit
specificBasis = null;			  -- Use a specific construction of basis
trainOrbitModel = false;		  -- fit orbit coefficients or load saved data
trainJacobiConstantModel = false;	  -- fit C coefficients or load saved data
trainDirectLeastSquaresModel = true;	  -- fit the direct least-squares model or not
evaluateModel = true;			  -- evaluate the model or not
minID = 400;				  -- min data ID (in the range of 1 to 1037)
maxID = 600;				  -- max data ID (in the range of 1 to 1037)
jacobiConstantDegree = 3;		  -- degree of Jacobi constant polynomial model
modelDegree = 6;			  -- degree of orbit polynomial model
evaluateByOptimization = false;		  -- evaluate the model with optimization distance or not
elapsedTime needs (fittingExperimentDirectory | "Polynomials-In-C-Fitting.m2");
