-----------------------------------------------------------------------------
-------- start M2 in "../../NavigationProblemsInCR3BP.m2" ---------------
-----------------------------------------------------------------------------
restart
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Halo";			   -- Type of orbit to fit
specificBasis = null;			   -- Use a specific construction of basis
retrainOrbitModel = false;		   -- fit orbit coefficients or load saved data
retrainJacobiConstantModel = false;	   -- fit C coefficients or load saved data
evaluateModel = true;			   -- evaluate the model or not
minID = 400;				   -- min data ID (in the range of 1 to 1146)
maxID = 600;				   -- max data ID (in the range of 2 to 1147)
jacobiConstantDegree = 3;		   -- degree of Jacobi constant polynomial model
modelDegree = 4;			   -- degree of orbit polynomial model
evaluateByOptimization = false;		   -- evaluate the model with optimization distance or not
elapsedTime needs (fittingExperimentDirectory | "Polynomials-In-C-Fitting.m2");
