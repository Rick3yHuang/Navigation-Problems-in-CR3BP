-----------------------------------------------------------------------------
-------- start M2 in "../../NavigationProblemsInCR3BP.m2" ---------------
-----------------------------------------------------------------------------
restart
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Lyapunov";					  -- Type of orbit to fit
specificBasis = null;					  -- Use a specific construction of basis
-- Gradient descent fitting parameters
maxSteps = 1000;					  -- max number of steps for gradient descent
learningRate = 1e-1;					  -- learning rate for gradient descent
relativeErrorTol = 1e-7;				  -- relative error tolerance for gradient descent convergence
initialGuessIDRange = {400,600};			  -- initial guess for the coefficients (if null, will use random initial guess)

evaluateModel = false;					  -- evaluate the model or not
evaluateOnly = false;					  -- only evaluate the model, do not fit the coefficients
linearGD = true;					  -- use linear gradient descent or not
minID = 400;						  -- min data ID (in the range of 1 to 1037)
maxID = 600;						  -- max data ID (in the range of 1 to 1037)
jacobiConstantDegree = 3;				  -- degree of Jacobi constant polynomial model
modelDegree = 4;					  -- degree of orbit polynomial model
evaluateByOptimization = false;				  -- evaluate the model with optimization distance or not
elapsedTime needs (fittingExperimentDirectory | "Gradient-Descent-Fitting.m2");
