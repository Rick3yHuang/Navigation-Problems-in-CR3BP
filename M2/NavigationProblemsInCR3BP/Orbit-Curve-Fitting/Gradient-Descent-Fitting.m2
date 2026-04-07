-----------------------------------------------------------------------------
---- Shortcuts---------------------------------------------------------------
-----------------------------------------------------------------------------
modelSaveDirectory = savedModelsDirectory | "Jacobi-Constant-Degree-" | jacobiConstantDegree | "-Model-Degree-" | modelDegree | "/" | orbitType;
coeffGDFilePath = modelSaveDirectory | "-C-Polynomial-Coeff-Output-From-GD-" | minID | "-to-" | maxID | ".txt";
heightCoeffGDFilePath = modelSaveDirectory | "-C-Polynomial-Height-Coeff-Output-From-GD-" | minID | "-to-" | maxID | ".txt";
InitialGuessFilePath = modelSaveDirectory | "-C-Polynomial-Coeff-Output-" | initialGuessIDRange_0 | "-to-" | initialGuessIDRange_1 | ".txt";
isHalo = orbitType == "Halo";
isLyapunov = orbitType == "Lyapunov";
-----------------------------------------------------------------------------
---- Data Preparation -------------------------------------------------------
-----------------------------------------------------------------------------
IDRange = {(minID-1)..(maxID-1)};
fittingDirectory = fittingExperimentDirectory | orbitType | "-Orbits/";
-----------------------------------------------------------------------------
---- Load Data --------------------------------------------------------------
-----------------------------------------------------------------------------
elapsedTime needs (fittingExperimentDirectory | "Load-Data.m2");
initialGuess = transpose matrix{flatten entries transpose readPowersInCResults InitialGuessFilePath};

-----------------------------------------------------------------------------
---- fit with gradient descent ----------------------------------------------
-----------------------------------------------------------------------------
if linearGD then (
    GDOutput = gradientDescentFit(RR[x,y,c],posAndScaledCList_IDRange,modelDegree,jacobiConstantDegree,
	MaxGDSteps => maxSteps,GDLearningRate => learningRate,GDRelativeErrorTol => relativeErrorTol,
	InitialGuess => initialGuess);
    )
<< "-- Gradient descent loss: " << GDOutput#"loss" << endl;
<< "-- Theta values at final step: " << GDOutput#"theta" << endl;
CPolynomialCoeffMatrix = reshape(RR^(floor((#(GDOutput#"theta"))/(jacobiConstantDegree+1))),RR^(jacobiConstantDegree+1),matrix{GDOutput#"theta"});
-- store the results
storePowersInCResults(coeffGDFilePath,CPolynomialCoeffMatrix);


if evaluateOnly then (
    -- read the results directly
    CPolynomialCoeffMatrix = readPowersInCResults coeffGDFilePath;
    if isHalo then (
	CPolynomialHeightCoeffMatrix = readPowersInCResults heightCoeffGDFilePath;
	);
    )

-----------------------------------------------------------------------------
-------------- Prepare for evaluating the fitted equation -------------------
-----------------------------------------------------------------------------
if evaluateModel then (
    -- Evaluate with coefficients approximated by polynomials in C
    	CPolynomialDesignMatrix = constructCPolynomialDesignMatrix(scaledCList_IDRange, jacobiConstantDegree, UseSpecificBasis=>specificBasis);
	allCoefficientsApprox = CPolynomialDesignMatrix * (transpose CPolynomialCoeffMatrix);
	coeffListApprox = apply(numRows allCoefficientsApprox, i -> transpose allCoefficientsApprox^{i});
    if orbitType == "Lyapunov" then (
	elapsedTime evaluateLeastSquaresFit(positionMatrixList_IDRange,coeffListApprox,{dataIDList_IDRange,CList_IDRange,scaledCList_IDRange},modelDegree,SolveByOptimization => evaluateByOptimization);
	);
    if orbitType == "Halo" then (
	allHeightCoefficientsApprox = CPolynomialDesignMatrix * (transpose CPolynomialHeightCoeffMatrix);
	heightCoeffListApprox = apply(numRows allHeightCoefficientsApprox, i -> transpose allHeightCoefficientsApprox^{i});
	elapsedTime evaluateLeastSquaresFit(projectedPositionMatrixList_IDRange,coeffListApprox,{dataIDList_IDRange,CList_IDRange,scaledCList_IDRange},modelDegree,SolveByOptimization => evaluateByOptimization);
	elapsedTime evaluateLeastSquaresFit(transformedPositionMatrixList_IDRange,heightCoeffListApprox,{dataIDList_IDRange,CList_IDRange,scaledCList_IDRange},modelDegree,EquationType => "Surface",SolveByOptimization => evaluateByOptimization);
	);
    )
