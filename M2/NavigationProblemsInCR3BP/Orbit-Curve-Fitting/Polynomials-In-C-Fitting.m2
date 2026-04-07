-----------------------------------------------------------------------------
---- Shortcuts---------------------------------------------------------------
-----------------------------------------------------------------------------
modelSaveDirectory = savedModelsDirectory | "Jacobi-Constant-Degree-" | jacobiConstantDegree | "-Model-Degree-" | modelDegree | "/" | orbitType;
coeffFilePath = modelSaveDirectory | "-C-Polynomial-Coeff-Output-" | minID | "-to-" | maxID | ".txt";
heightCoeffFilePath = modelSaveDirectory | "-C-Polynomial-Height-Coeff-Output-" | minID | "-to-" | maxID | ".txt";
directLeastSquaresCoeffFilePath = modelSaveDirectory | "-Direct-Least-Squares-Coeff-Output-" | minID | "-to-" | maxID | ".txt";
allCoefficientsFilePath = modelSaveDirectory | "/" | orbitType | "-All-Coefficients.txt";
allHeightCoefficientsFilePath = modelSaveDirectory | "/" | orbitType | "-All-Height-Coefficients.txt";
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
if not trainJacobiConstantModel then(
    elapsedTime needs (fittingExperimentDirectory | "Load-Data.m2");
    ) else (
    scaledCList = loadScaledCList (dataDirectory | orbitType | "-Data/ID_Jacobi_Constants.txt");
    );
-----------------------------------------------------------------------------
---- If asked, train the coefficients of orbit equation -------------------
-----------------------------------------------------------------------------
if trainOrbitModel or trainDirectLeastSquaresModel then (
    needs (fittingDirectory | orbitType | "-Init.m2")
    ) else (
    allCoefficients = readAllCoefficients allCoefficientsFilePath;
    if isHalo then (
	allHeightCoefficients = readAllCoefficients allHeightCoefficientsFilePath;
	);
    )
-----------------------------------------------------------------------------
---- If asked, fit polynomials in C -----------------------------------------
-----------------------------------------------------------------------------
if trainJacobiConstantModel then (
    CPolynomialCoeffMatrix = solveForCCoefficients(scaledCList_IDRange,allCoefficients^IDRange,jacobiConstantDegree,UseSpecificBasis=>specificBasis);
    -- store the results
    storePowersInCResults(coeffFilePath,CPolynomialCoeffMatrix);
    if isHalo then (
	CPolynomialHeightCoeffMatrix = solveForCCoefficients(scaledCList_IDRange,allHeightCoefficients^IDRange,jacobiConstantDegree,UseSpecificBasis=>specificBasis);
	-- store the results
	storePowersInCResults(heightCoeffFilePath,CPolynomialHeightCoeffMatrix);
	);
    ) else if (not trainDirectLeastSquaresModel) then (
    -- read the results directly
    CPolynomialCoeffMatrix = readPowersInCResults coeffFilePath;
    if isHalo then (
	CPolynomialHeightCoeffMatrix = readPowersInCResults heightCoeffFilePath;
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
