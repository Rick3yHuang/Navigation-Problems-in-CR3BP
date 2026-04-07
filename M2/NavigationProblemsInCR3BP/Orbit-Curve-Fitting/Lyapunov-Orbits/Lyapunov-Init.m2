if trainDirectLeastSquaresModel then (
    data = fold((X,Y) -> X || Y, posAndScaledCList_IDRange);
    numXYVariables = numRows exponentMatrix;
    exponentMatrixListIncludingC = apply(jacobiConstantDegree+1, i -> exponentMatrix | transpose matrix{toList(numXYVariables:i)});
    exponentMatrixIncludingC = fold((X,Y) -> X || Y, exponentMatrixListIncludingC);
    << "-- The exponent matrix including C is: " << exponentMatrixIncludingC << endl;
    coeff = first solveForOrbitCoefficients({data}, exponentMatrixIncludingC);
    CPolynomialCoeffMatrix = reshape(RR^(numXYVariables),RR^(jacobiConstantDegree+1),coeff);
    -- store the results
    storePowersInCResults(directLeastSquaresCoeffFilePath,CPolynomialCoeffMatrix);
    ) else (
    -----------------------------------------------------------------------------
    --------------- Fix C, find coefficients of orbit equations -----------------
    -----------------------------------------------------------------------------
    coeffList = solveForOrbitCoefficients(positionMatrixList, exponentMatrix);

    -----------------------------------------------------------------------------
    ---------------- Peek the orbit equation coefficients -----------------------
    -----------------------------------------------------------------------------
    peekOrbitEquationCoefficients(scaledCList, coeffList);

    -----------------------------------------------------------------------------
    --------------------- Evaluate orbit equation fit ---------------------------
    -----------------------------------------------------------------------------
    --evaluateLeastSquaresFit(positionMatrixList,coeffList,{dataIDList,CList,scaledCList},modelDegree,SolveByOptimization => evaluateByOptimization);

    -- each row is the coefficients of the orbit equation for a specific C
    allCoefficients = transpose listOfVectorsToMatrix(coeffList);
    storeAllCoefficients((savedModelsDirectory | "Jacobi-Constant-Degree-" | jacobiConstantDegree | "-Model-Degree-" | modelDegree| "/" | "Lyapunov-All-Coefficients.txt"),allCoefficients);
    )

