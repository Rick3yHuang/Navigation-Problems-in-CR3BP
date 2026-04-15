-----------------------------------------------------------------------------
-------------- Fix C, find coefficients of f(u,v) = 0  ----------------------
-----------------------------------------------------------------------------
uvCoeffList = solveForOrbitCoefficients(projectedPositionMatrixList, exponentMatrix);
-----------------------------------------------------------------------------
-------------- Fix C, find coefficients of g(u,v) = w  ----------------------
-----------------------------------------------------------------------------
wExponentMatrix := matrix{{0,0}} || exponentMatrix;
heightCoeffList = solveForHeightCoefficients(projectedPositionMatrixList, wPositionMatrixList, wExponentMatrix);

-----------------------------------------------------------------------------
------------------------- Peek orbit coefficients ---------------------------
-----------------------------------------------------------------------------
peekOrbitEquationCoefficients(scaledCList, uvCoeffList);
peekOrbitEquationCoefficients(scaledCList, heightCoeffList);

-----------------------------------------------------------------------------
--------------------- Evaluate orbit equations fit --------------------------
-----------------------------------------------------------------------------
evaluateLeastSquaresFit(projectedPositionMatrixList,uvCoeffList,{dataIDList,CList,scaledCList},modelDegree,SolveByOptimization => evaluateByOptimization);
evaluateLeastSquaresFit(transformedPositionMatrixList,heightCoeffList,{dataIDList,CList,scaledCList},modelDegree,EquationType => "Surface",SolveByOptimization => evaluateByOptimization);

-- each row is the coefficients of the orbit equation for a specific C
allCoefficients = transpose listOfVectorsToMatrix(uvCoeffList);
storeAllCoefficients((savedModelsDirectory | "Halo-All-Coefficients.txt"),allCoefficients);

allHeightCoefficients = transpose listOfVectorsToMatrix(heightCoeffList);
storeAllCoefficients((savedModelsDirectory | "Halo-All-Height-Coefficients.txt"),allHeightCoefficients);
