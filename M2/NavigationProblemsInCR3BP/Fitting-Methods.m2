--------------------------------------------------------------------------------
-- Least Square Methods for fitting orbit equations and C-polynomials ----------
--------------------------------------------------------------------------------
leastSquareSolution = method()
leastSquareSolution (Matrix,Matrix) := (designMatrix, constantTerm) -> (
    designMatrixT := transpose designMatrix;
    --    pseudoInv := (designMatrixT * designMatrix)^-1 * designMatrixT;
    --    pseudoInv * constantTerm
    solve(designMatrixT * designMatrix, designMatrixT * constantTerm)
    )

--------------------------------------------------------------------------------
-- Solve for orbit coefficients and C coefficients -----------------------------
--------------------------------------------------------------------------------
solveForOrbitCoefficients = method()
solveForOrbitCoefficients (List, Matrix) := (positionMatrixList, exponentMatrix) -> (
    << "-- Constructing design matrix ..." << endl;
    elapsedTime orbitDesignMatrixList := apply(positionMatrixList, positionMatrix -> constructDesignMatrix(positionMatrix, exponentMatrix));
    coeffList := apply(orbitDesignMatrixList, designMatrix -> (
	    ones := transpose matrix {toList (numRows designMatrix : 1)};
	    leastSquareSolution(designMatrix, ones))
	)
    )

solveForHeightCoefficients = method()
solveForHeightCoefficients (List,List,Matrix) := (projectedPositionMatrixList, wPositionMatrixList, exponentMatrix) -> (
    orbitDesignMatrixList := apply(projectedPositionMatrixList, positionMatrix -> constructDesignMatrix(positionMatrix, exponentMatrix));
    coeffList := apply((#orbitDesignMatrixList), i -> (
	    designMatrix := orbitDesignMatrixList_i;
	    w := wPositionMatrixList_i;
	    leastSquareSolution(designMatrix, w))
	)
    )

solveForCCoefficients = method(Options => {UseSpecificBasis => null})
solveForCCoefficients (List, Matrix, ZZ) := o -> (scaledCList, allCoefficients, maxDegree) -> (
    CPolynomialDesignMatrix := constructCPolynomialDesignMatrix(scaledCList, maxDegree, UseSpecificBasis=>o#UseSpecificBasis);
    CPolynomialCoeffList := apply(numColumns allCoefficients, j -> (
	    leastSquareSolution(CPolynomialDesignMatrix, matrix allCoefficients_j)
	    ));
    transpose listOfVectorsToMatrix(CPolynomialCoeffList)
    )

--------------------------------------------------------------------------------
-- Construct Design Matrices for fitting orbit equations and C-polynomials -----
--------------------------------------------------------------------------------
constructDesignMatrix = method()
constructDesignMatrix (Matrix, Matrix) := (dataMatrix, exponentMatrix) -> (
    numVars := numColumns dataMatrix;
    numPoints := numrows dataMatrix;
    numExponents := numrows exponentMatrix;

    if numColumns exponentMatrix != numVars then error "Number of variables and number of exponents do not match";

    variableList := apply(numVars, i -> flatten entries dataMatrix_i);

    designMatrix := mutableMatrix(RR,numPoints,numExponents);
    scan(numExponents, i -> (
	    scan(numPoints, j -> (
		    designMatrix_(j,i) = product(apply(numVars, k -> (variableList_k_j)^(exponentMatrix_(i,k))))
		    )
		)
	    )
	);
    matrix designMatrix
    )

constructCPolynomialDesignMatrix = method(Options => {UseSpecificBasis => null})
constructCPolynomialDesignMatrix (List, ZZ) := o -> (CList, maxDegree) -> (
    designMatrix := mutableMatrix(RR, #CList, maxDegree + 1);
    scan(#CList, i -> (
	    specificBasisValues := null;
	    if instance(o#UseSpecificBasis,String) then specificBasisValues = chebyshevRow(CList_i, maxDegree);
	    designMatrix_(i,0) = 1;
	    scan(maxDegree, j -> (
		    designMatrix_(i,j+1) = (CList_i)^(j+1);
		    if instance(o#UseSpecificBasis,String) then designMatrix_(i,j+1) = specificBasisValues_(j+1)
		    ))
	    ));
    matrix designMatrix    
    )

--------------------------------------------------------------------------------
-- Using other Basis for Design Matrix Construction ----------------------------
--------------------------------------------------------------------------------
-- Chebyshev Basis Functions
-- output: list of Chebyshev polynomial values at C up to maxDegree, each entry is a constant T_k(C)
chebyshevRow = method()
chebyshevRow (RR, ZZ) := (C, maxDegree) -> (
    if maxDegree == 0 then out := {1}
    else if maxDegree == 1 then out {1, C}
    else (
	oldT0 := 1; oldT1 := C; out = {oldT0,oldT1};
	currentT0 := oldT0; currentT1 := oldT1;
	for k from 2 to maxDegree do (
	    currentT1 = 2*C*oldT1 - oldT0;
	    out = append(out, currentT1);
	    currentT0 = oldT1;
	    oldT1 = currentT1;
	    oldT0 = currentT0;
	    );
	out
	)
    )

--------------------------------------------------------------------------------
-- Gradient Descent Methods for Fitting C Polynomials --------------------------
--------------------------------------------------------------------------------

gradientDescentFit = method(Options => {MaxGDSteps => 1000, GDLearningRate => 0.01, GDRelativeErrorTol => 1e-4, InitialGuess => null})
gradientDescentFit(Ring, List, ZZ, ZZ) := o -> (R, orbitDataList, modelDegree, jacobiConstantDegree) -> (
    (x,y,C) := toSequence flatten entries vars R;
    (maxSteps, lr, relativeErrorTol) := (o#MaxGDSteps, o#GDLearningRate, o#GDRelativeErrorTol);
    
    b := returnAllModelBasis({x,y,C}, modelDegree, jacobiConstantDegree);
    <<"-- Basis: " << b << endl;
    
    m := #b;

    --------------------------------------------------------------------------
    -- initial guess
    --------------------------------------------------------------------------
    if class o#InitialGuess === Matrix then (
	theta := o#InitialGuess;
	if numRows theta != m then error "Initial guess has wrong number of rows";
	if numColumns theta != 1 then error "Initial guess should be a column vector";
	) else (
	theta = transpose matrix{toList(m : 0.)};
	);

    data := fold((X,Y) -> X || Y, orbitDataList);
    n := numRows data;
    onesVec := transpose matrix{toList((numRows data) : 1.)};
    --------------------------------------------------------------------------
    -- loss and gradient
    --------------------------------------------------------------------------
    -- Define design matrix A
    exponentMatrix := findModelExponentMatrix modelDegree;
    numXYVariables := numRows exponentMatrix;
    exponentMatrixListIncludingC := apply(jacobiConstantDegree+1, i -> exponentMatrix | transpose matrix{toList(numXYVariables:i)});
    exponentMatrixIncludingC := fold((X,Y) -> X || Y, exponentMatrixListIncludingC);
    A := constructDesignMatrix(data, exponentMatrixIncludingC);
    
    
    lossValue := theta -> (
	r := A * theta - onesVec;
	(1.0/n) * (norm_2(r))^2
	);

    gradientValue := theta -> (
	r := A * theta - onesVec;
	(2.0/n) * transpose(A) * r
	);

    oldLoss := lossValue theta;
    << "-- initial loss = " << oldLoss << endl;

    converged := false;
    stepsTaken := 0;

    scan(maxSteps, i -> (
	    if not converged then (
		grad := gradientValue theta;
		theta = theta - lr * grad;
		newLoss := lossValue theta;

		stepsTaken = i + 1;

		relativeError := abs(newLoss - oldLoss) / (abs(oldLoss));
		if relativeError < relativeErrorTol then (
		    converged = true;
		    << "-- Converged at step " << i+1 << " with relative error " << relativeError << endl;
		    );
		if i % 50 == 0 then (
		    << "-- Step " << i+1 << ": loss = " << newLoss << ", relative error: " << relativeError << endl;
		    );

		oldLoss = newLoss;
		)
	    )
	);
    
    new HashTable from {
	"theta" => flatten entries theta,
	"loss" => oldLoss
	}
    )


-- TODO: Need to take the Halo orbit into account for the height polynomial and the corresponding coefficients
linearGradientDescentFit = method(Options => {MaxGDSteps => 1000, GDLearningRate => 0.01, GDRelativeErrorTol => 1e-4})
linearGradientDescentFit(Ring, List, ZZ, ZZ) := o -> (R, orbitDataList, modelDegree,jacobiConstantDegree) -> (
    (x,y,C) := toSequence flatten entries vars R;
    (maxSteps, lr, relativeErrorTol) := (o#MaxGDSteps, o#GDLearningRate, o#GDRelativeErrorTol);
    
    b := returnAllModelBasis({x,y,C}, modelDegree, jacobiConstantDegree);
    <<"-- Basis: " << b << endl;
    
    m := #b;

    --------------------------------------------------------------------------
    -- initial guess
    --------------------------------------------------------------------------
    theta := transpose matrix{toList(m : 0.)}; -- change this to a better guess (e.g. the solution of the two-step least squares)

    data := fold((X,Y) -> X || Y, orbitDataList);
    onesVec := transpose matrix{toList((numRows data) : 1.)};
    --------------------------------------------------------------------------
    -- loss and gradient
    --------------------------------------------------------------------------
    lossValue := (theta,data) -> (
	--------------------------------------------------------------------------
	-- build design matrix A once using substitution
	-- A_(i,j) = b_j evaluated at data point i
	--------------------------------------------------------------------------
	n := numRows data;
	A := matrix apply(n, i -> (
		p := data^{i};
		apply(b, f -> sub(f,p))
		));
	r := A * theta - onesVec;
	(1.0/n) * (norm_2(r))^2
	);

    gradientValue := (theta,data) -> (
	n := numRows data;
	A := matrix apply(n, i -> (
		p := data^{i};
		apply(b, f -> sub(f,p))
		));
	r := A * theta - onesVec;
	(2.0/n) * transpose(A) * r
	);

    oldLoss := lossValue(theta,data);
    << "-- initial loss = " << oldLoss << endl;

    eps := 1e-12;
    converged := false;
    stepsTaken := 0;

    scan(maxSteps, i -> (
	    grad := gradientValue(theta,data);
	    theta = theta - lr * grad;
	    newLoss := lossValue(theta,data);

	    stepsTaken = i + 1;

	    if i % 50 == 0 then (
		<< "-- Step " << i+1 << ": loss = " << newLoss               
		);

	    oldLoss = newLoss;
	    ));

    new HashTable from {
	"theta" => flatten entries theta,
	"loss" => oldLoss
	}
    )

-- Gradient Descent with geometric loss
