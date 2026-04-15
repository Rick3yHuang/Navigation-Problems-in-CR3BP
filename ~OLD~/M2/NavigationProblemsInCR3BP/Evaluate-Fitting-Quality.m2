needsPackage "NumericalAlgebraicGeometry"

evaluateLeastSquaresFit = method(Options => {EquationType => "Curve",SolveByOptimization => false})
evaluateLeastSquaresFit (List,List,List,ZZ) := o -> (positionMatrixList,coeffList,IDandCList,modelDegree) -> (
    (IDList, scaledCList, CList) := toSequence IDandCList;
    << "-*" << endl;
    assert (#scaledCList > 1);
    n := floor(#scaledCList/2);
    << "Printing the evaluation results for the two orbits in the middle:" << endl;
    scan(2, i -> (
	    pointMatrix := positionMatrixList_(n+i);
	    << "--------------------------------------------------------" << endl;
	    << "Data ID: " << IDList_(n+i) << ", scaled C: " << CList_(n+i)  << ", C: " << scaledCList_(n+i) << endl;
	    << "Number of data points: " << numRows pointMatrix << endl;
	    (meanApproxDistance, maxApproxDistance,RMSE) := computeAllApproxDistanceStats(pointMatrix, matrix coeffList_(n+i), modelDegree,EquationType => o#EquationType);
	    << "Mean approximated distance: " << meanApproxDistance << endl;
	    << "Max approximated distance: " << maxApproxDistance << endl;
	    << "Root Mean Squared Error from orbit: " << RMSE << endl;
	    if o#SolveByOptimization then (
		(meanOptimizationDistance, maxOptimizationDistance) := computeAllOptimizationDistanceStats(pointMatrix, matrix coeffList_(n+i), modelDegree,EquationType => o#EquationType);
		<< "Mean distance from optimization: " << meanOptimizationDistance << endl;
		<< "Max distance from optimization: " << maxOptimizationDistance << endl;
		);
	    ));
    << "--------------------------------------------------------" << endl;
    << "*-" << endl;
    )

computeApproxDistance = method(Options => {EquationType => "Curve"})
computeApproxDistance (Matrix,Matrix,ZZ) := o -> (pointMatrix,coeffMatrix,modelDegree) -> (
    xVal := pointMatrix_(0,0);
    yVal := pointMatrix_(0,1);
    x := local x;
    y := local y;
    R := RR[x,y];
    f := returnModelBasis({x,y}, modelDegree);
    fx := apply(f, monom -> diff(x, monom));
    fy := apply(f, monom -> diff(y, monom));
    fM := matrix{apply(f,monom -> sub(monom,matrix{{xVal,yVal}}))};
    fxM := matrix{apply(fx,monom -> sub(monom,matrix{{xVal,yVal}}))};
    fyM := matrix{apply(fy,monom -> sub(monom,matrix{{xVal,yVal}}))};
    if o#EquationType == "Curve" then (
	out := abs(first flatten entries (fM*coeffMatrix)-1)/sqrt((first flatten entries (fxM*coeffMatrix))^2 + (first flatten entries (fyM*coeffMatrix))^2);
	) else if o#EquationType == "Surface" then (
	zVal := pointMatrix_(0,2);
	out = abs(first flatten entries ((matrix{{1}}|fM)*coeffMatrix)-zVal)/sqrt((first flatten entries ((matrix{{1}}|fxM)*coeffMatrix))^2 + (first flatten entries ((matrix{{1}}|fyM)*coeffMatrix))^2 + 1);
	);
    out
    )

computeDistanceFromOptimization = method(Options => {EquationType => "Curve"})
computeDistanceFromOptimization (Matrix,Matrix,ZZ) := o -> (pointMatrix,coeffMatrix,modelDegree) -> (
    xVal := pointMatrix_(0,0);
    yVal := pointMatrix_(0,1);
    x := local x;
    y := local y;
    R := RR[x,y];
    fBasis := returnModelBasis({x,y}, modelDegree);
    fxBasis := apply(fBasis, monom -> diff(x, monom));
    fyBasis := apply(fBasis, monom -> diff(y, monom));
    f := first flatten entries ((matrix{fBasis})*coeffMatrix) - 1;
    fx := first flatten entries ((matrix{fxBasis})*coeffMatrix);
    fy := first flatten entries ((matrix{fyBasis})*coeffMatrix);
    if o#EquationType == "Surface" then (
	z := local z;
	zVal := pointMatrix_(0,2);
	use RR[x,y,z];
	f = first flatten entries ((matrix{fBasis})*coeffMatrix) - z;
	fx = first flatten entries ((matrix{fxBasis})*coeffMatrix);
	fy = first flatten entries ((matrix{fyBasis})*coeffMatrix);
	fz := -1;
	);
    lambda := local lambda;
    if o#EquationType == "Curve" then (
	R1 := RR[x,y,lambda];
	(f,fx,fy) = toSequence apply({f,fx,fy}, expr -> sub(expr,R1));
	system := {f, 2*(x-xVal)-lambda*fx, 2*(y-yVal)-lambda*fy};
	(xSol,ySol,lambdaSol) := toSequence coordinates first solveSystem system;
	out := sqrt((xSol-xVal)^2 + (ySol-yVal)^2);
	) else if o#EquationType == "Surface" then (
	R2 := RR[x,y,z,lambda];
	(f,fx,fy,fz) = toSequence apply({f,fx,fy,fz}, expr -> sub(expr,R2));
	system = {f, 2*(x-xVal)-lambda*fx, 2*(y-yVal)-lambda*fy, 2*(z-zVal)-lambda*fz};
	zSol := null;
	(xSol,ySol,zSol,lambdaSol) = toSequence coordinates first solveSystem system;
	out = sqrt((xSol-xVal)^2 + (ySol-yVal)^2 + (zSol-zVal)^2);
	);
    out
    )

computeSquaredError = method(Options => {EquationType => "Curve"})
computeSquaredError (Matrix,Matrix,ZZ) := o -> (positionMatrix,coeffMatrix,modelDegree) -> (
    xVal := positionMatrix_(0,0);
    yVal := positionMatrix_(0,1);
    x := local x;
    y := local y;
    R := RR[x,y];
    f := returnModelBasis({x,y}, modelDegree);
    fM := matrix{apply(f,monom -> sub(monom,matrix{{xVal,yVal}}))};
    if o#EquationType == "Curve" then (
	out := (first flatten entries (fM*coeffMatrix)-1)^2;
	) else if o#EquationType == "Surface" then (
	out = (first flatten entries ((matrix{{1}}|fM)*coeffMatrix)-positionMatrix_(0,2))^2;
	);
    out
    )

computeAllApproxDistanceStats = method(Options => {EquationType => "Curve"})
computeAllApproxDistanceStats (Matrix,Matrix,ZZ) := o -> (positionMatrix,coeffMatrix,modelDegree) ->(
    approxDistanceList := apply(numRows positionMatrix, i -> computeApproxDistance(positionMatrix^{i}, coeffMatrix,modelDegree,EquationType => o#EquationType));
    residualList := apply(numRows positionMatrix, i -> computeSquaredError(positionMatrix^{i}, coeffMatrix,modelDegree,EquationType => o#EquationType));
    RMSE := sqrt((sum residualList)/#residualList);
    numDataPoints := #approxDistanceList;
    meanApproxDistance := (sum approxDistanceList) / numDataPoints;
    maxApproxDistance := max approxDistanceList;
    (meanApproxDistance, maxApproxDistance, RMSE)
    )

computeAllOptimizationDistanceStats = method(Options => {EquationType => "Curve"})
computeAllOptimizationDistanceStats (Matrix,Matrix,ZZ) := o -> (positionMatrix,coeffMatrix,modelDegree) ->(
    optimizationDistanceList := apply(numRows positionMatrix, i -> computeDistanceFromOptimization(positionMatrix^{i}, coeffMatrix,modelDegree,EquationType => o#EquationType));
    numDataPoints := #optimizationDistanceList;
    meanOptimizationDistance := (sum optimizationDistanceList) / numDataPoints;
    maxOptimizationDistance := max optimizationDistanceList;
    (meanOptimizationDistance, maxOptimizationDistance)
    )

peekOrbitEquationCoefficients = method()
peekOrbitEquationCoefficients (List,List) := (scaledCList,coeffList) -> (
    assert (#scaledCList > 1);
    n := floor(#scaledCList/2);
    << "-*" << endl;
    << "Printing the coefficients of the orbit equations for the two orbits in the middle:" << endl;
    scan(2, i -> (
	    << "--------------------------------------------------------" << endl;
	    << "Jacobi Constant: " << scaledCList_(n+i) << " => " <<  coeffList_(n+i) << endl;
	    )
	)
    << "--------------------------------------------------------" << endl;
    << "*-" << endl;
    )