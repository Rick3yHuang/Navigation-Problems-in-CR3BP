------------------------------------------------------------------------------------------
------------------------------------ Data Loader -----------------------------------------
------------------------------------------------------------------------------------------
loadPositionDataForFitting = method(Options => {OrbitType => "Lyapunov"})
loadPositionDataForFitting (String,RR) :=  o -> (fileName,xTranslation) -> (
    rawData := lines get fileName;
    C := value rawData_0;
    stateData := drop(rawData,1);
    positionMatrix := matrix apply(stateData, s -> {(value s)_0 - xTranslation, (value s)_1});
    if o#OrbitType == "Halo" then (
	positionMatrix = matrix apply(stateData, s -> {(value s)_0 - xTranslation, (value s)_1, (value s)_2});
	);
    C, positionMatrix
    )

loadVelocityDataForFitting = method(Options => {OrbitType => "Lyapunov"})
loadVelocityDataForFitting String :=  o -> fileName -> (
    rawData := lines get fileName;
    velocityMatrix := matrix apply(rawData, s -> {(value s)_2, (value s)_3});
    if o#OrbitType == "Halo" then (
	velocityMatrix = matrix apply(rawData, s -> {(value s)_3, (value s)_4, (value s)_5});
	);
    velocityMatrix
    )

loadScaledCList = method()
loadScaledCList String := fileName -> (
    rawData := lines get fileName;
    CList := apply(#rawData, i -> (value rawData_i)_1);
    scaledCList := rescaleCBetweenMinus1And1 CList
    )

------------------------------------------------------------------------------------------
------------------------------------ C Rescaling -----------------------------------------
------------------------------------------------------------------------------------------
rescaleCBetweenMinus1And1 = method()
rescaleCBetweenMinus1And1 List := CList -> (
    Cmin := min CList;
    Cmax := max CList;
    apply(CList, C -> 2*(C - Cmin)/(Cmax - Cmin) - 1)
    )

------------------------------------------------------------------------------------------
------------------------- Coefficient Storage and Retrieval ------------------------------
------------------------------------------------------------------------------------------
storeAllCoefficients = method()
storeAllCoefficients (String,Matrix) := (fileName,allCoefficients) -> (
    outFile := openOut fileName;
    outFile << entries allCoefficients << endl;
    outFile << close;
    )

readAllCoefficients = method()
readAllCoefficients String := fileName -> (
    text := lines get fileName;
    matrix value first text
    )

storePowersInCResults = method()
storePowersInCResults (String,Matrix) := (fileName,CPolynomialCoeffMatrix) -> (
    outFile := openOut fileName;
    outFile << entries CPolynomialCoeffMatrix << endl;
    outFile << close;
    )

readPowersInCResults = method()
readPowersInCResults String := fileName -> (
    text := first lines get fileName;
    matrix value text
    )

------------------------------------------------------------------------------------------
--------------------------- Load data for minimal problems -------------------------------
------------------------------------------------------------------------------------------

loadEnergyIntervalsAndModels = method()
loadEnergyIntervalsAndModels (List,String,List,List) := (IDIntervals,orbitType,scaledCList,maxDegreeList) -> (
    scaledCIntervals := {};
    CPolynomialCoeffMatrixList := {};
    heightCPolynomialCoeffMatrixList := {};
    (jacobiConstantDegree,modelDegree) := toSequence maxDegreeList;
    scan(IDIntervals, IDInterval -> (
	    (minID,maxID) := toSequence IDInterval;
	    scaledCIntervals = append(scaledCIntervals, interval(scaledCList_(minID-1),scaledCList_(maxID-1)));
	    CPolynomialCoeffMatrixList = append(CPolynomialCoeffMatrixList, readPowersInCResults (savedModelsDirectory | "Jacobi-Constant-Degree-" | jacobiConstantDegree | "-Model-Degree-" | modelDegree | "/" | orbitType | "-C-Polynomial-Coeff-Output-" | minID | "-to-" | maxID | ".txt"));
	    if orbitType == "Halo" then (
		heightCPolynomialCoeffMatrixList = append(heightCPolynomialCoeffMatrixList, readPowersInCResults (savedModelsDirectory | "Jacobi-Constant-Degree-" | jacobiConstantDegree | "-Model-Degree-" | modelDegree | "/" | orbitType | "-C-Polynomial-Height-Coeff-Output-" | minID | "-to-" | maxID | ".txt"));
	    );
	));
    out := new MutableHashTable from {
        "scaledIntervals" => scaledCIntervals,
        "CPolynomialCoeffMatrixList" => CPolynomialCoeffMatrixList,
        "heightCPolynomialCoeffMatrixList" => null
    };
    if orbitType == "Halo" then out#"heightCPolynomialCoeffMatrixList" = heightCPolynomialCoeffMatrixList;
    out
    )
