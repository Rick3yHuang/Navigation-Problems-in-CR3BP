-----------------------------------------------------------------------------
------------------------------ Read Data ------------------------------------
-----------------------------------------------------------------------------
if isLyapunov then (
    inputDirectory = lyapunovDataDirectory
    ) else if isHalo then (
    inputDirectory = haloDataDirectory
    ) else(
    error("Unsupported orbit type: " | orbitType)
    );

if isLyapunov then (
    dataIDList = toList(1..1037)
    ) else if isHalo then (
    dataIDList = toList(1..1147);
    );
    
CList = {};
positionMatrixList = {}
massRatioEarthMoon = 1.215058560962404e-2
elapsedTime scan(dataIDList, dataID -> (
	fileName = (inputDirectory | orbitType | "_Orbit_Data_ID_" | (toString dataID) | ".txt");
	(C, positionMatrix) = loadPositionDataForFitting(fileName,-massRatioEarthMoon,OrbitType=>orbitType);
	CList = append(CList,C);
	positionMatrixList = append(positionMatrixList,positionMatrix);
	));
<< "-*" << endl;
<< "   Loaded " << #CList << " data sets." << endl;
<< "   Jacobi Constant range: [" << min CList << ", " << max CList << "]." << endl;
<< "   Sample data sizes: (ID " << dataIDList_0 << ") "
<< numRows positionMatrixList_0 << " points; (ID " << dataIDList_(-1) << ") "
<< numRows positionMatrixList_(-1) << " points." << endl;
<< "*-" << endl;
-----------------------------------------------------------------------------
------------------------------ Prepare Data ---------------------------------
-----------------------------------------------------------------------------
-- Rescale CList to make it in between -1 and 1
scaledCList = rescaleCBetweenMinus1And1 CList;
-- concatenate position matrices with scaled C values
posAndScaledCList = apply(#dataIDList, i -> positionMatrixList_i | transpose matrix{toList((numRows(positionMatrixList_i)):scaledCList_i)});

exponentMatrix = findModelExponentMatrix modelDegree;

<< "-*" << endl;
<< "The exponent matrix is: " << endl;
<< exponentMatrix << endl;
<< "*-" << endl;

if orbitType == "Halo" then (
    -- Project to z = -x plane
    normalVector = transpose matrix{{1,0,1}}; -- w
    changeOfBasisTransformation = transpose matrix{{-1/sqrt(2),0,1/sqrt(2)},{0,1,0},{1/sqrt(2),0,1/sqrt(2)}}; -- uvw = R * xyz
    transformedPositionMatrixList = apply(positionMatrixList, positionMatrix ->  transpose (changeOfBasisTransformation * (transpose positionMatrix)));

    -- Split into uv positions and w positions
    projectedPositionMatrixList = apply(transformedPositionMatrixList, positionMatrix -> positionMatrix_{0..1});
    wPositionMatrixList = apply(transformedPositionMatrixList, positionMatrix -> positionMatrix_{2})
    );
