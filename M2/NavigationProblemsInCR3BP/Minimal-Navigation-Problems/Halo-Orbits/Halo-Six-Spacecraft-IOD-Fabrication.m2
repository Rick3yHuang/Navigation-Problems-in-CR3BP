-----------------------------------------------------------------------------
--------------------- Load the saved coeffs ---------------------------------
-----------------------------------------------------------------------------
scaledCList = loadScaledCList (dataDirectory | orbitType | "-Data/ID_Jacobi_Constants.txt");
(scaledCIntervals,CPolynomialCoeffMatrixList,heightCPolynomialCoeffMatrixList) = loadEnergyIntervalsAndModels(IDIntervals,orbitType,scaledCList,{jacobiConstantDegree,modelDegree});

massRatioEarthMoon = 1.215058560962404e-2

(CATrue,positionMatrixA) = loadPositionDataForFitting(dataDirectory | orbitType | "-Data/" | orbitType | "_Orbit_Data_ID_" | (toString orbitIDA) | ".txt",-massRatioEarthMoon,OrbitType=>orbitType);
(CBTrue,positionMatrixB) = loadPositionDataForFitting(dataDirectory | orbitType | "-Data/" | orbitType | "_Orbit_Data_ID_" | (toString orbitIDB) | ".txt",-massRatioEarthMoon,OrbitType=>orbitType);
(CCTrue,positionMatrixC) = loadPositionDataForFitting(dataDirectory | orbitType | "-Data/" | orbitType | "_Orbit_Data_ID_" | (toString orbitIDC) | ".txt",-massRatioEarthMoon,OrbitType=>orbitType);
(CDTrue,positionMatrixD) = loadPositionDataForFitting(dataDirectory | orbitType | "-Data/" | orbitType | "_Orbit_Data_ID_" | (toString orbitIDD) | ".txt",-massRatioEarthMoon,OrbitType=>orbitType);
(CETrue,positionMatrixE) = loadPositionDataForFitting(dataDirectory | orbitType | "-Data/" | orbitType | "_Orbit_Data_ID_" | (toString orbitIDE) | ".txt",-massRatioEarthMoon,OrbitType=>orbitType);
(CFTrue,positionMatrixF) = loadPositionDataForFitting(dataDirectory | orbitType | "-Data/" | orbitType | "_Orbit_Data_ID_" | (toString orbitIDF) | ".txt",-massRatioEarthMoon,OrbitType=>orbitType);
(scaledCATrue,scaledCBTrue,scaledCCTrue,scaledCDTrue,scaledCETrue,scaledCFTrue) = toSequence scaledCList_{orbitIDA-1, orbitIDB-1, orbitIDC-1,orbitIDD-1, orbitIDE-1, orbitIDF-1};
CTrueList = {CATrue,CBTrue,CCTrue,CDTrue,CETrue,CFTrue};
scaledCTrueList = {scaledCATrue,scaledCBTrue,scaledCCTrue,scaledCDTrue,scaledCETrue,scaledCFTrue};
positionMatrixList = {positionMatrixA,positionMatrixB,positionMatrixC,positionMatrixD,positionMatrixE,positionMatrixF};
pointIndexList = {pointIndexA,pointIndexB,pointIndexC,pointIndexD,pointIndexE,pointIndexF};
<< "------------------------------------------------------------------------------" << endl;
<< "-- The selected orbit ID are " << orbitIDA << ", " << orbitIDB << ", " << orbitIDC << ", " << orbitIDD << ", " << orbitIDE << " and " << orbitIDF << endl;
<< "-- The Jacobi constant C for these six orbits are " << CATrue << ", " << CBTrue << ", " << CCTrue << ", " << CDTrue << ", " << CETrue << " and " << CFTrue <<endl;
<< "-- which are " << scaledCATrue << ", " << scaledCBTrue << ", " << scaledCCTrue << ", " << scaledCDTrue << ", " << scaledCETrue << " and " << scaledCFTrue << " after rescaling" << endl;
-- Assume we know the xy-coordinates of the observed Point along the orbit:
changeOfBasisTransformation = transpose matrix{{-1/sqrt(2),0,1/sqrt(2)},{0,1,0},{1/sqrt(2),0,1/sqrt(2)}}
transformedObservedPointList = apply(#pointIndexList, i -> (
	positionMatrix = positionMatrixList_i;
	pointIndex = pointIndexList_i;
	transpose (changeOfBasisTransformation*(transpose (positionMatrix^{pointIndex})))
	));
if fabricateFromPointsOnModels then (
    modelList = apply(6, i -> modelPolynomial(RR[u,v,c],CPolynomialCoeffMatrixList_i,jacobiConstantDegree,modelDegree));
    heightModelList = apply(6, i -> modelPolynomial(RR[u,v,c],heightCPolynomialCoeffMatrixList_i,jacobiConstantDegree,modelDegree,FittingHeight=>true));
    vList = flatten entries (0.1*random(RR^6,RR^1));
    uList = apply(#modelList, i -> (
	    model := modelList_i;
	    v := vList_i;
	    realPart first select(roots sub(model,sub(matrix{{u, v, scaledCTrueList_i}},RR[u])), r->imaginaryPart r < 1e-6)
	    ));
    wList = apply(6, i -> (
	    heightModel := heightModelList_i;
	    u := uList_i;
	    v := vList_i;
	    sub(heightModel,sub(matrix{{u, v, scaledCTrueList_i}},RR))
	    ));
    transformedObservedPointA = matrix{{uList_0,vList_0,wList_0}};
    transformedObservedPointB = matrix{{uList_1,vList_1,wList_1}};
    transformedObservedPointC = matrix{{uList_2,vList_2,wList_2}};
    transformedObservedPointD = matrix{{uList_3,vList_3,wList_3}};
    transformedObservedPointE = matrix{{uList_4,vList_4,wList_4}};
    transformedObservedPointF = matrix{{uList_5,vList_5,wList_5}};
    )
distAB = norm_2(transformedObservedPointA - transformedObservedPointB);
distAC = norm_2(transformedObservedPointA - transformedObservedPointC);
distBC = norm_2(transformedObservedPointB - transformedObservedPointC);
distAD = norm_2(transformedObservedPointA - transformedObservedPointD);
distBC = norm_2(transformedObservedPointB - transformedObservedPointC);
distBD = norm_2(transformedObservedPointB - transformedObservedPointD);
distCD = norm_2(transformedObservedPointC - transformedObservedPointD);

distAE = norm_2(transformedObservedPointA - transformedObservedPointE);
distBE = norm_2(transformedObservedPointB - transformedObservedPointE);
distCE = norm_2(transformedObservedPointC - transformedObservedPointE);

distAF = norm_2(transformedObservedPointA - transformedObservedPointF);
distBF = norm_2(transformedObservedPointB - transformedObservedPointF);
distCF = norm_2(transformedObservedPointC - transformedObservedPointF);
<< "------------------------------------------------------------------------------" << endl;
if fabricateFromPointsOnModels then (
	<< "-- The observed points are fabricated from the fitted models:" << endl;
	);
<< "-- The measurements provided are as follows: " << endl;
<< "---- distance between A and B: " << distAB << endl;
<< "---- distance between A and C: " << distAC << endl;
<< "---- distance between A and D: " << distAD << endl;
<< "---- distance between B and C: " << distBC << endl;
<< "---- distance between B and D: " << distBD << endl;
<< "---- distance between C and D: " << distCD << endl;

<< "---- distance between A and E: " << distAE << endl;
<< "---- distance between B and E: " << distBE << endl;
<< "---- distance between C and E: " << distCE << endl;

<< "---- distance between A and F: " << distAF << endl;
<< "---- distance between B and F: " << distBF << endl;
<< "---- distance between C and F: " << distCF << endl;
<< "------------------------------------------------------------------------------" << endl;
<< "-- True Solutions are:" << endl;
<< "-- orrbit A scaled Jacobi constant: " << scaledCATrue << endl;
<< "-- orrbit B scaled Jacobi constant: " << scaledCBTrue << endl;
<< "-- orrbit C scaled Jacobi constant: " << scaledCCTrue << endl;
<< "-- orrbit D scaled Jacobi constant: " << scaledCDTrue << endl;
<< "-- orrbit E scaled Jacobi constant: " << scaledCETrue << endl;
<< "-- orrbit F scaled Jacobi constant: " << scaledCFTrue << endl;
<< "-- Spacecraft A Position (transformed): " << transformedObservedPointA << endl;
<< "-- Spacecraft B Position (transformed): " << transformedObservedPointB << endl;
<< "-- Spacecraft C Position (transformed): " << transformedObservedPointC << endl;
<< "-- Spacecraft D Position (transformed): " << transformedObservedPointD << endl;
<< "-- Spacecraft E Position (transformed): " << transformedObservedPointE << endl;
<< "-- Spacecraft F Position (transformed): " << transformedObservedPointF << endl;
<< "------------------------------------------------------------------------------" << endl;
trueSols = {scaledCATrue,scaledCBTrue,scaledCCTrue,scaledCDTrue,scaledCETrue,scaledCFTrue} |
flatten entries transformedObservedPointA | flatten entries transformedObservedPointB |
flatten entries transformedObservedPointC | flatten entries transformedObservedPointD |
flatten entries transformedObservedPointE | flatten entries transformedObservedPointF;

