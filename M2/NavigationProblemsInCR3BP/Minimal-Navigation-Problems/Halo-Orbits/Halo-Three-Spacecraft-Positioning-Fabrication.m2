-----------------------------------------------------------------------------
--------------------- Load the saved coeffs ---------------------------------
-----------------------------------------------------------------------------
scaledCList = loadScaledCList (dataDirectory | orbitType | "-Data/ID_Jacobi_Constants.txt");
(scaledCIntervals,CPolynomialCoeffMatrixList,heightCPolynomialCoeffMatrixList) = loadEnergyIntervalsAndModels (IDIntervals,orbitType,scaledCList,{jacobiConstantDegree,modelDegree});

massRatioEarthMoon = 1.215058560962404e-2
(CATrue,positionMatrixA) = loadPositionDataForFitting(dataDirectory | orbitType | "-Data/" | orbitType | "_Orbit_Data_ID_" | (toString orbitIDA) | ".txt",-massRatioEarthMoon,OrbitType=>orbitType);
(CBTrue,positionMatrixB) = loadPositionDataForFitting(dataDirectory | orbitType | "-Data/" | orbitType | "_Orbit_Data_ID_" | (toString orbitIDB) | ".txt",-massRatioEarthMoon,OrbitType=>orbitType);
(CCTrue,positionMatrixC) = loadPositionDataForFitting(dataDirectory | orbitType | "-Data/" | orbitType | "_Orbit_Data_ID_" | (toString orbitIDC) | ".txt",-massRatioEarthMoon,OrbitType=>orbitType);
(scaledCATrue,scaledCBTrue,scaledCCTrue) = toSequence scaledCList_{orbitIDA-1, orbitIDB-1, orbitIDC-1};
scaledCTrueList = {scaledCATrue, scaledCBTrue, scaledCCTrue};
<< "------------------------------------------------------------------------------" << endl;
<< "-- The selected orbit ID are " << orbitIDA << ", " << orbitIDB << " and " << orbitIDC << endl;
<< "-- The Jacobi constant C for these two orbits are " << CATrue << ", " << CBTrue << " and " << CCTrue <<endl;
<< "-- which are " << scaledCATure << ", " << scaledCBTrue << " and " << scaledCCTrue << " after rescaling" << endl;

-- Assume we know the xy-coordinates of the observed Point along the orbit:
observedPointA = positionMatrixA^{pointIndexA};
observedPointB = positionMatrixB^{pointIndexB};
observedPointC = positionMatrixC^{pointIndexC};

if fabricateFromPointsOnModels then (
    modelList = apply(3, i -> modelPolynomial(RR[x,y,c],CPolynomialCoeffMatrixList_i,jacobiConstantDegree,modelDegree));
    heightModelList = apply(3, i -> modelPolynomial(RR[u,v,c],heightCPolynomialCoeffMatrixList_i,jacobiConstantDegree,modelDegree,FittingHeight=>true));
    yList = flatten entries (0.1*random(RR^3,RR^1));
    xList = apply(#yList, i -> (
	    y = yList_i;
	    model = modelList_i;
	    scaledCTrue = scaledCTrueList_i;
	    eqn = sub(model,sub(matrix{{x,y,scaledCTrue}},RR[x]));
	    realPart first select(roots sub(model,sub(matrix{{x,y,scaledCTrue}},RR[x])), r->imaginaryPart r < 1e-6)
	    ));
    zList = apply(#yList, i -> (
	    y = yList_i;
	    x = xList_i;
    	    heightModel = heightModelList_i;
	    scaledCTrue = scaledCTrueList_i;
	    sub(heightModel,matrix{{x,y,scaledCTrue}})
	    ));
    observedPointA = matrix{{xList_0,yList_0,zList_0}};
    observedPointB = matrix{{xList_1,yList_1,zList_1}};
    observedPointC = matrix{{xList_2,yList_2,zList_2}};
    )
distAB = norm_2(observedPointA - observedPointB);
distAC = norm_2(observedPointA - observedPointC);
distBC = norm_2(observedPointB - observedPointC);
<< "------------------------------------------------------------------------------" << endl;
if fabricateFromPointsOnModels then (
	<< "-- The observed points are fabricated from the fitted models:" << endl;
	);
<< "-- The measurements provided are as follows: " << endl;
<< "---- distance between A and B: " << distAB << endl;
<< "---- distance between A and C: " << distAC << endl;
<< "---- distance between B and C: " << distBC << endl;
<< "---- Jacobi constant for Spacecraft A: " << scaledCATrue << endl;
<< "---- Jacobi constant for Spacecraft B: " << scaledCBTrue << endl;
<< "---- Jacobi constant for Spacecraft C: " << scaledCCTrue << endl;
<< "------------------------------------------------------------------------------" << endl;
<< "-- True Solutions are:" << endl;
<< "-- Spacecraft A Position: " << observedPointA << endl;
<< "-- Spacecraft B Position: " << observedPointB << endl;
<< "-- Spacecraft C Position: " << observedPointC << endl;
<< "------------------------------------------------------------------------------" << endl;
trueSols = flatten entries observedPointA | flatten entries observedPointB | flatten entries observedPointC;
