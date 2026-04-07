-----------------------------------------------------------------------------
--------------------- Load the saved coeffs ---------------------------------
-----------------------------------------------------------------------------
scaledCList = loadScaledCList (dataDirectory | orbitType | "-Data/ID_Jacobi_Constants.txt");
(scaledCIntervals,CPolynomialCoeffMatrixList,heightCPolynomialCoeffMatrixList) = loadEnergyIntervalsAndModels (IDIntervals,orbitType,scaledCList,{jacobiConstantDegree,modelDegree});

scaledCATrue = scaledCList#(orbitIDA-1);
scaledCBTrue = scaledCList#(orbitIDB-1);
scaledCTrueList = {scaledCATrue,scaledCBTrue};
massRatioEarthMoon = 1.215058560962404e-2
(CATrue,positionMatrixA) = loadPositionDataForFitting(dataDirectory | orbitType | "-Data/" | orbitType | "_Orbit_Data_ID_" | (toString orbitIDA) | ".txt",-massRatioEarthMoon,OrbitType=>orbitType);
(CBTrue,positionMatrixB) = loadPositionDataForFitting(dataDirectory | orbitType | "-Data/" | orbitType | "_Orbit_Data_ID_" | (toString orbitIDB) | ".txt",-massRatioEarthMoon,OrbitType=>orbitType);
<< "------------------------------------------------------------------------------" << endl;
<< "-- The selected orbit ID are " << orbitIDA << " and " << orbitIDB << endl;
<< "-- The Jacobi constant C for these two orbits are " << CATrue << " and " << CBTrue <<endl;
<< "-- which are " << scaledCATrue << " and " << scaledCBTrue << " after rescaling" << endl;
-- If pick directly from the data
observedPointA1 = positionMatrixA^{pointIndicesA_0};
observedPointA2 = positionMatrixA^{pointIndicesA_1};
observedPointB1 = positionMatrixB^{pointIndicesB_0};
observedPointB2 = positionMatrixB^{pointIndicesB_1};
-- If fabricate from the fitted models
if fabricateFromPointsOnModels then (
    modelA = modelPolynomial(RR[x,y,c],CPolynomialCoeffMatrixList_0,jacobiConstantDegree,modelDegree);
    modelB = modelPolynomial(RR[x,y,c],CPolynomialCoeffMatrixList_1,jacobiConstantDegree,modelDegree);
    modelList = {modelA,modelB};
    heightModelList = apply(2, i -> modelPolynomial(RR[u,v,c],heightCPolynomialCoeffMatrixList_i,jacobiConstantDegree,modelDegree,FittingHeight=>true));
    yList = flatten entries (0.1*random(RR^4,RR^1));
    xList = apply(#yList, i -> (
	    y = yList_i;
	    model = modelList_(floor(i/2));
	    scaledCTrue = scaledCTrueList#(floor(i/2));
	    eqn = sub(model,sub(matrix{{x,y,scaledCTrue}},RR[x]));
	    realPart first select(roots sub(model,sub(matrix{{x,y,scaledCTrue}},RR[x])), r->imaginaryPart r < 1e-6)
	    ));
    zList = apply(#yList, i -> (
	    y = yList_i;
	    x = xList_i;
	    heightModel = heightModelList_(floor(i/2));
	    scaledCTrue = scaledCTrueList#(floor(i/2));
	    sub(heightModel,matrix{{x,y,scaledCTrue}})
	    ));
    (A1,A2,B1,B2) = toSequence apply(#yList,i -> matrix{{xList_i,yList_i,zList_i}});
    )
M = matrix{flatten entries (0.1*random(RR^3,RR^1))};
distA1M = norm_2(A1 - M);
distB1M = norm_2(B1 - M);
distA1B1 = norm_2(A1 - B1);
distA2M = norm_2(A2 - M);
distB2M = norm_2(B2 - M);
distA2B2 = norm_2(A2 - B2);
<< "------------------------------------------------------------------------------" << endl;
if fabricateFromPointsOnModels then (
	<< "-- The observed points are fabricated from the fitted models:" << endl;
	);
<< "-- The measurements provided are as follows: " << endl;
<< "---- The position of the mothership: " << flatten entries M << endl;
<< "---- mothership-spacecraft distance of spacecraft A1: " << distA1M << endl;
<< "---- mothership-spacecraft distance of spacecraft B1: " << distB1M << endl;
<< "---- interspacecraft distance between spacecraft A1 and B1: " << distA1B1 << endl;
<< "---- mothership-spacecraft distance of spacecraft A2: " << distA2M << endl;
<< "---- mothership-spacecraft distance of spacecraft B2: " << distB2M
<< "---- interspacecraft distance between spacecraft A2 and B2: " << distA2B2 << endl;
<< "------------------------------------------------------------------------------" << endl;
<< "-- True Solutions are:" << endl;
<< "---- C for the orbit A: " << scaledCATrue << endl;
<< "---- C for the orbit B: " << scaledCBTrue << endl;
<< "---- The positions of the spacecraft A1: " << flatten entries A1 << endl;
<< "---- The positions of the spacecraft B1: " << flatten entries B1 << endl;
<< "---- The positions of the spacecraft A2: " << flatten entries A2 << endl;
<< "---- The positions of the spacecraft B2: " << flatten entries B2 << endl;
<< "------------------------------------------------------------------------------" << endl;
trueSols = {scaledCATrue,scaledCBTrue} | flatten entries A1 | flatten entries B1 | flatten entries A2 | flatten entries B2;
