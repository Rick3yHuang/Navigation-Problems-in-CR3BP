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
------------------------------------------------------------------------------
observedPointA1 = positionMatrixA^{pointIndicesA_0};
observedPointA2 = positionMatrixA^{pointIndicesA_1};
observedPointB1 = positionMatrixB^{pointIndicesB_0};
observedPointB2 = positionMatrixB^{pointIndicesB_1};
------------------------------------------------------------------------------
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
    zList = apply(#yList, i-> (
	    y = yList_i;
	    x = xList_i;
	    heightModel = heightModelList_(floor(i/2));
	    scaledCTrue = scaledCTrueList#(floor(i/2));
	    sub(heightModel,matrix{{x,y,scaledCTrue}})
	    ));
    (A1,A2,B1,B2) = toSequence apply(#yList,i -> matrix{{xList_i,yList_i,zList_i}});
    )
distAB1 = norm_2(A1 - B1);
distAB2 = norm_2(A2 - B2);
(SA1,SA2,SB1,SB2) = toSequence (norm_2 \ {A1,A2,B1,B2});
(dirA1,dirA2,dirB1,dirB2) = toSequence (normalize \ {A1,A2,B1,B2});
<< "------------------------------------------------------------------------------" << endl;
if fabricateFromPointsOnModels then (
	<< "-- The observed points are fabricated from the fitted models:" << endl;
	);
<< "-- The measurements provided are as follows: " << endl;
<< "---- direction of LOS to A1: " << dirA1 << endl;
<< "---- direction of LOS to A2: " << dirA2 << endl;
<< "---- direction of LOS to B1: " << dirB1 << endl;
<< "---- direction of LOS to B2: " << dirB2 << endl;
<< "---- distance between A1 and B1: " << distAB1 << endl;
<< "---- distance between A2 and B2: " << distAB2 << endl;
<< "------------------------------------------------------------------------------" << endl;
<< "-- True Solutions are:" << endl;
<< "---- C for A1 and A2: " << scaledCList#(orbitIDA-1) << endl;
<< "---- C for B1 and B2: " << scaledCList#(orbitIDB-1) << endl;
<< "---- parameters sA1, sA2 for A1, A2: " << SA1 << ", " << SA2 << endl;
<< "---- parameters sB1, sB2 for B1, B2: " << SB1 << ", " << SB2 << endl;
<< "------------------------------------------------------------------------------" << endl;
trueSols = {scaledCList#(orbitIDA-1),scaledCList#(orbitIDB-1),SA1,SA2,SB1,SB2};
