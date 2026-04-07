-----------------------------------------------------------------------------
--------------------- Load the saved coeffs ---------------------------------
-----------------------------------------------------------------------------
scaledCList = loadScaledCList (dataDirectory | orbitType | "-Data/ID_Jacobi_Constants.txt");
(scaledCIntervals,CPolynomialCoeffMatrixList,heightCPolynomialCoeffMatrixList) = loadEnergyIntervalsAndModels(IDIntervals,orbitType,scaledCList,{jacobiConstantDegree,modelDegree});

massRatioEarthMoon = 1.215058560962404e-2
scaledCTrue = scaledCList#(orbitID-1);
(CTrue,positionMatrix) = loadPositionDataForFitting(dataDirectory | orbitType | "-Data/" | orbitType | "_Orbit_Data_ID_" | (toString orbitID) | ".txt",-massRatioEarthMoon,OrbitType=>orbitType);
<< "------------------------------------------------------------------------------" << endl;
<< "-- The selected orbit ID for the spacecraft is " << orbitID << endl;
<< "-- The Jacobi constant C for this orbit is " << CTrue << endl;
<< "-- which is " << scaledCTrue << " after rescaling" << endl;
-- Assume we know the xy-coordinates of the observed Point along the orbit:
changeOfBasisTransformation = transpose matrix{{-1/sqrt(2),0,1/sqrt(2)},{0,1,0},{1/sqrt(2),0,1/sqrt(2)}}
transformedObservedPointList = apply(#pointIndexList, i -> (
	pointIndex = pointIndexList_i;
	transpose (changeOfBasisTransformation*(transpose (positionMatrix^{pointIndex})))
	));
if fabricateFromPointsOnModels then (
    model = modelPolynomial(RR[u,v,c],CPolynomialCoeffMatrixList_0,jacobiConstantDegree,modelDegree);
    heightModel = modelPolynomial(RR[u,v,c],heightCPolynomialCoeffMatrixList_0,jacobiConstantDegree,modelDegree,FittingHeight=>true);
    vList = flatten entries (0.1*random(RR^2,RR^1));
    uList = apply(#vList, i -> (
	    v := vList_i;
	    vModel := sub(model,sub(matrix{{u, v, scaledCTrue}},RR[u]));
	    rt := roots vModel;
	    realRt := select(rt, r-> abs (imaginaryPart r) < 1e-6);
	    out := realPart first realRt;
	    out
	    ));
    wList = apply(#vList, i -> (
	    u := uList_i;
	    v := vList_i;
	    sub(heightModel,sub(matrix{{u, v, scaledCTrue}},RR))
	    ));
    transformedObservedPointA = matrix{{uList_0,vList_0,wList_0}};
    transformedObservedPointB = matrix{{uList_1,vList_1,wList_1}};
    )
M = matrix{flatten entries (0.1*random(RR^3,RR^1))};
distAB = norm_2(transformedObservedPointA - transformedObservedPointB);
distAM = norm_2(transformedObservedPointA - M);
distBM = norm_2(transformedObservedPointB - M);
<< "------------------------------------------------------------------------------" << endl;
if fabricateFromPointsOnModels then (
    << "-- The observed points are fabricated from the fitted models:" << endl;
    );
<< "-- The measurements provided are as follows: " << endl;
<< "---- The position of the mothership: " << flatten entries M << endl;
<< "---- mothership-spacecraft distance of spacecraft A: " << distAM << endl;
<< "---- mothership-spacecraft distance of spacecraft B: " << distBM << endl;
<< "---- interspacecraft distance between spacecraft A and B: " << distAB << endl;
<< "------------------------------------------------------------------------------" << endl;
<< "-- True Solutions are:" << endl;
<< "---- C for the orbit including both spacecraft: " << scaledCTrue << endl;
<< "---- The positions of the spacecraft A (transformed): " << flatten entries transformedObservedPointA << endl;
<< "---- The positions of the spacecraft B (transformed): " << flatten entries transformedObservedPointB << endl;
<< "------------------------------------------------------------------------------" << endl;
trueSols = {scaledCTrue} | flatten entries transformedObservedPointA | flatten entries transformedObservedPointB;
