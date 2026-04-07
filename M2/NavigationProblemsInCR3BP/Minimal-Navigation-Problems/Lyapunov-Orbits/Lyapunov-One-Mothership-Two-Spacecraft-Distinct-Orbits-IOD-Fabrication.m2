-----------------------------------------------------------------------------
--------------------- Load the saved coeffs ---------------------------------
-----------------------------------------------------------------------------
scaledCList = loadScaledCList (dataDirectory | orbitType | "-Data/ID_Jacobi_Constants.txt");
(scaledCIntervals,CPolynomialCoeffMatrixList) = loadEnergyIntervalsAndModels (IDIntervals,orbitType,scaledCList,{jacobiConstantDegree,modelDegree});

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
A = (positionMatrixA^{pointIndexA})_{0..1};
B = (positionMatrixB^{pointIndexB})_{0..1};
-- If fabricate from the fitted models
if fabricateFromPointsOnModels then (
    modelA = modelPolynomial(RR[x,y,c],CPolynomialCoeffMatrixList_0,jacobiConstantDegree,modelDegree);
    modelB = modelPolynomial(RR[x,y,c],CPolynomialCoeffMatrixList_1,jacobiConstantDegree,modelDegree);
    modelList = {modelA,modelB};
    yList = flatten entries (0.1*random(RR^2,RR^1));
    xList = apply(#yList, i -> (
	    y = yList_i;
	    model = modelList_i;
	    scaledCTrue = scaledCTrueList#i;
	    eqn = sub(model,sub(matrix{{x,y,scaledCTrue}},RR[x]));
	    realPart first select(roots sub(model,sub(matrix{{x,y,scaledCTrue}},RR[x])), r->imaginaryPart r < 1e-6)
	    ));
    (A,B) = toSequence apply(#yList,i -> matrix{{xList_i,yList_i}});
    vyList = flatten entries (0.1*random(RR^2,RR^1));
    vxList = apply(#vyList, i -> (2*findEffectivePotential(transpose ({A,B}_i),massRatioEarthMoon,OrbitScenario => orbitType)-scaledCTrueList_i-vyList_i));
    (vA,vB) = toSequence apply(#vyList,i -> matrix{{vxList_i,vyList_i}});
    )
M = matrix{flatten entries (0.1*random(RR^2,RR^1))};
earthMoonDist = flatten apply({A,B}, point -> findEarthMoonDistance(point,massRatioEarthMoon,OrbitScenario => orbitType));
vM = matrix{{0,0}};
distAM = norm_2(A - M);
distBM = norm_2(B - M);
distAB = norm_2(A - B);
distDerivativeAM = first flatten entries ((A - M)*transpose(vA-vM))/distAM;
distDerivativeBM = first flatten entries ((B - M)*transpose(vB-vM))/distBM;
distDerivativeAB = first flatten entries ((A - B)*transpose(vA-vB))/distAB;
<< "------------------------------------------------------------------------------" << endl;
if fabricateFromPointsOnModels then (
	<< "-- The observed points are fabricated from the fitted models:" << endl;
	);
<< "-- The measurements provided are as follows: " << endl;
<< "---- The position of the mothership: " << flatten entries M << endl;
<< "---- mothership-spacecraft distance of spacecraft A: " << distAM << endl;
<< "---- mothership-spacecraft distance of spacecraft B: " << distBM << endl;
<< "---- interspacecraft distance between spacecraft A and B: " << distAB << endl;
<< "---- the first derivative of mothership-spacecraft distance of spacecraft A: " << distDerivativeAM << endl;
<< "---- the first derivative of mothership-spacecraft distance of spacecraft B: " << distDerivativeBM << endl;
<< "---- the first derivative of interspacecraft distance between spacecraft A and B:" << distDerivativeAB << endl;
<< "------------------------------------------------------------------------------" << endl;
<< "-- True Solutions are:" << endl;
<< "---- C for the orbit A: " << scaledCATrue << endl;
<< "---- C for the orbit B: " << scaledCBTrue << endl;
<< "---- The positions of the spacecraft A: " << flatten entries A << endl;
<< "---- The positions of the spacecraft B: " << flatten entries B << endl;
<< "---- The velocities of the spacecraft A: " << flatten entries vA << endl;
<< "---- The velocities of the spacecraft B: " << flatten entries vB << endl;
<< "------------------------------------------------------------------------------" << endl;
trueSols = {scaledCATrue,scaledCBTrue} | flatten entries A | flatten entries vA | flatten entries B | flatten entries vB | earthMoonDist
