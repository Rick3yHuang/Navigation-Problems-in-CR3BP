-----------------------------------------------------------------------------
--------------------- Load the saved coeffs ---------------------------------
-----------------------------------------------------------------------------
scaledCList = loadScaledCList (dataDirectory | orbitType | "-Data/ID_Jacobi_Constants.txt");
(scaledCIntervals,CPolynomialCoeffMatrixList,heightCPolynomialCoeffMatrixList) = loadEnergyIntervalsAndModels (IDIntervals,orbitType,scaledCList,{jacobiConstantDegree,modelDegree});

scaledCTrue = scaledCList#(orbitID-1);
massRatioEarthMoon = 1.215058560962404e-2
(CTrue,positionMatrix) = loadPositionDataForFitting(dataDirectory | orbitType | "-Data/" | orbitType | "_Orbit_Data_ID_" | (toString orbitID) | ".txt",-massRatioEarthMoon,OrbitType=>orbitType);
<< "------------------------------------------------------------------------------" << endl;
<< "-- The selected orbit ID for the spacecraft is " << orbitID << endl;
<< "-- The Jacobi constant C for this orbit is " << CTrue << endl;
<< "-- which is " << scaledCTrue << "after rescaling" << endl;
-- If pick directly from the data
observedPoint = positionMatrix^{pointIndices_0};
-- If fabricate from the fitted models
if fabricateFromPointsOnModels then (
    model = modelPolynomial(RR[x,y,c],CPolynomialCoeffMatrixList_0,jacobiConstantDegree,modelDegree);
    heightModel = modelPolynomial(RR[u,v,c],heightCPolynomialCoeffMatrixList_0,jacobiConstantDegree,modelDegree,FittingHeight=>true);
    Ay = first flatten entries (0.1*random(RR^1,RR^1));
    Ax = realPart first select(roots sub(model,sub(matrix{{x,y,scaledCTrue}},RR[x])), r->imaginaryPart r < 1e-6);
    Az = sub(heightModel,matrix{{Ax,Ay,scaledCTrue}});
    A = matrix{{Ax,Ay,Az}};
    )
mothershipFabrication = flatten entries (0.1*random(RR^6,RR^1))
M1 = matrix{mothershipFabrication_{0,1,2}};
M2 = matrix{mothershipFabrication_{3,4,5}};
distAM1 = norm_2(A - M1);
distAM2 = norm_2(A - M2);
<< "------------------------------------------------------------------------------" << endl;
if fabricateFromPointsOnModels then (
	<< "-- The observed points are fabricated from the fitted models:" << endl;
	);
<< "-- The measurements provided are as follows: " << endl;
<< "---- The position of the mothership 1: " << flatten entries M1 << endl;
<< "---- The position of the mothership 2: " << flatten entries M2 << endl;
<< "---- mothership-spacecraft distance of mothership 1: " << distAM1 << endl;
<< "---- mothership-spacecraft distance of mothership 2: " << distAM2 << endl;
<< "------------------------------------------------------------------------------" << endl;
<< "-- True Solutions are:" << endl;
<< "---- C for the orbit including the spacecraft: " << scaledCTrue << endl;
<< "---- The positions of the spacecraft: " << flatten entries A << endl;
<< "------------------------------------------------------------------------------" << endl;
trueSols = {scaledCTrue} | flatten entries A;
