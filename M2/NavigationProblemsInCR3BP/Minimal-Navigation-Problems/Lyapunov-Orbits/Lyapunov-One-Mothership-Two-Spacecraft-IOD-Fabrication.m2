-----------------------------------------------------------------------------
--------------------- Load the saved coeffs ---------------------------------
-----------------------------------------------------------------------------
scaledCList = loadScaledCList (dataDirectory | orbitType | "-Data/ID_Jacobi_Constants.txt");
(scaledCIntervals,CPolynomialCoeffMatrixList) = loadEnergyIntervalsAndModels (IDIntervals,orbitType,scaledCList,{jacobiConstantDegree,modelDegree});

scaledCTrue = scaledCList#(orbitID-1);
massRatioEarthMoon = 1.215058560962404e-2
(CTrue,positionMatrix) = loadPositionDataForFitting(dataDirectory | orbitType | "-Data/" | orbitType | "_Orbit_Data_ID_" | (toString orbitID) | ".txt",-massRatioEarthMoon,OrbitType=>orbitType);
<< "------------------------------------------------------------------------------" << endl;
<< "-- The selected orbit ID for the spacecraft is " << orbitID << endl;
<< "-- The Jacobi constant C for this orbit is " << CTrue << endl;
<< "-- which is " << scaledCTrue << " after rescaling" << endl;
-- If pick directly from the data
observedPointA = positionMatrix^{pointIndexA};
observedPointB = positionMatrix^{pointIndexB};
-- If fabricate from the fitted models
if fabricateFromPointsOnModels then (
    model = modelPolynomial(RR[x,y,c],CPolynomialCoeffMatrixList_0,jacobiConstantDegree,modelDegree);
    (Ay,By) = toSequence flatten entries (0.1*random(RR^2,RR^1));
    Ax = realPart first select(roots sub(model,sub(matrix{{x,Ay,scaledCTrue}},RR[x])), r->abs (imaginaryPart r) < 1e-6);
    Bx = realPart first select(roots sub(model,sub(matrix{{x,By,scaledCTrue}},RR[x])), r->abs (imaginaryPart r) < 1e-6);
    A = matrix{{Ax,Ay}};
    B = matrix{{Bx,By}};
    if fabricateOneExtraPoint then (
	Py = 0.1*random(RR);
	Px = realPart first select(roots sub(model,sub(matrix{{x,Py,scaledCTrue}},RR[x])), r->abs (imaginaryPart r) < 1e-6);
	P = matrix{{Px,Py}};
	);
    )
M = matrix{flatten entries (0.1*random(RR^2,RR^1))};
distAM = norm_2(A - M);
distBM = norm_2(B - M);
distAB = norm_2(A - B);
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
<< "---- The positions of the spacecraft A: " << flatten entries A << endl;
<< "---- The positions of the spacecraft B: " << flatten entries B << endl;
<< "------------------------------------------------------------------------------" << endl;
trueSols = {scaledCTrue} | flatten entries A | flatten entries B;
