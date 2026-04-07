-- Pick orbit ID 500
orbitID = 300;
massRatioEarthMoon = 1.215058560962404e-2
(CTrue,positionMatrix) = loadPositionDataForFitting(dataDirectory | orbitType | "-Data/" | orbitType | "_Orbit_Data_ID_" | (toString orbitID) | ".txt",-massRatioEarthMoon,OrbitType=>orbitType);
<< "------------------------------------------------------------------------------" << endl;
<< "-- The selected orbit ID is " << orbitID << endl;
<< "-- The Jacobi constant C for this orbit is " << CTrue << ", which is " << scaledCList#(orbitID-1) << " after rescaling" << endl;
-- Assume we know the xy-coordinates of the observed Point along the orbit:
observedPoint = positionMatrix^{100};
<< "--" << observedPoint << endl;
<< "------------------------------------------------------------------------------" << endl;
