-----------------------------------------------------------------------------
---------------- Fabricate an observed point from known data ----------------
-----------------------------------------------------------------------------
scaledCList = loadScaledCList (dataDirectory | orbitType | "-Data/ID_Jacobi_Constants.txt");
-----------------------------------------------------------------------------
--------------------- Load the saved coeffs ---------------------------------
-----------------------------------------------------------------------------
(scaledCIntervals,CPolynomialCoeffMatrixList) = loadEnergyIntervalsAndModels(IDIntervals,orbitType,scaledCList);
needs (minimalProblemDirectory | orbitType | "-Orbits/Lyapunov-One-Spacecraft-IOD-Fabrication.m2")
end

-----------------------------------------------------------------------------
---------- start M2 in "../../NavigationProblemsInCR3BP.m2" -------------
-----------------------------------------------------------------------------
restart
needsPackage "NavigationProblemsInCR3BP"
orbitType = "Lyapunov";		       -- Type of orbit to fit
IDIntervals = {{1,200},{200,400},{400,600},{600,800},{800,1037}};
needs (minimalProblemDirectory | orbitType | "-Orbits/Lyapunov-One-Spacecraft-IOD.m2")
solveForCWithAllModels(observedPoint,scaledCIntervals,CPolynomialCoeffMatrixList,CC[c])
