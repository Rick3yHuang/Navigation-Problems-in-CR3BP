-*- coding: utf-8 -*-
newPackage("NavigationProblemsInCR3BP",
    AuxiliaryFiles => true, 
    Version => "0.1",
    Date => "April 24, 2025",
    Authors => {
	{Name => "Ruiqi (Rickey) Huang", Email => "rickey.huang747@gmail.com"}
	},
    Headline => "Minimal Navigation Problems in CR3BP Package",
    Keywords => {"Three Body Problem", "Circular Restricted Three Body Problem"},
    PackageImports => {"NumericalAlgebraicGeometry"},
    DebuggingMode => true
    )

export {
    ----------------------------
    -- Fitting Orbits ----------
    ----------------------------
    -- Data IO
    "loadPositionDataForFitting",
    "loadScaledCList",
    "storePowersInCResults",
    "readPowersInCResults",
    "storeAllCoefficients",
    "readAllCoefficients",
    -- Rescaling C
    "rescaleCBetweenMinus1And1",
    -- Fitting Methods
    "leastSquareSolution",
    "constructOrbitDesignMatrix",
    "constructCPolynomialDesignMatrix",
    "solveForOrbitCoefficients",
    "solveForHeightCoefficients",
    "solveForCCoefficients",
    "gradientDescentFit",
    "linearGradientDescentFit",
    -- Evaluation Methods
    "evaluateLeastSquaresFit",
    "modelPolynomial",
    "baseModelPolynomial",
    -- Minimal Problems
    "loadEnergyIntervalsAndModels",
    "findDegree",
    -- Optional inputs
    "UseSpecificBasis",
    "OrbitType",
    "FittingHeight",
    "MinimalProblemConstructor",
    "EquationType",
    "SolveByOptimization",
    "OrbitScenario",
    "DerivativesOfDistances",
    "MaxGDSteps",
    "GDLearningRate",
    "GDRelativeErrorTol",
    "InitialGuess",
    -- Helper methods
    "conditionNum",
    "listOfVectorsToMatrix",
    "normalize",
    "findEffectivePotential",
    "findEarthMoonDistance",
    "returnAllModelBasis",
    "findModelExponentMatrix",
    "returnModelBasis",
    ----------------------------
    -- Directory Shortcuts -----
    ----------------------------
    "repoDirectory",
    "baseDirectory",
    "dataDirectory",
    "lyapunovDataDirectory",
    "haloDataDirectory",
    "scratchDirectory",
    "savedModelsDirectory",
    "fittingExperimentDirectory",
    "minimalProblemDirectory"
    }
--------------------------------------------------------------------
----- Directories
--------------------------------------------------------------------
repoDirectory = currentDirectory() | "../"
baseDirectory = repoDirectory | "M2/NavigationProblemsInCR3BP/"
dataDirectory = repoDirectory | "Data/"
lyapunovDataDirectory = dataDirectory | "Lyapunov-Data/"
haloDataDirectory = dataDirectory | "Halo-Data/"
savedModelsDirectory = repoDirectory | "Saved-Models/"
orbitFittingDirectory = baseDirectory | "Orbit-Curve-Fitting/"
minimalProblemDirectory = baseDirectory | "Minimal-Navigation-Problems/"

--------------------------------------------------------------------
----- CODE
--------------------------------------------------------------------
load(baseDirectory | "Data-IO.m2")
load(baseDirectory | "Fitting-Methods.m2")
load(baseDirectory | "Evaluate-Fitting-Quality.m2")
load(baseDirectory | "Construct-Minimal-Problems.m2")
load(baseDirectory | "Helper-Methods.m2")
--------------------------------------------------------------------
----- TESTS
--------------------------------------------------------------------

--------------------------------------------------------------------
----- DOCUMENTATION
--------------------------------------------------------------------
beginDocumentation()
load(baseDirectory | "Documentation.m2")
end

--------------------------------------------------------------------
----- SCRATCH SPACE
--------------------------------------------------------------------
restart
uninstallPackage "NavigationProblemsInCR3BP"
installPackage "NavigationProblemsInCR3BP"
needsPackage "NavigationProblemsInCR3BP"
debug "NavigationProblemsInCR3BP"
check "NavigationProblemsInCR3BP"
