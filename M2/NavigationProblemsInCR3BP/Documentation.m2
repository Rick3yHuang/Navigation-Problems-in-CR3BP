doc ///
  Key
   NavigationProblemsInCR3BP
  Headline
     A Macaulay2 package to work with fitting algebraic equations for CR3BP orbits
  Usage
     needsPackage "NavigationProblemsInCR3BP"
  Description
   Text
     {\em NavigationProblemsInCR3BP} fits algebraic equations to orbits in the Circular Restricted Three Body Problem (CR3BP). It provides methods for loading position data, fitting cubic splines, solving for orbit coefficients, and evaluating the fitted equations.
  Subnodes
    loadPositionDataForFitting
    storePowersInCResults
    leastSquareSolution
    solveForOrbitCoefficients
    solveForCCoefficients
    constructMinimalProblems
    findDegree
///


doc ///
  Key
   loadPositionDataForFitting
   (loadPositionDataForFitting,String,RR)
   [loadPositionDataForFitting,OrbitType]
  Headline
     Load Jacobi Constant and Position Data for Fitting algebraic orbit curves and translate the x coordinates of the positions if needed
  Usage
    loadPositionDataForFitting(fileName,xTranslation)
  Inputs
    fileName: String
       the file name containing the position data
    xTranslation: RR
       a real number indicating the translation to be applied to the x-coordinates of the positions
    OrbitType => String
       a string indicating the type of orbit, either "Lyapunov" or "Halo". Default is "Lyapunov". If "Halo" is chosen, the z-coordinates of the positions will also be loaded.
  Outputs
    C: RR
       the Jacobi constant for the orbit
    positionMatrix: Matrix
       a matrix where each row represents the (x,y) or (x,y,z) coordinates of a point on the orbit after applying the x-translation
  Description
   Text
     This function loads the Jacobi constant and position data from a specified file. It applies a translation to the x-coordinates of the positions based on the provided xTranslation value. The function supports both "Lyapunov" and "Halo" orbit types, adjusting the output accordingly.
///

doc ///
  Key
   storePowersInCResults
   (storePowersInCResults,String,Matrix)
  Headline
    Store the powers used in fitting the Jacobi constant C into a specified results file
  Usage
    storePowersInCResults(fileName,CPolynomialCoeffMatrix)
  Inputs
    fileName: String
       the file name to store the powers, typically contains the min and max values for C
    CPolynomialCoeffMatrix: Matrix
       a matrix containing the coefficients of the polynomial fitted to the Jacobi constant C; each row corresponds to the coefficients for a specific power of C to generate the coefficients of the orbit equations
  Outputs
    
  Description
   Text
     This function stores the powers used in fitting the Jacobi constant C into a specified results file. The powers are derived from the provided coefficient matrix, where each row corresponds to a specific power of C. This information is essential for reconstructing the fitted orbit equations later.
///

doc ///
  Key
    leastSquareSolution
    (leastSquareSolution,Matrix,Matrix)
  Headline
     Solve for the least squares solution to an overdetermined system of equations
  Usage
    leastSquareSolution(designMatrix, constantTerm)
  Inputs
    designMatrix: Matrix
       the design matrix representing the left-hand side of the equations by evaluating basis functions at given data points
    constantTerm: Matrix
       the constant term matrix representing the right-hand side of the equations
  Outputs
    solution: Matrix
       the matrix containing the least squares solution to the system of equations by minimizing the residual sum of squares between the observed and predicted values
  Description
   Text
     This function computes the least squares solution to an overdetermined system of equations represented by the design matrix and constant term. It minimizes the residual sum of squares between the observed and predicted values, providing an optimal fit for the given data.
  SeeAlso
///

doc ///
  Key
    solveForOrbitCoefficients
    (solveForOrbitCoefficients,List,Matrix)
  Headline
     Solve for the orbit coefficients using least squares fitting for a list of position data
  Usage
    solveForOrbitCoefficients(positionMatrixList, exponentMatrix)
  Inputs
    positionMatrixList: Matrix
       a list of matrices where each matrix contains the position data for an orbit; each row represents the (x,y) or (x,y,z) coordinates of a point on the orbit
    exponentMatrix: Matrix
       a matrix where each row contains the exponents for x and y (and z if applicable) used to construct the design matrix for fitting the orbit equations
  Outputs
    coeffList: List
       a list of column matrices where each matrix contains the coefficients of the fitted orbit equations for the given exponent basis
  Description
   Text
     This function computes the orbit coefficients for a list of position data using least squares fitting. It constructs the design matrix based on the provided exponent matrix and solves for the coefficients that best fit the orbit equations to the given position data.
  SeeAlso
///

doc ///
  Key
    solveForCCoefficients
    (solveForCCoefficients,List,Matrix,ZZ)
    [solveForCCoefficients,UseSpecificBasis]
  Headline
     Solve for the coefficients of the powers of the Jacobi constant C polynomial using least squares fitting
  Usage
    solveForCCoefficients(scaledCList, allCoefficients, maxDegree)
  Inputs
    scaledCList: List
       a list of scaled Jacobi constant C values corresponding to the orbits for which the coefficients are to be fitted
    allCoefficients: Matrix
       a matrix where each column contains the coefficients of the fitted orbit equations for each orbit in the scaledCList
    maxDegree: ZZ
       an integer specifying the maximum degree of the polynomial to be fitted for the Jacobi constant C
    UseSpecificBasis => String
       an optional string specifying a particular basis to be used for fitting the C polynomial; if not provided, the default power basis will be used.
  Outputs
    CPolynomialCoeffMatrix: Matrix
       a matrix when multiplied by the power basis of C up to maxDegree gives the coefficients for the orbit equations
  Description
   Text
     This function computes the coefficients of the powers of the Jacobi constant C polynomial using least squares fitting. It constructs a design matrix based on the provided scaled C values and solves for the coefficients that best fit the orbit equations to the given data. An optional specific basis can be used for fitting if provided.
  SeeAlso
///

doc ///
  Key
    constructMinimalProblems
   (constructMinimalProblems,List,List,Ring,List)
  Headline
    Construct the Lyapunov minimal problems based on orbit and distance constraints
  Usage
    constructMinimalProblems(orbitConstraintsTuples, distanceConstraintsTuples, R)
  Inputs
    orbitConstraintsTuples: List
       a list of tuples where each tuple contains a point matrix, a C polynomial coefficient matrix, and an energy value for the orbit constraints
    distanceConstraintsTuples: List
       a list of tuples where each tuple contains two point matrices and a distance value for the distance constraints
    R: Ring
       the polynomial ring in which the constraints are constructed
    maxDegreeList: List
	   a list of integers specifying the maximum degrees for the jacobi constant C and modeling polynomials
  Outputs
    minimalProblems: List
       a list of polynomials defining the minimal problem in the specified polynomial ring R
  Description
   Text
     This function constructs the Lyapunov minimal problems based on the provided orbit and distance constraints. It generates a list of polynomials that define the minimal problem in the specified polynomial ring R, allowing for further analysis and solution of the system.
  SeeAlso
    findDegree
///

doc ///
  Key
   findDegree
   (findDegree,List,List,Ring,List)
   [findDegree,OrbitScenario]
  Headline
    Find the degree of the minimal problem in navigation based on measurements and C polynomial coefficient matrices
  Usage
    findDegree(measurements, CPolynomialCoeffMatrixList, R)
  Inputs
    measurements: List
      a list of measurements including direction vectors and distances used to define the navigation problem
    CPolynomialCoeffMatrixList: List
      a list of C polynomial coefficient matrices corresponding to the orbits involved in the navigation problem
    R: Ring
      the polynomial ring in which the minimal problem is constructed
    maxDegreeList: List
	   a list of integers specifying the maximum degrees for the jacobi constant C and modeling polynomials      
    OrbitScenario => String
      an optional method used to construct the minimal problem; defaults to constructMinimalProblems if not provided
  Outputs
    degree: ZZ
      an integer representing the degree of the minimal problem in navigation
  Description
   Text
     This function computes the degree of the minimal problem in navigation based on the provided measurements and C polynomial coefficient matrices. It constructs the minimal problem using the specified constructor method and determines its degree, which is essential for understanding the complexity of the navigation solution.
  SeeAlso
    constructMinimalProblems
///

doc ///
  Key
    OrbitType
  Headline
     Optional input to specify the type of orbit for data loading and processing
  Description
   Text
     The optional input {\em OrbitType} allows users to specify the type of orbit they are working with when loading position or velocity data. The accepted values are "Lyapunov" and "Halo". Depending on the chosen orbit type, the data loading functions will adjust their behavior, such as including z-coordinates for "Halo" orbits. The default value is "Lyapunov".
  SeeAlso
     loadPositionDataForFitting
///
    
doc ///
  Key
    UseSpecificBasis
  Headline
     Optional input to specify a particular basis for fitting the Jacobi constant C polynomial
  Description
   Text
     The optional input {\em UseSpecificBasis} allows users to specify a particular basis for fitting the Jacobi constant C polynomial when using the {\em solveForCCoefficients} function. By providing a string that indicates the desired basis, users can customize the fitting process beyond the default power basis. If this option is not provided, the function will use the standard power basis for fitting.
  SeeAlso
     solveForCCoefficients
///
