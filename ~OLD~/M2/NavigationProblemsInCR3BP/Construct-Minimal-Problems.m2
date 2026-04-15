---------------------------------------------------------------------------------------------
-- Method that constructs the minimal problems:
-- Lyapunov: unknowns: x,y
-- Halo: unknowns: u,v,w
-- Given: orbitConstraintsTuples = { (pointA,CPolynomialCoeffMatrixA,energyA), ... }
--        distanceConstraintsTuples = { (pointA,pointB,distAB), ... }
--        R is the polynomial ring where the contraints are constructed
-- Returns: List of polynomials defining the minimal problem in R
---------------------------------------------------------------------------------------------
constructMinimalProblems = method(Options => {DerivativesOfDistances => false, OrbitScenario => "Lyapunov"})
constructMinimalProblems (List,List,Ring,List) := o -> (orbitConstraintsTuples,distanceConstraintsTuples,R,maxDegreeList) -> (
    (jacobiConstantDegree,modelDegree) := toSequence maxDegreeList;
    orbitConstraints := apply(orbitConstraintsTuples, orbitConstraintsTuple -> (
	    if not o#DerivativesOfDistances then (point,CPolynomialCoeffMatrix,energy) := toSequence orbitConstraintsTuple_{0..2};
	    if o#DerivativesOfDistances then (point,CPolynomialCoeffMatrix,energy) = toSequence orbitConstraintsTuple_{0,2,3};
	    powerBasis := transpose matrix{apply(jacobiConstantDegree+1, i -> energy^i)};
	    orbitCoefficients := CPolynomialCoeffMatrix*powerBasis;
	    x := local x;
	    y := local y;
	    S := R[x,y];
	    modelBasis := returnModelBasis({x,y},modelDegree);
	    if #orbitConstraintsTuple == 3 then (
		out := (first flatten entries(sub(matrix{modelBasis},transpose point^{0..1})*sub(orbitCoefficients,R))) - 1;
		) else if #orbitConstraintsTuple == 4 then (
		heightCPolynomialCoeffMatrix := orbitConstraintsTuple_3;
		heightOrbitCoefficients := heightCPolynomialCoeffMatrix*powerBasis;
		z := local z;
		S = R[x,y,z];
		out = {};
		out = append(out, (first flatten entries(sub(matrix{modelBasis},transpose point^{0..1})*sub(orbitCoefficients,R))) - 1);
		out = append(out,first flatten entries((matrix{{1}}|sub(sub(matrix{modelBasis},S),transpose point))*sub(heightOrbitCoefficients,R))-point_(2,0));
		out
		);
	    out
	    )
	);
    energyConstraints := apply(orbitConstraintsTuples, orbitConstraintsTuple -> (
	    if o#DerivativesOfDistances then (
		(point,velocity) := toSequence orbitConstraintsTuple_{0..1};
		(energy,mu,sqrtDist) := toSequence orbitConstraintsTuple_{-3..-1};
		(findEffectivePotentialConstraints((point||velocity),sqrtDist,{energy,mu},OrbitScenario => o#OrbitScenario)) | (findEarthMoonDistanceConstraints(point,sqrtDist,mu,OrbitScenario => o#OrbitScenario))
		) else {}
	    )
	);
    if #(flatten distanceConstraintsTuples) != 0 then (
	distanceConstraints := apply(distanceConstraintsTuples, distanceConstraintsTuple -> (
		if not o#DerivativesOfDistances then (
		    (pointA,pointB,distAB) := toSequence distanceConstraintsTuple_{0..2};
		    ) else if o#DerivativesOfDistances then (
		    (velocityA,velocityB,distDerivativeAB) := (0,0,0);
		    (pointA,pointB,velocityA,velocityB,distAB,distDerivativeAB) = toSequence distanceConstraintsTuple;
		    );
		if #distanceConstraintsTuple == 3 then (
		    out := (sum apply(flatten entries (pointA - pointB), x -> x^2)) - distAB^2;
		    ) else if #distanceConstraintsTuple == 5 then (
		    (cosTheta,sinTheta) := toSequence distanceConstraintsTuple_{3..4};
		    out = {cosTheta^2 + sinTheta^2 - 1};
		    );
		if o#DerivativesOfDistances then (
		    out = {out} | {first flatten entries ((transpose (pointA - pointB))*(velocityA - velocityB)) - distAB*distDerivativeAB};
		    );
		out
		)
	    );
	) else distanceConstraints = {};
    (flatten orbitConstraints) | (flatten energyConstraints) | (flatten distanceConstraints)
    )

findEarthMoonDistanceConstraints = method(Options => {OrbitScenario => "Lyapunov"})
findEarthMoonDistanceConstraints (Matrix,Matrix,RR) := o -> (pos,sqrtDist,mu) -> (
    (r1,r2) := toSequence flatten entries sqrtDist;
    if o#OrbitScenario == "Lyapunov" then (
	(x,y) := toSequence flatten entries pos;
	r1Constraint := r1^2 - (x + mu)^2 - y^2;
	r2Constraint := r2^2 - (x - (1 - mu))^2 - y^2;
	) else if o#OrbitScenario == "Halo" then (
	z := 0;
	(x,y,z) = toSequence flatten entries pos;
	r1Constraint = r1^2 - (x + mu)^2 - y^2 - z^2;
	r2Constraint = r2^2 - (x - (1 - mu))^2 - y^2 - z^2;
	);
    {r1Constraint,r2Constraint}
    )

findEffectivePotentialConstraints = method(Options => {OrbitScenario => "Lyapunov"})
findEffectivePotentialConstraints (Matrix,Matrix,List) := o -> (stateMatrix,sqrtDist,energyAndMu) -> (
    (energy,mu) := toSequence energyAndMu;
    (r1,r2) := toSequence flatten entries sqrtDist;
    if o#OrbitScenario == "Lyapunov" then (
	(x,y,vx,vy) := toSequence flatten entries stateMatrix;
	out := r1*r2*(vx+vy) - r1*r2*(x^2+y^2) - 2*r2*(1 - mu) - 2*r1*mu + r1*r2*energy;
	) else if o#OrbitScenario == "Halo" then (
	z := 0;
	vz := 0;
	(x,y,z,vx,vy,vz) = toSequence flatten entries stateMatrix;
	out = r1*r2*(vx+vy+vz) - r1*r2*(x^2+y^2+z^2) - 2*r2*(1 - mu) - 2*r1*mu + r1*r2*energy;
	);
    {out}
    )
---------------------------------------------------------------------------------------------
-- Method that computes the degree of minimal problem in navigation:
-- Given: measurements = {dirA1,dirA2,dirB1,dirB2,dist1,dist2}
--        CPolynomialCoeffMatrixList = {CPolynomialCoeffMatrixA,CPolynomialCoeffMatrixB
--        R = QQ[c1,c2,s1,s2,t1,t2]
-- Returns: List of polynomials defining the minimal problem
---------------------------------------------------------------------------------------------
findDegree = method(Options => {DerivativesOfDistances => false, OrbitScenario => "Lyapunov"})
findDegree (List,List,Ring,List) := o -> (orbitConstraintsTuples,distanceConstraintTuples,R,maxDegreeList) -> (
    constraints := constructMinimalProblems(orbitConstraintsTuples,distanceConstraintTuples,R,maxDegreeList,DerivativesOfDistances => o#DerivativesOfDistances,OrbitScenario => o#OrbitScenario);
    I := ideal constraints;
    << "--Dimension of I: " << dim I << endl;
    assert ((dim I) == 0);
    degree I
    )
