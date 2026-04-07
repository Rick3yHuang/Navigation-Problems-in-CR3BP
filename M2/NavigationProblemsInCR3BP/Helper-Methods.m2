conditionNum = method()
conditionNum Matrix := M -> (
    if (ring M === RRi_53) then (
	-- M = matrix applyTable(entries M, i->midpoint i)
	-*
	v := random(RR^(numColumns M), RR^1);    
	norm_2 (M * v)
	w := random(RRi^(numRows M), RRi^1);    
	solve(transpose M * M, transpose M * w) 
	*-
	"max diameter" => max flatten applyTable(entries M, diameter) 
	) else (
	singVal := first SVD sub(M,RR);
	"condition #" => (max singVal)/(min singVal)
	)
    )

listOfVectorsToMatrix = method()
listOfVectorsToMatrix List := l -> (
    nestedList := apply(#l, i -> flatten entries l_i);
    transpose matrix nestedList
    )

normalize = method()
normalize Matrix := M -> (
    assert(numColumns M == 1 or numRows M == 1);
    (1/(norm_2 M)) * M
    )

modelPolynomial = method(Options => {FittingHeight => false})
modelPolynomial (Ring,Matrix,ZZ,ZZ) := o -> (R,CPolynomialCoeffMatrix,jacobiConstantDegree,modelDegree) -> (
    (x,y,c) := toSequence gens R;
    jacobiConstantBasis := transpose matrix{apply(jacobiConstantDegree+1, i -> c^i)};
    coeffs := (sub(CPolynomialCoeffMatrix,R))*jacobiConstantBasis;
    monoms := matrix{returnModelBasis({x,y},modelDegree)};
    if not o#FittingHeight then out := sum flatten entries(monoms*coeffs) - 1;
    if o#FittingHeight then (
	monoms = matrix{{1}} | monoms;
	out = sum flatten entries(monoms*coeffs);
	);
    out
    )

baseModelPolynomial = method(Options => {FittingHeight => false})
baseModelPolynomial (Ring,Matrix,ZZ,ZZ) := o -> (R,AllCoeffMatrix,orbitID,modelDegree) -> (
    (x,y) := toSequence gens R;
    coeffs := transpose AllCoeffMatrix^{orbitID-1};
    monoms := matrix{returnModelBasis({x,y},modelDegree)};
    out := sum flatten entries(monoms*coeffs) - 1;
    if o#FittingHeight then (
	monoms = matrix{{1}} | monoms;
	out = sum flatten entries(monoms*coeffs);
	);
    out
    )

findEffectivePotential = method(Options => {OrbitScenario => "Lyapunov"})
findEffectivePotential (Matrix,RR) := o ->  (pos,mu) -> (
    (r1,r2) := toSequence findEarthMoonDistance(pos,mu,OrbitScenario => o#OrbitScenario);
    if o#OrbitScenario == "Lyapunov" then (
	(x,y) := toSequence flatten entries pos;
	out := 0.5*(x^2 + y^2) + (1 - mu)/r1 + mu/r2;
	) else if o#OrbitScenario == "Halo" then (
	z := 0;
	(x,y,z) = toSequence flatten entries pos;
	out = 0.5*(x^2 + y^2 + z^2) + (1 - mu)/r1 + mu/r2;
	);
    out
    )

findEffectivePotential (Matrix,Matrix,RR) := o -> (pos,sqrtDist,mu) -> (
    if o#OrbitScenario == "Lyapunov" then (
	(x,y) := toSequence flatten entries pos;
	(r1,r2) := toSequence sqrtDist;
	out := 0.5*(x^2 + y^2) + (1 - mu)/r1 + mu/r2
	) else if o#OrbitScenario == "Halo" then (
	z := 0;
	(x,y,z) = toSequence flatten entries pos;
	out = 0.5*(x^2 + y^2 + z^2) + (1 - mu)/r1 + mu/r2;
	);
    out
    )

findAcceleration = method()
findAcceleration (Matrix,Matrix,RR) := (pos,vel,mu) -> (
    (x,y) := toSequence flatten entries pos;
    (vx,vy) := toSequence flatten entries vel;
    r1 := sqrt((x + mu)^2 + y^2);
    r2 := sqrt((x - (1 - mu))^2 + y^2);
    ax := 2*vy + x - (1 - mu)*(x + mu)/(r1^3) - mu*(x - (1 - mu))/(r2^3);
    ay := -2*vx + y - (1 - mu)*y/(r1^3) - mu*y/(r2^3);
    matrix{{ax},{ay}}
    )

findEarthMoonDistance = method(Options => {OrbitScenario => "Lyapunov"})
findEarthMoonDistance (Matrix,RR) := o -> (pos,mu) -> (
    if o#OrbitScenario == "Lyapunov" then (
	(x,y) := toSequence flatten entries pos;
	r1 := sqrt((x + mu)^2 + y^2);
	r2 := sqrt((x - (1 - mu))^2 + y^2);
	) else if o#OrbitScenario == "Halo" then (
	z := 0;
	(x,y,z) = toSequence flatten entries pos;
	r1 = sqrt((x + mu)^2 + y^2 + z^2);
	r2 = sqrt((x - (1 - mu))^2 + y^2 + z^2);
	);
    {r1,r2}
    )

returnModelBasis = method()
returnModelBasis (List,ZZ) := (variables,maxDegree) -> (
    assert(#variables == 2);    
    originalRing := ring variables_0;
    R := ZZ[variables_1,variables_0,MonomialOrder => Lex];
    (y,x) := toSequence flatten entries vars R;
    modelBasis := {};
    scan(maxDegree+1, i-> (
	    scan(floor ((maxDegree-i)/2+1), j -> (
		    if i+j != 0 then modelBasis = append(modelBasis, x^(i)*y^(2*j));
		    )))
	);
    use originalRing;
    apply(sort(modelBasis),a -> sub(a,originalRing))
    )

returnAllModelBasis = method()
returnAllModelBasis (List,ZZ,ZZ) := (variables,modelDegree,jacobiConstantDegree) -> (
    JCBasis := transpose matrix{apply(jacobiConstantDegree+1, i -> (variables_2)^i)};
    modelBasis := returnModelBasis(variables_{0..1},modelDegree);
    flatten entries (JCBasis * (matrix {modelBasis}))
    )

findModelExponentMatrix = method()
findModelExponentMatrix ZZ := modelDegree -> (
    modelBasis := {};
    scan(floor(modelDegree/2)+1, j-> (
	    scan(modelDegree-2*j+1, i -> (
		    if i+j != 0 then modelBasis = append(modelBasis, {i,2*j});
		    )
		)
	    )
	);
    matrix modelBasis
    )
