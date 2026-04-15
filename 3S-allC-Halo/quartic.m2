-- Problem: 3S all unknown C (Halo quartic)
-- Three Halo spacecraft A, B, D on DISTINCT orbits with unknown
-- Jacobi constants CA, CB, CD.  LOS from A to B and from A to D
-- are known (3D directions), along with ranges dAB, dAD.
-- Positions: (uB,vB,wB) = (uA+duB, vA+dvB, wA+dwB)  and similarly D.
-- wA is determined by h(uA,vA,CA); substituting gives h-constraints.
-- Unknowns:  (uA, vA, CA, CB, CD)
-- Equations: g(uA,vA,CA)=1,
--            g(uA+duB, vA+dvB, CB)=1,
--            h(uA+duB, vA+dvB, CB) - h(uA,vA,CA) = dwB,
--            g(uA+duD, vA+dvD, CD)=1,
--            h(uA+duD, vA+dvD, CD) - h(uA,vA,CA) = dwD
-- Expected degree: 3024  (WARNING: very slow)
needs "../degree-computation.m2"

R = kk[uA, vA, CA, CB, CD]
duB = rc(); dvB = rc(); dwB = rc()
duD = rc(); dvD = rc(); dwD = rc()
uB = uA + duB; vB = vA + dvB
uD = uA + duD; vD = vA + dvD
-- A's orbit model (CA unknown)
cGA = makeCubics(#haloGQMons(uA, vA), CA)
cHA = makeCubics(#haloHQMons(uA, vA), CA)
gA = buildGcoeffs(haloGQMons(uA, vA), cGA)
hA = buildHcoeffs(haloHQMons(uA, vA), cHA)
-- B's orbit model (CB unknown)
cGB = makeCubics(#haloGQMons(uB, vB), CB)
cHB = makeCubics(#haloHQMons(uB, vB), CB)
gB = buildGcoeffs(haloGQMons(uB, vB), cGB)
hB = buildHcoeffs(haloHQMons(uB, vB), cHB)
-- D's orbit model (CD unknown)
cGD = makeCubics(#haloGQMons(uD, vD), CD)
cHD = makeCubics(#haloHQMons(uD, vD), CD)
gD = buildGcoeffs(haloGQMons(uD, vD), cGD)
hD = buildHcoeffs(haloHQMons(uD, vD), cHD)
eqs = {gA, gB, hB - hA - dwB, gD, hD - hA - dwD}
assert(degree ideal eqs == 3024)
