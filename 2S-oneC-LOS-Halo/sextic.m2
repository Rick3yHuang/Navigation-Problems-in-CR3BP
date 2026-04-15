-- Problem: 2S one unknown C, range + LOS (Halo sextic)
-- Spacecraft A has KNOWN CA; spacecraft B has UNKNOWN CB.
-- Range dAB and 3D LOS direction (du, dv, dw) from A to B known.
-- B's position: (uB,vB,wB) = (uA+du, vA+dv, wA+dw).
-- wA = h(uA,vA; CA); substituting gives the h-constraint below.
-- Unknowns:  (uA, vA, CB)
-- Equations: g(uA,vA; CA known)=1,
--            g(uA+du, vA+dv; CB)=1,
--            h(uA+du, vA+dv; CB) - h(uA,vA; CA known) = dw
-- Expected degree: 108
needs "../degree-computation.m2"

R = kk[uA, vA, CB]
du = rc(); dv = rc(); dw = rc()
uB = uA + du; vB = vA + dv
-- A's model: CA known, plain field coefficients
gA = buildGfixed(haloGSMons(uA, vA))
hA = buildHfixed(haloHSMons(uA, vA))
-- B's model: CB unknown, cubic coefficients in CB
cGB = makeCubics(#haloGSMons(uB, vB), CB)
cHB = makeCubics(#haloHSMons(uB, vB), CB)
gB = buildGcoeffs(haloGSMons(uB, vB), cGB)
hB = buildHcoeffs(haloHSMons(uB, vB), cHB)
eqs = {gA, gB, hB - hA - dw}
assert(degree ideal eqs == 108)
