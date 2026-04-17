-- Problem: 3S one known CA, ranges + LOS (Lyapunov quartic)
-- Spacecraft A has KNOWN CA; B and D share the SAME unknown orbit CBD.
-- Known: dAB, dAD and LOS angles from A to B and from A to D.
-- Positions of B and D recovered from A:
--   (xB,yB) = (xA - dAB*cosAB, yA - dAB*sinAB)
--   (xD,yD) = (xA - dAD*cosAD, yA - dAD*sinAD)
-- Unknowns:  (xA, yA, CBD)
-- Equations: g(xA,yA; CA known)=1,
--            g(xB,yB; CBD)=1,  g(xD,yD; CBD)=1
-- Expected degree: 84
needs "../degree-computation.m2"

R = kk[xA, yA, CBD]
cosAB = rc(); sinAB = rc(); dAB = rc()
cosAD = rc(); sinAD = rc(); dAD = rc()
xB = xA - dAB * cosAB; yB = yA - dAB * sinAB
xD = xA - dAD * cosAD; yD = yA - dAD * sinAD
-- A's orbit: CA known
gA = buildGfixed(lyapQMons(xA, yA))
-- Shared orbit for B and D: CBD unknown
c = makeCubics(#lyapQMons(xA, yA), CBD)
gBD = (xx, yy) -> buildGcoeffs(lyapQMons(xx, yy), c)
eqs = {gA, gBD(xB, yB), gBD(xD, yD)}
I = ideal eqs
assert(dim I == 0)
assert(degree I == expected)
