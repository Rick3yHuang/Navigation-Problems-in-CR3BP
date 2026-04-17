-- Problem: 2S known C, range + LOS (Lyapunov sextic)
-- Two spacecraft A, B with KNOWN Jacobi constants.
-- Inter-spacecraft range dAB and LOS angle theta_AB are known,
-- so B's position is determined by A's:
--   (xB, yB) = (xA - dAB*cos(theta), yA - dAB*sin(theta))
-- Unknowns:  (xA, yA)
-- Equations: g(xA,yA)=1 (CA known),  g(xB,yB)=1 (CB known)
-- Expected degree: 36
needs "../degree-computation.m2"

R = kk[xA, yA]
cosT = rc(); sinT = rc(); dAB = rc()
xB = xA - dAB * cosT
yB = yA - dAB * sinT
eqs = {
    buildGfixed(lyapSMons(xA, yA)),
    buildGfixed(lyapSMons(xB, yB))}
I = ideal eqs
assert(dim I == 0)
assert(degree I == expected)

