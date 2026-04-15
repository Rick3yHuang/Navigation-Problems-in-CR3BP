-- Problem: M+3S all ranges (Lyapunov sextic)
-- One mothership M; three spacecraft A, B, C on DISTINCT Lyapunov
-- orbits with unknown Jacobi constants CA, CB, CC.
-- Complete measurement graph: M-A, M-B, M-C, A-B, B-C, A-C.
-- Unknowns:  (CA, xA, yA, CB, xB, yB, CC, xC, yC)
-- Equations: 3 orbit + 3 mothership-ranges + 3 pairwise ranges
-- Expected degree: 864  (WARNING: very slow)
needs "../degree-computation.m2"

R = kk[CA, xA, yA, CB, xB, yB, CC, xC, yC]
xM = rc(); yM = rc()
dMAsq = rc(); dMBsq = rc(); dMCsq = rc()
dABsq = rc(); dBCsq = rc(); dACsq = rc()
cA = makeCubics(#lyapSMons(xA, yA), CA)
cB = makeCubics(#lyapSMons(xB, yB), CB)
cC = makeCubics(#lyapSMons(xC, yC), CC)
gA = buildGcoeffs(lyapSMons(xA, yA), cA)
gB = buildGcoeffs(lyapSMons(xB, yB), cB)
gC = buildGcoeffs(lyapSMons(xC, yC), cC)
eqs = {
    gA, gB, gC,
    (xA-xM)^2 + (yA-yM)^2 - dMAsq,
    (xB-xM)^2 + (yB-yM)^2 - dMBsq,
    (xC-xM)^2 + (yC-yM)^2 - dMCsq,
    (xA-xB)^2 + (yA-yB)^2 - dABsq,
    (xB-xC)^2 + (yB-yC)^2 - dBCsq,
    (xA-xC)^2 + (yA-yC)^2 - dACsq}
assert(degree ideal eqs == 864)
