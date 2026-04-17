-- Problem: 3S all C (Lyapunov quartic)
-- three spacecraft A, B, C on DISTINCT Lyapunov
-- orbits with known Jacobi constants CA, CB, CC.
-- Complete measurement graph: A-B, B-C, A-C.
-- Unknowns:  (xA, yA, xB, yB, xC, yC)
-- Equations: 3 orbit + 3 pairwise ranges
-- Expected degree: 256  (WARNING: slow)
needs "../degree-computation.m2"

xB = symbol xB; yB = symbol yB;
R = kk[xA, yA, xB, yB, xC, yC]
dABsq = rc(); dBCsq = rc(); dACsq = rc()
gA = buildGfixed(lyapQMons(xA, yA))
gB = buildGfixed(lyapQMons(xB, yB))
gC = buildGfixed(lyapQMons(xC, yC))
eqs = {
    gA, gB, gC,
    (xA-xB)^2 + (yA-yB)^2 - dABsq,
    (xB-xC)^2 + (yB-yC)^2 - dBCsq,
    (xA-xC)^2 + (yA-yC)^2 - dACsq}
I = ideal eqs
assert(dim I == 0)
assert(degree I == expected)
