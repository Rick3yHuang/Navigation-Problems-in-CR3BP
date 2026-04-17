-- Problem: 3S all C (Halo quartic)
-- three spacecraft A, B, C on DISTINCT Lyapunov
-- orbits with known Jacobi constants CA, CB, CC.
-- Complete measurement graph: A-B, B-C, A-C.
-- Unknowns:  (uA, vA, uB, vB, uC, vC)
-- Equations: 3 orbit + 3 pairwise ranges
-- Expected degree: 8192  (WARNING: slow)
needs "../degree-computation.m2"

uB = symbol uB; vB = symbol vB;
R = kk[uA, vA, uB, vB, uC, vC]
dABsq = rc(); dBCsq = rc(); dACsq = rc()
gA = buildGfixed(haloGQMons(uA, vA))
hA = buildHfixed(haloHQMons(uA, vA))
gB = buildGfixed(haloGQMons(uB, vB))
hB = buildHfixed(haloHQMons(uB, vB))
gC = buildGfixed(haloGQMons(uC, vC))
hC = buildHfixed(haloHQMons(uC, vC))
eqs = {
    gA, gB, gC,
    (uA-vB)^2 + (uA-vB)^2 + (hA-hB)^2 - dABsq,
    (uB-vC)^2 + (uB-vC)^2 + (hB-hC)^2 - dBCsq,
    (uA-vC)^2 + (uA-vC)^2 + (hA-hC)^2 - dACsq}
I = ideal eqs
assert(dim I == 0)
assert(degree I == expected)
