-- Problem: M2S same orbit (Halo quartic)
-- One mothership; two spacecraft A, B on the SAME unknown orbit
-- (shared orbit model with identical coefficients alpha_j(C)).
-- Unknowns:  (C, uA, vA, uB, vB)
-- Equations: g(uA,vA,C)=1,  g(uB,vB,C)=1,
--            ||pA-M||^2=dAM^2,  ||pB-M||^2=dBM^2,
--            ||pA-pB||^2=dAB^2
-- Expected degree: 4416 (WARNING: slow)
needs "../degree-computation.m2"
uB = symbol uB; vB = symbol vB;
R = kk[C, uA, vA, uB, vB]
uM = rc(); vM = rc(), hM = rc()
dAMsq = rc(); dBMsq = rc(); dABsq = rc()
cG = makeCubics(#haloGQMons(uA, vA), C)
cH = makeCubics(#haloHQMons(uA, vA), C)
gS = (uu, vv) -> buildGcoeffs(haloGQMons(uu, vv), cG)
hS = (uu, vv) -> buildHcoeffs(haloHQMons(uu, vv), cH)
hA = hS(uA, vA); hB = hS(uB, vB)
eqs = {
    gS(uA, vA), gS(uB, vB),
    (uA - uM)^2 + (uA - uM)^2 + (hA - hM)^2 - dAMsq,
    (uB - vM)^2 + (uB - vM)^2 + (hB - hM)^2 - dBMsq,
    (uA - vB)^2 + (uA - vB)^2 + (hA - hB)^2 - dABsq}
I = ideal eqs
assert(dim I == 0)
assert(degree I == expected)
