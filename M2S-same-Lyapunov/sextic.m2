-- Problem: M2S same orbit (Lyapunov sextic)
-- One mothership; two spacecraft A, B on the SAME unknown orbit
-- (shared orbit model with identical coefficients alpha_j(C)).
-- Unknowns:  (C, xA, yA, xB, yB)
-- Equations: g(xA,yA,C)=1,  g(xB,yB,C)=1,
--            ||pA-M||^2=dAM^2,  ||pB-M||^2=dBM^2,
--            ||pA-pB||^2=dAB^2
-- Expected degree: 132
needs "../degree-computation.m2"

R = kk[C, xA, yA, xB, yB]
xM = rc(); yM = rc()
dAMsq = rc(); dBMsq = rc(); dABsq = rc()
c = makeCubics(#lyapSMons(xA, yA), C)
gS = (xx, yy) -> buildGcoeffs(lyapSMons(xx, yy), c)
eqs = {
    gS(xA, yA), gS(xB, yB),
    (xA - xM)^2 + (yA - yM)^2 - dAMsq,
    (xB - xM)^2 + (yB - yM)^2 - dBMsq,
    (xA - xB)^2 + (yA - yB)^2 - dABsq}
I = ideal eqs
assert(dim I == 0)
assert(degree I == expected)
