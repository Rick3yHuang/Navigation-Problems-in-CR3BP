-- Problem: M2S LOS two time instances (Lyapunov quartic)
-- Mothership M at known positions at times t=1,2; spacecraft A, B
-- on DISTINCT unknown orbits CA, CB.  At each time, M measures
-- LOS to A and B (each LOS direction + unknown scalar range gives
-- the spacecraft position) and the inter-spacecraft range.
-- Unknowns:  (CA, CB, dMA1, dMA2, dMB1, dMB2)
-- Equations: g(pA^t, CA)=1 for t=1,2;  g(pB^t, CB)=1 for t=1,2;
--            ||pA^t - pB^t||^2 = dAB_t^2 for t=1,2
-- Expected degree: 1152  (WARNING: slow)
needs "../degree-computation.m2"

R = kk[CA, CB, dMA1, dMA2, dMB1, dMB2]
-- Mothership positions and LOS unit vectors at each time (known)
xM1 = rc(); yM1 = rc();  cA1 = rc(); sA1 = rc();  cB1 = rc(); sB1 = rc()
xM2 = rc(); yM2 = rc();  cA2 = rc(); sA2 = rc();  cB2 = rc(); sB2 = rc()
dAB1sq = rc(); dAB2sq = rc()
-- Spacecraft positions as polynomials in the scalar-range unknowns
xA1 = xM1 + dMA1*cA1;  yA1 = yM1 + dMA1*sA1
xA2 = xM2 + dMA2*cA2;  yA2 = yM2 + dMA2*sA2
xB1 = xM1 + dMB1*cB1;  yB1 = yM1 + dMB1*sB1
xB2 = xM2 + dMB2*cB2;  yB2 = yM2 + dMB2*sB2
-- Distinct orbit models for A and B
cA = makeCubics(#lyapQMons(xA1, yA1), CA)
cB = makeCubics(#lyapQMons(xB1, yB1), CB)
gA = (xx, yy) -> buildGcoeffs(lyapQMons(xx, yy), cA)
gB = (xx, yy) -> buildGcoeffs(lyapQMons(xx, yy), cB)
eqs = {
    gA(xA1, yA1), gA(xA2, yA2),
    gB(xB1, yB1), gB(xB2, yB2),
    (xA1-xB1)^2 + (yA1-yB1)^2 - dAB1sq,
    (xA2-xB2)^2 + (yA2-yB2)^2 - dAB2sq}
I = ideal eqs
assert(dim I == 0)
assert(degree I == expected)
