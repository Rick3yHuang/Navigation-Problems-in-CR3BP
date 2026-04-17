-- Problem: 2M+1S (Halo sextic)
-- Two motherships at known 3D positions; one spacecraft on an
-- unknown Halo orbit.  The height w = h(u,v,C) is substituted
-- into the 3D distance equations.
-- Unknowns:  (u, v, C)
-- Equations: g(u,v,C)=1,
--            ||p-M1||^2 = d1^2  (with w = h(u,v,C)),
--            ||p-M2||^2 = d2^2  (with w = h(u,v,C))
-- Expected degree: 72
needs "../degree-computation.m2"

R = kk[u, v, C]
uM1 = rc(); vM1 = rc(); wM1 = rc()
uM2 = rc(); vM2 = rc(); wM2 = rc()
d1sq = rc(); d2sq = rc()
g1  = buildOrbit(haloGSMons(u, v), C)
hw  = buildH(haloHSMons(u, v), C)     -- w = h(u,v,C)
eqs = {
    g1,
    (u-uM1)^2 + (v-vM1)^2 + (hw-wM1)^2 - d1sq,
    (u-uM2)^2 + (v-vM2)^2 + (hw-wM2)^2 - d2sq}
I = ideal eqs
assert(dim I == 0)
assert(degree I == expected)
