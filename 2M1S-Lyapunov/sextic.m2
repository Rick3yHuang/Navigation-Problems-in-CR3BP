-- Problem: 2M+1S (Lyapunov sextic)
-- Two motherships at known positions; one spacecraft on an
-- unknown Lyapunov orbit.
-- Unknowns:  (x, y, C)
-- Equations: orbit model, range to M1, range to M2
-- Expected degree: 6
needs "../degree-computation.m2"

R = kk[x, y, C]
xM1 = rc(); yM1 = rc(); xM2 = rc(); yM2 = rc()
d1sq = rc(); d2sq = rc()
eqs = {
    buildOrbit(lyapSMons(x, y), C),
    (x - xM1)^2 + (y - yM1)^2 - d1sq,
    (x - xM2)^2 + (y - yM2)^2 - d2sq}
I = ideal eqs
assert(dim I == 0)
assert(degree I == expected)

