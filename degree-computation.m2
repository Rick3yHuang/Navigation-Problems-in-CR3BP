-- ================================================================
-- degree-computation.m2  (shared helpers)
--
-- Orbit models and polynomial builders for minimal-problem
-- degree computations.  Loaded via  needs "../degree-computation.m2"
-- from each problem subdirectory.
--
-- Orbit models (even powers of y / v only):
--   Lyapunov quartic  g(x,y,C)-1  -- 8 monomials
--   Lyapunov sextic   g(x,y,C)-1  -- 15 monomials
--   Halo quartic      g(u,v,C)-1  -- 8 monomials (same as Lyap Q)
--                     h(u,v,C)    -- 9 monomials (constant added)
--   Halo sextic       g           -- 15 monomials
--                     h           -- 16 monomials (constant added)
--
-- Each physical coefficient alpha_j(C) is a random cubic in C
-- (matching the paper's two-stage fitting model).
-- ================================================================

-- the prime for the finite field could be changed to a larger prime if
-- needed to show the correctness of the degree computations.  The current
-- prime is small enough to make the computations fast, but large enough
-- to avoid collisions in the random coefficients. (We tried 7772777, and
-- the degree computations were still correct.)
kk = ZZ/32003
--kk = ZZ/7772777

-- Random element of kk
rc = () -> random kk

-- ----------------------------------------------------------------
-- Monomial lists
-- ----------------------------------------------------------------

-- Lyapunov quartic: 8 monomials, even y-powers, total deg <= 4
lyapQMons = (x, y) -> {x, x^2, x^3, x^4, y^2, x*y^2, x^2*y^2, y^4}

-- Lyapunov sextic: 15 monomials, even y-powers, total deg <= 6
lyapSMons = (x, y) -> {
    x, x^2, x^3, x^4, x^5, x^6,
    y^2, x*y^2, x^2*y^2, x^3*y^2, x^4*y^2,
    y^4, x*y^4, x^2*y^4,
    y^6}

-- Halo g-monomials: same structure as Lyapunov (u plays x, v plays y)
haloGQMons = (u, v) -> lyapQMons(u, v)
haloGSMons = (u, v) -> lyapSMons(u, v)

-- Halo h-monomials: constant prepended  (quartic: 9, sextic: 16)
haloHQMons = (u, v) -> prepend(1_(ring u), lyapQMons(u, v))
haloHSMons = (u, v) -> prepend(1_(ring u), lyapSMons(u, v))

-- ----------------------------------------------------------------
-- Polynomial builders
-- ----------------------------------------------------------------

-- g(x,y,C)-1 with C as a ring variable.
-- Each coefficient is an independent random cubic in C.
buildOrbit = (mons, C) -> (
    n := #mons;
    sum(0..n-1, i -> (sum(0..3, j -> (rc()) * C^j)) * mons_i) - 1)

-- h(u,v,C)  (no -1; caller subtracts w as needed)
buildH = (mons, C) -> (
    n := #mons;
    sum(0..n-1, i -> (sum(0..3, j -> (rc()) * C^j)) * mons_i))

-- g(x,y)-1 with C already substituted (coefficients are field elements)
buildGfixed = mons -> sum(0..#mons-1, i -> (rc()) * mons_i) - 1

-- h(u,v) with C already substituted
buildHfixed = mons -> sum(0..#mons-1, i -> (rc()) * mons_i)

-- n independent random cubics in C  (for shared orbit models)
makeCubics = (n, C) -> for i from 0 to n-1 list sum(0..3, j -> (rc()) * C^j)

-- g from pre-computed coefficient list
buildGcoeffs = (mons, coeffs) -> sum(0..#mons-1, i -> coeffs_i * mons_i) - 1

-- h from pre-computed coefficient list  (caller subtracts w)
buildHcoeffs = (mons, coeffs) -> sum(0..#mons-1, i -> coeffs_i * mons_i)
