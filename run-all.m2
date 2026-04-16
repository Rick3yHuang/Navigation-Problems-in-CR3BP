-- ================================================================
-- run-all.m2
--
-- Loads and runs every degree computation in the subdirectories.
-- Usage (from the DEGREE-COMPUTATION directory):
--
--   M2 --script run-all.m2
--
-- Each problem file contains `needs "../degree-computation.m2"`.
-- To make that path resolve correctly when files are loaded from
-- here, each problem subdirectory is prepended to M2's search
-- path so that "../degree-computation.m2" normalises to
-- DEGREE-COMPUTATION/degree-computation.m2.
-- ================================================================

dir = currentFileDirectory   -- absolute path of DEGREE-COMPUTATION/

-- Add every problem subdirectory to the search path
path = {
    dir | "2M1S-Lyapunov/",
    dir | "M2S-same-Lyapunov/",
    dir | "2S-knownC-LOS-Lyapunov/",
    dir | "2S-oneC-LOS-Halo/",
    dir | "3S-knownCA-LOS-Lyapunov/",
    dir | "M2S-LOS-twice-Lyapunov/",
    dir | "3S-allC-Halo/",
    dir | "2M1S-Halo/",
    dir | "3S-allrange-Lyapunov/",
    dir | "M2S-same-Halo/"} | path

-- Load shared helpers once (subsequent `needs` calls in problem
-- files resolve to the same absolute path and are no-ops)
needs(dir | "degree-computation.m2")

-- ----------------------------------------------------------------
-- Problem files ordered by increasing expected computation time
-- ----------------------------------------------------------------
problemFiles = {
    -- fast (seconds)
    ("2M1S-Lyapunov/quartic.m2",          6),
    ("2M1S-Lyapunov/sextic.m2",           6),
    ("2S-knownC-LOS-Lyapunov/quartic.m2", 16),
    ("2S-knownC-LOS-Lyapunov/sextic.m2",  36),
    ("2S-oneC-LOS-Halo/quartic.m2",       96),
    ("2S-oneC-LOS-Halo/sextic.m2",        216),
    ("2M1S-Halo/quartic.m2",              48),
    ("2M1S-Halo/sextic.m2",               72),
    ("M2S-same-Lyapunov/quartic.m2",      84),
    ("M2S-same-Lyapunov/sextic.m2",       132),
    ("3S-knownCA-LOS-Lyapunov/quartic.m2",84),
    ("3S-knownCA-LOS-Lyapunov/sextic.m2", 198),
    ("3S-allrange-Lyapunov/quartic.m2",   256),
    ("3S-allrange-Lyapunov/sextic.m2",    864),
    -- moderate (minutes)
    ("M2S-LOS-twice-Lyapunov/quartic.m2", 1152),
    -- slow (many minutes)
    ("M2S-LOS-twice-Lyapunov/sextic.m2",  2592),
    ("3S-allC-Halo/quartic.m2",           3024),
    ("M2S-same-Halo/quartic.m2",	  4416),
    ("3S-allrange-Halo/quartic.m2",       8192)}

for f'expected in problemFiles do (
    (f, expected) := f'expected;
    print("--- " | f | "  (expected degree " | toString expected | ") ---");
    time load(dir | f);
    print "  OK";
    )

print "=== All assertions passed ==="
