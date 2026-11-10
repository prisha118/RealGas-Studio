#!/usr/bin/env python3
"""
test_cubic_eos.py -- run the whole verification suite.

    python test_cubic_eos.py

The individual test modules live in ``tests/`` and are written for pytest.
This wrapper runs them with pytest when it is installed and falls back to a
self-contained set of checks when it is not, so the suite can always be
executed with nothing but the standard library plus the project itself.
"""

from __future__ import annotations

import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)


def _run_with_pytest() -> int:
    import pytest
    return pytest.main(["-q", os.path.join(HERE, "tests")])


def _run_fallback() -> int:
    """Minimal self-checking suite used when pytest is unavailable."""
    import math

    from pvt_eos import (R, compute, list_eos, pressure, residual_enthalpy,
                         residual_entropy, residual_gibbs,
                         saturation_pressure, solve_cubic, solve_state,
                         solve_Z)
    from pvt_eos.eos import EOS_MODELS
    from pvt_eos.database import get_component

    passed = failed = 0

    def check(name, condition):
        nonlocal passed, failed
        if condition:
            passed += 1
            print(f"  PASS  {name}")
        else:
            failed += 1
            print(f"  FAIL  {name}")

    print("Verification of cubic_eos\n" + "=" * 60)

    roots = solve_cubic(-6.0, 11.0, -6.0)
    check("cubic solver finds the roots 1, 2, 3",
          all(abs(a - b) < 1e-9 for a, b in zip(roots, [1.0, 2.0, 3.0])))

    for key in list_eos():
        st = solve_state("nitrogen", key, 300.0, 1.0)
        check(f"{key}: Z -> 1 as P -> 0", abs(st.Z - 1.0) < 1e-5)

    for key in list_eos():
        st = solve_state("CO2", key, 400.0, 30e5)
        check(f"{key}: P(T, V) round-trip",
              abs(pressure("CO2", key, 400.0, st.V) - 30e5) < 1e-4)

    for key, Zc in [("vdw", 0.375), ("rk", 1 / 3), ("srk", 1 / 3),
                    ("pr", 0.30740130493414073)]:
        c = get_component("methane")
        st = solve_state(c, key, c.Tc, c.Pc_pa)
        check(f"{key}: Zc = {Zc:.4f} at the critical point", abs(st.Z - Zc) < 1e-8)

    check("three real roots inside the two-phase envelope",
          len(solve_Z("propane", "pr", 300.0, 9.9e5)) == 3)
    check("the liquid root is denser than the vapour root",
          solve_state("propane", "pr", 300.0, 9.9e5, phase="liquid").V
          < solve_state("propane", "pr", 300.0, 9.9e5, phase="vapour").V)

    for key in list_eos():
        T, P, h = 320.0, 80e5, 0.032
        g = lambda t: residual_gibbs(solve_state("CO2", key, t, P)) / (R * t)
        num = -R * T * T * (g(T + h) - g(T - h)) / (2 * h)
        ana = residual_enthalpy(solve_state("CO2", key, T, P))
        check(f"{key}: H_res satisfies Gibbs-Helmholtz",
              abs(ana - num) / abs(num) < 1e-5)

        gp = residual_gibbs(solve_state("CO2", key, T + h, P))
        gm = residual_gibbs(solve_state("CO2", key, T - h, P))
        num_s = -(gp - gm) / (2 * h)
        ana_s = residual_entropy(solve_state("CO2", key, T, P))
        check(f"{key}: S_res = -d(G_res)/dT",
              abs(ana_s - num_s) / abs(num_s) < 1e-5)

    # ln(phi) against the pressure integral of (Z-1)/P
    T, P, n = 300.0, 60e5, 400
    dP, total = P / n, 0.0
    for i in range(n + 1):
        Pi = i * dP if i else 1e-6
        w = 1 if i in (0, n) else (4 if i % 2 else 2)
        total += w * (solve_state("CO2", "pr", T, Pi).Z - 1.0) / Pi
    numeric = total * dP / 3.0
    analytic = math.log(compute("CO2", "pr", T, P).phi)
    check("ln(phi) equals the integral of (Z-1)/P dP",
          abs(analytic - numeric) < 1e-3)

    for species, T, exp in [("water", 373.15, 1.0133), ("propane", 300.0, 9.9742),
                            ("ammonia", 300.0, 10.613)]:
        p = saturation_pressure(species, "pr", T) / 1e5
        check(f"Psat({species}, {T:.2f} K) = {exp:.4f} bar within 10 %",
              abs(p - exp) / exp < 0.10)

    r = compute("methane", "pr", 200.0, 50e5)
    check("H = H_ideal + H_res", abs(r.H - (r.H_ig + r.H_res)) < 1e-9)
    check("S = S_ideal + S_res", abs(r.S - (r.S_ig + r.S_res)) < 1e-12)
    check("f = phi * P", abs(r.f - r.phi * 50e5) < 1e-6)

    print("=" * 60)
    print(f"{passed} passed, {failed} failed")
    return 1 if failed else 0


if __name__ == "__main__":
    try:
        import pytest                                   # noqa: F401
    except ImportError:
        raise SystemExit(_run_fallback())
    raise SystemExit(_run_with_pytest())
