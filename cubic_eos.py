#!/usr/bin/env python3
"""
cubic_eos.py -- top-level entry point of the project.

    PVT behaviour and thermodynamic properties (H, S) and pure-component
    fugacity of real gases from cubic equations of state.

    CHL2002 Chemical Engineering Thermodynamics
    Prisha Mahawar  --  2025CH12177

The numerical engine lives in the ``pvt_eos`` package, which is organised as

    pvt_eos/constants.py    physical constants and unit conversion
    pvt_eos/database.py     critical constants and ideal-gas Cp (Tables B.1, C.1)
    pvt_eos/eos.py          the generalised cubic EoS and its analytical solver
    pvt_eos/properties.py   residual H and S, absolute H and S, fugacity
    pvt_eos/saturation.py   vapour pressure from the equal-fugacity criterion
    pvt_eos/cli.py          the command-line interface

This file simply exposes that engine under the single, familiar name used
throughout the report, so the whole program can be driven with

    python cubic_eos.py point -s methane -e pr -T 200 -P 50
    python cubic_eos.py point -s CO2 -T 320 -P 80 --all-eos
    python cubic_eos.py table -s ethane -e srk -T 350 --p-range 1 100 20
    python cubic_eos.py sat   -s propane -e pr -T 300
    python cubic_eos.py list
    python cubic_eos.py                      (menu-driven interactive mode)
"""

from pvt_eos import (                                   # noqa: F401
    R, Component, EOSState, PropertyResult,
    compute, cp_ideal, enthalpy_ideal, entropy_ideal,
    fugacity, fugacity_coefficient, ln_fugacity_coefficient,
    get_component, get_eos, list_components, list_eos,
    pressure, residual_enthalpy, residual_entropy, residual_gibbs,
    saturation_pressure, saturation_properties,
    solve_cubic, solve_state, solve_Z,
    to_kelvin, to_pascal, from_pascal,
)
from pvt_eos.cli import main

__version__ = "1.0.0"

if __name__ == "__main__":
    raise SystemExit(main())
