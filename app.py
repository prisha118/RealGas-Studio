"""
app.py -- interactive web front end for the cubic equation-of-state program.

    streamlit run app.py

Everything shown here is computed by the ``pvt_eos`` package; this file only
collects the user's input and draws the results.

    CHL2002 Chemical Engineering Thermodynamics
    Prisha Mahawar  --  2025CH12177
"""

from __future__ import annotations

import math

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import streamlit as st

from pvt_eos import (R, compute, list_components, list_eos, pressure,
                     saturation_pressure, saturation_properties, solve_state,
                     to_kelvin, to_pascal)
from pvt_eos.database import get_component
from pvt_eos.eos import EOS_MODELS

COLORS = {"vdw": "#B3446C", "rk": "#E8A33D", "srk": "#2E8B57", "pr": "#1F4E99"}

st.set_page_config(page_title="Cubic EoS - PVT, H, S and fugacity",
                   page_icon="*", layout="wide")

plt.rcParams.update({"font.size": 8, "axes.grid": True, "grid.alpha": 0.3,
                     "legend.frameon": False, "figure.dpi": 130})


# ---------------------------------------------------------------------------
# Sidebar: all user input
# ---------------------------------------------------------------------------
st.sidebar.title("Input")

species = st.sidebar.selectbox("Species", sorted(list_components()),
                               index=sorted(list_components()).index("Methane"))
comp = get_component(species)

eos_key = st.sidebar.selectbox(
    "Equation of state", list_eos(),
    format_func=lambda k: EOS_MODELS[k].name, index=list_eos().index("pr"))

st.sidebar.markdown("---")
t_unit = st.sidebar.radio("Temperature unit", ["K", "C", "F"], horizontal=True)
T_in = st.sidebar.number_input(f"Temperature / {t_unit}", value=300.0, step=5.0,
                               format="%.4f")
p_unit = st.sidebar.radio("Pressure unit", ["bar", "atm", "kPa", "MPa", "psi"],
                          horizontal=True)
P_in = st.sidebar.number_input(f"Pressure / {p_unit}", value=50.0, step=1.0,
                               min_value=1e-9, format="%.4f")

phase = st.sidebar.selectbox(
    "Root to report when three roots exist", ["auto", "vapour", "liquid"],
    help="'auto' reports the vapour (largest) root. The middle root is "
         "mechanically unstable and is always discarded.")

T = to_kelvin(T_in, t_unit)
P = to_pascal(P_in, p_unit)

st.sidebar.markdown("---")
st.sidebar.caption(
    f"**{comp.name}** ({comp.formula})  \n"
    f"Mw = {comp.Mw:.3f} g/mol  \n"
    f"Tc = {comp.Tc:.2f} K, Pc = {comp.Pc:.2f} bar  \n"
    f"omega = {comp.omega:.3f}, Zc(exp) = {comp.Zc:.3f}")


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
st.title("PVT behaviour and thermodynamic properties of real gases")
st.caption("Cubic equations of state: van der Waals | Redlich-Kwong | "
           "Soave-Redlich-Kwong | Peng-Robinson  \n"
           "CHL2002 Chemical Engineering Thermodynamics - "
           "Prisha Mahawar (2025CH12177)")

if T <= 0:
    st.error("Temperature must be above absolute zero.")
    st.stop()

try:
    res = compute(comp, eos_key, T, P, phase=phase)
except Exception as exc:
    st.error(f"The equation of state could not be solved: {exc}")
    st.stop()

st.markdown(
    f"**{comp.name}** ({comp.formula}) &nbsp;|&nbsp; {EOS_MODELS[eos_key].name}"
    f" &nbsp;|&nbsp; T = {T:.2f} K (Tr = {res.state.Tr:.4f})"
    f" &nbsp;|&nbsp; P = {P/1e5:.4f} bar (Pr = {res.state.Pr:.4f})"
    f" &nbsp;|&nbsp; root reported: **{res.phase}** "
    f"({res.state.n_roots} real root{'s' if res.state.n_roots > 1 else ''})")

# ---------------------------------------------------------------------------
# Headline numbers
# ---------------------------------------------------------------------------
c1, c2, c3, c4, c5 = st.columns(5)
c1.metric("Z", f"{res.Z:.5f}", f"{100*(res.Z-1):+.2f} % vs ideal")
c2.metric("V / cm3 mol-1", f"{res.V*1e6:.3f}")
c3.metric("H / kJ mol-1", f"{res.H/1000:.4f}", f"H_res {res.H_res/1000:+.4f}")
c4.metric("S / J mol-1 K-1", f"{res.S:.4f}", f"S_res {res.S_res:+.4f}")
c5.metric("phi = f/P", f"{res.phi:.5f}", f"f = {res.f/1e5:.4f} bar")

tab1, tab2, tab3, tab4 = st.tabs(
    ["Full results", "All four equations", "Graphs", "Saturation"])

# ---------------------------------------------------------------------------
with tab1:
    left, right = st.columns(2)
    with left:
        st.subheader("Volumetric (PVT) behaviour")
        st.dataframe(pd.DataFrame({
            "quantity": ["Z", "V / cm3 mol-1", "density / kg m-3",
                         "V ideal / cm3 mol-1", "a / Pa m6 mol-2",
                         "b / cm3 mol-1", "alpha(Tr)", "A", "B"],
            "value": [f"{res.Z:.6f}", f"{res.V*1e6:.4f}", f"{res.rho_mass:.4f}",
                      f"{R*T/P*1e6:.4f}", f"{res.state.a:.6e}",
                      f"{res.state.b*1e6:.4f}", f"{res.state.alpha:.6f}",
                      f"{res.state.A:.6f}", f"{res.state.B:.6f}"],
        }), hide_index=True, use_container_width=True)

        st.subheader("Fugacity")
        st.dataframe(pd.DataFrame({
            "quantity": ["ln(phi)", "phi = f/P", f"f / {p_unit}"],
            "value": [f"{math.log(res.phi):.6f}", f"{res.phi:.6f}",
                      f"{res.f/to_pascal(1.0, p_unit):.6f}"],
        }), hide_index=True, use_container_width=True)

    with right:
        st.subheader("Enthalpy and entropy")
        st.dataframe(pd.DataFrame({
            "quantity": ["Cp ideal / J mol-1 K-1",
                         "H ideal / J mol-1", "H residual / J mol-1",
                         "H / J mol-1",
                         "S ideal / J mol-1 K-1", "S residual / J mol-1 K-1",
                         "S / J mol-1 K-1",
                         "G = H - TS / J mol-1", "U = H - PV / J mol-1"],
            "value": [f"{res.Cp_ig:.4f}", f"{res.H_ig:.3f}", f"{res.H_res:.3f}",
                      f"{res.H:.3f}", f"{res.S_ig:.4f}", f"{res.S_res:.4f}",
                      f"{res.S:.4f}", f"{res.G:.3f}", f"{res.U:.3f}"],
        }), hide_index=True, use_container_width=True)
        st.caption("Datum: ideal gas at 298.15 K and 1 bar, with H = 0 and S = 0. "
                   "Only differences in H and S are physically meaningful.")

# ---------------------------------------------------------------------------
with tab2:
    rows = []
    for key in list_eos():
        try:
            r = compute(comp, key, T, P, phase=phase)
            rows.append({"equation of state": EOS_MODELS[key].name,
                         "Z": r.Z, "V / cm3 mol-1": r.V*1e6,
                         "H_res / J mol-1": r.H_res,
                         "S_res / J mol-1 K-1": r.S_res,
                         "H / J mol-1": r.H, "S / J mol-1 K-1": r.S,
                         "phi": r.phi})
        except Exception as exc:
            st.warning(f"{EOS_MODELS[key].name}: {exc}")
    rows.append({"equation of state": "ideal gas", "Z": 1.0,
                 "V / cm3 mol-1": R*T/P*1e6, "H_res / J mol-1": 0.0,
                 "S_res / J mol-1 K-1": 0.0,
                 "H / J mol-1": res.H_ig, "S / J mol-1 K-1": res.S_ig, "phi": 1.0})
    st.dataframe(pd.DataFrame(rows).style.format({
        "Z": "{:.6f}", "V / cm3 mol-1": "{:.4f}", "H_res / J mol-1": "{:.2f}",
        "S_res / J mol-1 K-1": "{:.4f}", "H / J mol-1": "{:.2f}",
        "S / J mol-1 K-1": "{:.4f}", "phi": "{:.6f}"}),
        hide_index=True, use_container_width=True)
    st.caption("The ideal-gas row uses Z = 1, for which both residual "
               "properties vanish and phi = 1.")

# ---------------------------------------------------------------------------
with tab3:
    # The default must stay inside the slider's own range, otherwise
    # Streamlit raises for state points above 300 bar.
    pmax_default = min(600, max(10, int(round(2 * P / 1e5 / 10.0)) * 10))
    pmax = st.slider("Maximum pressure on the graphs / bar", 10, 600,
                     pmax_default, step=10)
    grid = np.linspace(0.5e5, pmax * 1e5, 300)

    g1, g2 = st.columns(2)

    with g1:
        fig, ax = plt.subplots(figsize=(4.4, 3.1))
        for key in list_eos():
            ax.plot(grid/1e5, [solve_state(comp, key, T, p).Z for p in grid],
                    lw=1.3, color=COLORS[key], label=EOS_MODELS[key].name)
        ax.axhline(1.0, color="k", ls=":", lw=0.9, label="ideal gas")
        ax.plot(P/1e5, res.Z, "o", color="crimson", ms=5, zorder=5)
        ax.set_xlabel("P / bar"); ax.set_ylabel("Z")
        ax.set_title(f"{comp.formula} at {T:.1f} K", fontsize=9)
        ax.legend(fontsize=6.5)
        st.pyplot(fig, use_container_width=True); plt.close(fig)

        fig, ax = plt.subplots(figsize=(4.4, 3.1))
        for key in list_eos():
            ax.plot(grid/1e5, [compute(comp, key, T, p).phi for p in grid],
                    lw=1.3, color=COLORS[key], label=EOS_MODELS[key].name)
        ax.axhline(1.0, color="k", ls=":", lw=0.9)
        ax.plot(P/1e5, res.phi, "o", color="crimson", ms=5, zorder=5)
        ax.set_xlabel("P / bar"); ax.set_ylabel("phi = f/P")
        ax.set_title("Fugacity coefficient", fontsize=9)
        ax.legend(fontsize=6.5)
        st.pyplot(fig, use_container_width=True); plt.close(fig)

    with g2:
        fig, ax = plt.subplots(figsize=(4.4, 3.1))
        for key in list_eos():
            from pvt_eos import residual_enthalpy
            ax.plot(grid/1e5, [residual_enthalpy(solve_state(comp, key, T, p))/1000
                               for p in grid],
                    lw=1.3, color=COLORS[key], label=EOS_MODELS[key].name)
        ax.axhline(0.0, color="k", ls=":", lw=0.9)
        ax.plot(P/1e5, res.H_res/1000, "o", color="crimson", ms=5, zorder=5)
        ax.set_xlabel("P / bar"); ax.set_ylabel("H residual / kJ mol-1")
        ax.set_title("Enthalpy departure", fontsize=9)
        ax.legend(fontsize=6.5)
        st.pyplot(fig, use_container_width=True); plt.close(fig)

        fig, ax = plt.subplots(figsize=(4.4, 3.1))
        for key in list_eos():
            from pvt_eos import residual_entropy
            ax.plot(grid/1e5, [residual_entropy(solve_state(comp, key, T, p))
                               for p in grid],
                    lw=1.3, color=COLORS[key], label=EOS_MODELS[key].name)
        ax.axhline(0.0, color="k", ls=":", lw=0.9)
        ax.plot(P/1e5, res.S_res, "o", color="crimson", ms=5, zorder=5)
        ax.set_xlabel("P / bar"); ax.set_ylabel("S residual / J mol-1 K-1")
        ax.set_title("Entropy departure", fontsize=9)
        ax.legend(fontsize=6.5)
        st.pyplot(fig, use_container_width=True); plt.close(fig)

    st.markdown("**P-V isotherms**")
    fig, ax = plt.subplots(figsize=(6.4, 3.4))
    for Ti in sorted({0.85*comp.Tc, 0.95*comp.Tc, comp.Tc, 1.1*comp.Tc, T}):
        V = np.logspace(np.log10(1.25*res.state.b), np.log10(3e-3), 500)
        ax.plot(V*1e6, [pressure(comp, eos_key, Ti, v)/1e5 for v in V],
                lw=1.2, label=f"{Ti:.1f} K")
    ax.plot(res.V*1e6, P/1e5, "o", color="crimson", ms=5, zorder=5)
    ax.set_xscale("log"); ax.set_ylim(0, max(1.6*comp.Pc, 1.4*P/1e5))
    ax.set_xlabel("V / cm3 mol-1"); ax.set_ylabel("P / bar")
    ax.set_title(f"{EOS_MODELS[eos_key].name} isotherms for {comp.formula}",
                 fontsize=9)
    ax.legend(fontsize=6.5, ncol=5)
    st.pyplot(fig, use_container_width=True); plt.close(fig)

# ---------------------------------------------------------------------------
with tab4:
    if T >= comp.Tc:
        st.info(f"T = {T:.2f} K is at or above the critical temperature of "
                f"{comp.name} (Tc = {comp.Tc:.2f} K), so no saturation state "
                f"exists. Lower the temperature to see the vapour pressure.")
    else:
        try:
            d = saturation_properties(comp, eos_key, T)
            s1, s2, s3 = st.columns(3)
            s1.metric("Psat / bar", f"{d['Psat_bar']:.5f}")
            s2.metric("dH vap / kJ mol-1", f"{d['dH_vap_J_per_mol']/1000:.4f}")
            s3.metric("dS vap / J mol-1 K-1", f"{d['dS_vap_J_per_mol_K']:.4f}")
            st.dataframe(pd.DataFrame({
                "quantity": ["Z liquid", "Z vapour", "V liquid / cm3 mol-1",
                             "V vapour / cm3 mol-1", "phi (both phases)"],
                "value": [f"{d['Z_liquid']:.6f}", f"{d['Z_vapour']:.6f}",
                          f"{d['V_liquid_cm3_per_mol']:.4f}",
                          f"{d['V_vapour_cm3_per_mol']:.4f}",
                          f"{d['phi_sat']:.6f}"]}),
                hide_index=True, use_container_width=True)
            st.caption("Psat is found by solving phi_liquid = phi_vapour, which "
                       "is the equal-fugacity (Maxwell) criterion.")

            Ts = np.linspace(max(0.45*comp.Tc, 0.5*T), comp.Tc - 0.5, 60)
            fig, ax = plt.subplots(figsize=(5.6, 3.2))
            for key in list_eos():
                try:
                    ax.semilogy(Ts, [saturation_pressure(comp, key, t)/1e5
                                     for t in Ts],
                                lw=1.3, color=COLORS[key],
                                label=EOS_MODELS[key].name)
                except Exception:
                    pass
            ax.plot(T, d["Psat_bar"], "o", color="crimson", ms=5, zorder=5)
            ax.plot(comp.Tc, comp.Pc, "k*", ms=9, label="critical point")
            ax.set_xlabel("T / K"); ax.set_ylabel("Psat / bar")
            ax.set_title(f"Vapour-pressure curve of {comp.formula}", fontsize=9)
            ax.legend(fontsize=6.5)
            st.pyplot(fig, use_container_width=True); plt.close(fig)
        except Exception as exc:
            st.warning(f"The saturation calculation did not converge: {exc}")

st.markdown("---")
st.caption("Critical constants and acentric factors from Table B.1, "
           "ideal-gas heat capacities from Table C.1 of Smith, Van Ness, "
           "Abbott and Swihart, *Introduction to Chemical Engineering "
           "Thermodynamics*.")
