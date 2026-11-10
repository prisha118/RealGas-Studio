# CHL2002 — PVT Behaviour and Thermodynamic Properties of Real Gases from Cubic Equations of State

**Prisha Mahawar — 2025CH12177**
CHL2002: Chemical Engineering Thermodynamics

A program that predicts the **PVT behaviour**, the **enthalpy and entropy**, and the
**pure-component fugacity** of real gases from four cubic equations of state:
**van der Waals**, **Redlich–Kwong**, **Soave–Redlich–Kwong** and **Peng–Robinson**.

**Live demo (no installation):** <https://PRISHA-GITHUB-USERNAME.github.io/CHL2002_Cubic_EOS_Project/>

---

## What the program does

| Part of the assignment | What the program computes |
|---|---|
| 1. PVT behaviour | Compressibility factor `Z`, molar volume `V`, density, all real roots of the cubic, and the phase each root belongs to |
| 2. Thermodynamic properties | Residual `H^R` and `S^R`, ideal-gas `H^ig` and `S^ig`, and absolute `H` and `S` (plus `G`, `U`, `Cp^ig`) |
| 3. Fugacity | `ln φ`, the fugacity coefficient `φ = f/P`, and the fugacity `f` |
| Extra | Vapour pressure and enthalpy of vaporisation from the equal-fugacity criterion |

Covers **45 species** with complete data, and every equation of state can be
applied to every species.

## Theory in one screen

All four equations are special cases of one generic cubic:

```
            R T              a(T)
  P  =  ---------  -  ------------------------
         V - b         (V + εb)(V + σb)

  a(T) = Ψ α(Tr) R²Tc²/Pc           b = Ω R Tc / Pc
```

| EoS | ε | σ | Ω | Ψ | α(Tr) | Zc |
|---|---|---|---|---|---|---|
| van der Waals | 0 | 0 | 1/8 | 27/64 | 1 | 0.375 |
| Redlich–Kwong | 0 | 1 | 0.08664 | 0.42748 | Tr^(−1/2) | 1/3 |
| Soave–RK | 0 | 1 | 0.08664 | 0.42748 | [1+m(1−Tr^½)]² | 1/3 |
| Peng–Robinson | 1−√2 | 1+√2 | 0.07780 | 0.45724 | [1+m(1−Tr^½)]² | 0.30740 |

With `A = aP/(RT)²` and `B = bP/(RT)` this becomes a cubic in `Z`:

```
Z³ + [(ε+σ−1)B − 1]Z² + [A + εσB² − (ε+σ)B(B+1)]Z − [AB + (B+1)εσB²] = 0
```

solved in closed form (Cardano) and polished by Newton's method. The residual
Gibbs energy then gives everything else:

```
G^R/RT = Z − 1 − ln(Z−β) − q·I            q = Ψα(Tr)/(Ω·Tr)
H^R/RT = Z − 1 + [dlnα/dlnTr − 1]·q·I     I = ln[(Z+σβ)/(Z+εβ)]/(σ−ε)
S^R/R  = ln(Z−β) + [dlnα/dlnTr]·q·I       ln φ = G^R/RT
```

Absolute properties use an ideal gas at 298.15 K and 1 bar as the datum, with
`Cp^ig/R = A + BT + CT² + DT⁻²` integrated analytically.

## Files

| File | Contents |
|---|---|
| `cubic_eos.py` | Top-level entry point — the whole program is driven from here |
| `app.py` | Streamlit web app (tables, four-EoS comparison, graphs, saturation) |
| `test_cubic_eos.py` | Runs the full verification suite |
| `requirements.txt` | Optional extras; the core program needs no third-party packages |
| `report.pdf` | Four-page project report: theory, method, validation, results |
| `pvt_eos/` | The engine (see below) |
| `tests/` | 290 pytest checks |
| `docs/index.html` | Browser version of the program, used for the live demo |
| `tools/` | Scripts that generate the component database, the web page and the screenshots |
| `examples/make_figures.py` | Regenerates every figure and table in the report |
| `figures/` | Figures used in the report |

Inside `pvt_eos/`:

| Module | Contents |
|---|---|
| `constants.py` | `R`, unit conversion, reference state |
| `database.py` | 45 species from Tables B.1 and C.1 |
| `_table_data.py` | Generated data file (do not edit by hand) |
| `eos.py` | The generic cubic, Table 3.1 parameters, the analytical root solver |
| `properties.py` | Ideal-gas `Cp`, `H`, `S`; residual `H`, `S`, `G`; fugacity |
| `saturation.py` | Vapour pressure and latent heat by equal fugacity |
| `cli.py` | Command-line interface |

## How to use the code

No installation is required for the core program — Python 3.8+ and the standard
library are enough.

```bash
git clone https://github.com/PRISHA-GITHUB-USERNAME/CHL2002_Cubic_EOS_Project.git
cd CHL2002_Cubic_EOS_Project
```

**All properties at one state point**

```bash
python cubic_eos.py point -s methane -e pr -T 200 -P 50
```

```
  Species            : Methane (CH4), Mw = 16.043 g/mol
  Equation of state  : Peng-Robinson
  State              : T = 200.00 K (1.0493 Tc), P = 50.0000 bar (1.0872 Pc)

  Volumetric (PVT) behaviour
    Z          = 0.522928
    V          = 173.9148 cm3/mol
    V_ideal    = 332.5785 cm3/mol  (deviation -47.71 %)

  Residual (departure) properties
    H_res      =    -2530.516 J/mol      =  -2.5305 kJ/mol
    S_res      =      -9.2930 J/(mol K)

  Absolute properties  (datum: ideal gas at 298.15 K, 1.00 bar, H = S = 0)
    H          =    -5654.298 J/mol
    S          =     -54.4405 J/(mol K)

  Fugacity
    phi = f/P  =     0.667602
    f          =    33.380102 bar
```

**Compare all four equations of state**

```bash
python cubic_eos.py point -s CO2 -T 320 -P 80 --all-eos
```

**Sweep temperature or pressure and export a CSV**

```bash
python cubic_eos.py table -s ethane -e srk -T 350 --p-range 1 100 20 --csv ethane.csv
python cubic_eos.py table -s CO2 -e pr --t-range 250 450 9 -P 60
```

**Saturation pressure and latent heat**

```bash
python cubic_eos.py sat -s propane -e pr -T 300
```

**Menu-driven mode and the species list**

```bash
python cubic_eos.py            # interactive menu
python cubic_eos.py list       # the 45-species database
```

**Units.** Temperature accepts `K`, `C`, `F` via `--t-unit`; pressure accepts
`Pa`, `kPa`, `bar`, `atm`, `MPa`, `psi` via `--p-unit`. Defaults are K and bar.
Use `--phase liquid` to report the liquid root where three roots exist.

**As a library**

```python
from pvt_eos import compute, saturation_pressure

r = compute("CO2", "pr", T=320.0, P=80e5)
print(r.Z, r.H, r.S, r.phi)
print(saturation_pressure("propane", "pr", 300.0) / 1e5)
```

## How to use the website

Two versions are provided.

![The browser interface](figures/webapp_full.png)

**1. The browser version — nothing to install.** Open `docs/index.html`, or
visit the live demo link at the top. It is a port of the same equations to
JavaScript and agrees with the Python program to machine precision
(worst relative difference 1.3 × 10⁻¹⁵ over the test points). Choose a species,
an equation of state, and a temperature and pressure; the page shows `Z`, `V`,
`H`, `S` and `φ`, a side-by-side table of all four equations against the ideal
gas, and a live `Z`-against-`P` graph.

**2. The Streamlit app** — richer graphs and a saturation tab:

```bash
pip install -r requirements.txt
streamlit run app.py
```

It displays the full property breakdown, the four-equation comparison, graphs of
`Z`, `φ`, `H^R` and `S^R` against pressure, P–V isotherms, and the vapour-pressure
curve with the critical point marked.

## Tests

```bash
python test_cubic_eos.py
```

**303 checks.** They run with pytest if it is installed and fall back to a
standard-library-only suite if it is not. What they verify:

- **Against the textbook.** Example 3.9 — molar volumes of saturated-vapour and
  saturated-liquid n-butane at 350 K and 9.4573 bar — is reproduced *exactly*
  for all four equations, along with the intermediate `Tr`, `Pr`, `q`, `β` and `Z`
  quoted in the solution.

  | cm³/mol | Exp. | vdW | RK | SRK | PR |
  |---|---|---|---|---|---|
  | V (vapour) — book | 2482 | 2667 | 2555 | 2520 | 2486 |
  | V (vapour) — program | | 2667.0 | 2555.5 | 2520.4 | 2486.5 |
  | V (liquid) — book | 115.0 | 191.0 | 133.3 | 127.8 | 112.6 |
  | V (liquid) — program | | 191.0 | 133.3 | 127.8 | 112.6 |

- **A second textbook benchmark.** Table 13.7 (after Example 13.5) gives `Z`,
  `H^R` and `S^R` for n-butane at 500 K and 50 bar for all four equations —
  so it validates the *departure functions*, not just the volumetric part.
  All twelve values match within 0.1 %, as do the intermediate
  `β = 0.09703`, `q = 3.8689` and `I = 0.13247` quoted in the solution.

  | n-butane, 500 K, 50 bar | vdW | RK | SRK | PR |
  |---|---|---|---|---|
  | `Z` book / program | 0.6608 / 0.6610 | 0.6850 / 0.6852 | 0.7222 / 0.7224 | 0.6907 / 0.6909 |
  | `H^R` book / program, J/mol | −3937 / −3936 | −4505 / −4503 | −4824 / −4822 | −4988 / −4985 |
  | `S^R` book / program, J/mol·K | −5.424 / −5.421 | −6.546 / −6.542 | −7.413 / −7.408 | −7.426 / −7.421 |

- **Thermodynamic consistency.** `H^R` is checked against
  `−RT² d(G^R/RT)/dT` and `S^R` against `−d(G^R)/dT` by numerical
  differentiation (worst error 5 × 10⁻⁷), and `ln φ` against
  `∫₀^P (Z−1)/P dP` by Simpson's rule (worst error 2 × 10⁻⁴).
- **Limits.** `Z → 1` and `φ → 1` as `P → 0`; residual properties vanish;
  each equation returns its own exact `Zc` at the critical point.
- **Vapour pressure** against experiment: mean absolute deviation
  **1.0 % for PR** and **2.1 % for SRK** over seven benchmarks.
- **Database integrity.** `Cp^ig/R` recomputed from the `A, B, C, D`
  coefficients matches the `Cp°₂₉₈/R` column of Table C.1 for every species
  (see the note below).

## Data

Critical constants and acentric factors are from **Table B.1** and ideal-gas
heat capacities from **Table C.1** of Smith, Van Ness, Abbott & Swihart,
*Introduction to Chemical Engineering Thermodynamics*, 8th ed.

The rows are **not typed in by hand**. `tools/extract_tables.py` parses both
tables straight out of the textbook PDF and writes `pvt_eos/_table_data.py`, so
transcription errors are impossible. A species is included only when it appears
in both tables, which gives 45 components with a complete data set.

```bash
python tools/extract_tables.py path/to/textbook.pdf > pvt_eos/_table_data.py
python tools/build_webapp.py          # rebuild docs/index.html from that data
python tools/make_screenshots.py      # re-render the screenshots (headless Chrome)
python examples/make_figures.py       # regenerate every figure in the report
```

Two things the data check turned up:

- **1,3-Butadiene.** Table C.1 prints `Cp°₂₉₈/R = 10.720`, but the `A`, `B`, `C`
  coefficients on the same row give **9.931**. The coefficients are the
  self-consistent pair — they match the literature heat capacity of butadiene
  at 298 K far better — so the printed value appears to be a misprint, and the
  program uses the coefficients. This is recorded as a test so it cannot be
  forgotten.
- **Chlorine** is the least self-consistent row of Table B.1:
  `Pc·Vc/(R·Tc)` gives 0.276 against the printed `Zc = 0.265`. Every other
  species agrees to better than 0.01.

## A note on numerical precision

At the critical point the cubic has a **triple** root, so an error `δ` in the
coefficients displaces `Z` by roughly `δ^(1/3)`. Using the textbook's five-figure
`Ω` and `Ψ` therefore gives `Zc = 0.3272` for RK instead of the exact `1/3` — an
error of 6 × 10⁻³ from inputs accurate to 6 × 10⁻⁷. The program stores `Ω` and
`Ψ` to full precision (`Ω_RK = (2^⅓−1)/3`, `Ψ_RK = 1/[9(2^⅓−1)]`), which
reproduces every `Zc` exactly.

## Report

`report.pdf` (4 pages) covers the generic cubic formulation and the parameter
assignments, the derivation of the departure functions and the fugacity
coefficient from `G^R`, the ideal-gas heat-capacity integrals, the numerical
method and root selection, the validation described above, and a discussion of
how the four equations compare.

## References

1. J. M. Smith, H. C. Van Ness, M. M. Abbott, M. T. Swihart, *Introduction to
   Chemical Engineering Thermodynamics*, 8th ed., McGraw-Hill — Table 3.1,
   Table B.1, Table C.1, Example 3.9, Example 13.5 and Table 13.7, Chapters 3, 6 and 13.
2. O. Redlich, J. N. S. Kwong, *Chem. Rev.*, **44**, 233 (1949).
3. G. Soave, *Chem. Eng. Sci.*, **27**, 1197 (1972).
4. D.-Y. Peng, D. B. Robinson, *Ind. Eng. Chem. Fundam.*, **15**, 59 (1976).
