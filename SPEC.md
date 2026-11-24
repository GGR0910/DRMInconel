# DRMInconel — Technical Specification (Phase 0 & Phase 1)

## 1. Project Scope
Macroscopic scientific simulator (Level 1) in Julia for Dry Reforming of Methane (DRM) over a self-supported structured **Inconel 600** monolith heated by electromagnetic induction inside a quartz tube reactor.

- **Primary reaction:** CH₄ + CO₂ ⇌ 2 H₂ + 2 CO (ΔH°₂₉₈ = +247 kJ/mol)
- **Phase 1 focus:** 0D thermodynamic equilibrium (Gibbs free energy minimization), carbon deposition boundary, gas expansion, and minimum thermal power demand (Watts).

---

## 2. Experimental Baseline Parameters

### 2.1 Feed Scenarios (Phase 1)
All scenarios are parameterized; custom ratios and flow rates are supported directly by the solver.

| Scenario | T (°C) | Feed Recipe (mL/min) | Molar Composition (%) | Total Flow (mL/min) | Experimental Reference |
| :--- | :---: | :--- | :--- | :---: | :--- |
| **Screening (`screening_800`)** | 800 | 50 (CH₄:CO₂:N₂ = 45:45:10) + 12 N₂ | 36.29% CH₄, 36.29% CO₂, 27.42% N₂ | 62.0 | MON-Ind-800-Air: >90% CH₄, >95% CO₂; MON-Conv-800-Air: ~80% CH₄, ~90% CO₂ |
| **Stability (`stability_700`)** | 700 | 160 (CH₄:CO₂:N₂ = 45:45:10) + 40 N₂ | 36.00% CH₄, 36.00% CO₂, 28.00% N₂ | 200.0 | 24 h TOS: H₂/CO ≈ 0.8; CO₂: 80% → 55%, CH₄: 60% → 39% |
| **Literature Benchmark (`benchmark_1to1`)** | 600–900 | 50 CH₄ + 50 CO₂ | 50.0% CH₄, 50.0% CO₂, 0.0% N₂ | 100.0 | Classical 1:1 DRM at 1 atm |

- **Role of O₂:** Not co-fed during DRM. Atmospheric O₂ (1 atm) is used exclusively during thermal surface activation (800 °C, 2 h or 15 min induction) prior to *in situ* H₂ reduction (800 °C, 2 h, 50 mL/min pure H₂). Retained in the thermodynamic database as an optional species for exploratory studies.
- **Temperature sweep:** 600 °C to 900 °C (step 10 °C), covering the 700 °C stability test, 800 °C screening, and the 750–850 °C operational window.

### 2.2 Configurable Fallbacks
- **Reactor pressure (`P_atm`):** `1.0 atm` (`1.01325 bar`)
- **Gas inlet temperature (`T_in_C`):** `25.0 °C`
- **Standard flow reference (`T_ref_K`):** `273.15 K` (STP; configurable to `298.15 K` if MFC calibration uses 25 °C)

### 2.3 Reactor & Monolith Geometry (Used for velocity in Phase 1; core of Phases 2–4)
- **Quartz tube:** `2.0 cm` inner diameter (configurable), `40.0 cm` height, uninsulated viewing zone (monitored by thermal camera with electrical power control).
- **Inconel 600 monolith:** 1 flat + 1 corrugated commercial foil (Goodfellow), `3 × 30 cm` each, `50 µm` thickness, `360 cm²` total geometric surface area, no washcoat (active Ni⁰ via exsolution). Estimated open porosity: `~0.88` (configurable).

---

## 3. Thermodynamic Model & Species Modes (Phase 1)
- **Assumptions:** Ideal gas phase, pure solid carbon phase (activity = 1), isothermal 0D equilibrium, zero heat losses (minimum theoretical energy floor).
- **Species modes:**
  - `:base` (8 species): `CH4`, `CO2`, `H2`, `CO`, `H2O`, `C_s`, `N2`, `O2`. Automatically captures DRM, RWGS, Boudouard, CH₄ cracking, CO/CO₂ hydrogenation to solid carbon, and CO/CO₂ methanation.
  - `:extended` (12 species): adds `C2H6`, `C2H4`, `CH3OH`, `CH3OCH3` (DME) to evaluate all 17 side reactions (oxidative coupling, methanol synthesis, DME routes).
- **Solid carbon tuning:** Optional `allow_coke::Bool` flag and `dG_carbon_J_mol::Float64` offset to compare ideal graphite (`0.0`) against filamentous carbon on Ni catalysts.

---

## 4. Validation Criteria (`test/runtests.jl`)
1. **Unit consistency:** Exact round-trip conversions (`°C ↔ K`, `atm ↔ bar ↔ Pa`, `mL/min ↔ mol/s`) and standard molar volume check (`22.414 L/mol` at STP).
2. **Thermodynamic data accuracy:** Shomate polynomials verified directly against NIST Chemistry WebBook tables at `298.15 K` and `1000 K`.
3. **Physical invariants:** Atomic conservation (C, H, O, N) with error `< 1e-10 mol` and optimal solver termination (`LOCALLY_SOLVED`) across all runs.
4. **Literature & external cross-check:** Reproduce 1:1 DRM benchmark (`X_CO2 > X_CH4`, `H2/CO < 1.0`, coke extinction near 850–900 °C) and cross-check 3 temperatures (`600`, `750`, `900 °C`) against **NASA CEA** (`< 0.5%` absolute deviation).
5. **Mode consistency:** Confirm that `:extended` species (`C2H6`, `C2H4`, `CH3OH`, `CH3OCH3`) remain negligible at 1 atm and 700–850 °C relative to `:base`.