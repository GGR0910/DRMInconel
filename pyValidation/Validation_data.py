import csv
import json
from pathlib import Path
import cantera as ct

CSV_PATH = Path(__file__).resolve().parent.parent / "output" / "equilibrium_reference.csv"

SCENARIOS = {
    "screening_800":  {"CH4": 22.5, "CO2": 22.5, "N2": 17.0},
    "stability_700":  {"CH4": 72.0, "CO2": 72.0, "N2": 56.0},
    "benchmark_1to1": {"CH4": 50.0, "CO2": 50.0},
}

GAS_SPECIES_NAMES = ["CH4", "CO2", "H2", "CO", "H2O", "N2", "O2"]
R_GAS = 8.31446261815324
T_STP_K = 273.15
ATM_TO_PA = 101325.0


def build_cantera_phases(ref_pressure_pa=100000.0):
    """Build independent gas and solid-graphite phases with explicit reference pressure."""
    all_gas = {s.name: s for s in ct.Species.list_from_file("nasa_gas.yaml")}
    gas_species = []
    for name in GAS_SPECIES_NAMES:
        sp_dict = dict(all_gas[name].input_data)
        sp_dict["thermo"] = dict(sp_dict["thermo"])
        sp_dict["thermo"]["reference-pressure"] = ref_pressure_pa
        gas_species.append(sp_dict)

    gas_yaml = json.dumps({
        "phases": [{
            "name": "gas",
            "thermo": "ideal-gas",
            "species": GAS_SPECIES_NAMES,
        }],
        "species": gas_species,
    })
    gas = ct.Solution(yaml=gas_yaml)

    all_cond = {s.name: s for s in ct.Species.list_from_file("nasa_condensed.yaml")}
    c_dict = dict(all_cond["C(gr)"].input_data)
    c_dict["thermo"] = dict(c_dict["thermo"])
    c_dict["thermo"]["reference-pressure"] = ref_pressure_pa
    c_dict["equation-of-state"] = {"model": "constant-volume", "density": 2260.0}
    graphite_yaml = json.dumps({
        "phases": [{
            "name": "graphite",
            "thermo": "fixed-stoichiometry",
            "species": ["C(gr)"],
        }],
        "species": [c_dict],
    })
    graphite = ct.Solution(yaml=graphite_yaml)
    return gas, graphite

def solve_cantera(gas, graphite, feed, T_C, P_atm, V_in_ml_min, allow_coke):
    T_K = T_C + 273.15
    P_Pa = P_atm * ct.one_atm

    total_feed = sum(feed.values())
    n_ch4_in = feed["CH4"] / total_feed
    n_co2_in = feed["CO2"] / total_feed

    # Inlet enthalpy at 25 °C (298.15 K) on a 1-mol basis
    gas.TPX = 298.15, P_Pa, feed
    h_in_J_mol = gas.enthalpy_mole / 1000.0
    m_in = gas.mean_molecular_weight

    # Equilibrate at (T_K, P_Pa)
    gas.TPX = T_K, P_Pa, feed
    if allow_coke:
        graphite.TP = T_K, P_Pa
        mix = ct.Mixture([(gas, 1.0), (graphite, 0.0)])
        mix.T = T_K
        mix.P = P_Pa
        mix.equilibrate("TP")
        n_gas_out = float(mix.phase_moles(0))
        n_coke_out = float(mix.phase_moles(1))
    else:
        gas.equilibrate("TP")
        n_gas_out = float(m_in / gas.mean_molecular_weight)
        n_coke_out = 0.0

    def n_out(sp):
        return n_gas_out * float(gas.X[gas.species_index(sp)])

    X_CH4_pct = 100.0 * (n_ch4_in - n_out("CH4")) / n_ch4_in
    X_CO2_pct = 100.0 * (n_co2_in - n_out("CO2")) / n_co2_in
    H2_CO = n_out("H2") / n_out("CO")
    coke_yield = n_coke_out / n_ch4_in

    # Thermal power (Watts)
    h_out_J_mol = (
        n_gas_out * (gas.enthalpy_mole / 1000.0)
        + n_coke_out * (graphite.enthalpy_mole / 1000.0)
    )
    F_tot_in = (ATM_TO_PA * (V_in_ml_min * 1e-6 / 60.0)) / (R_GAS * T_STP_K)
    Q_total_W = F_tot_in * (h_out_J_mol - h_in_J_mol)

    return {
        "X_CH4_pct": X_CH4_pct,
        "X_CO2_pct": X_CO2_pct,
        "H2_CO_ratio": H2_CO,
        "coke_mol_per_mol_CH4_in": coke_yield,
        "Q_total_W": Q_total_W,
    }


def main():
    if not CSV_PATH.exists():
        raise FileNotFoundError(
            f"Arquivo nao encontrado: {CSV_PATH}\n"
            "Rode antes: julia --project=. scripts/run_phase1.jl"
        )

    gas, graphite = build_cantera_phases()
    metrics = [
        ("X_CH4_pct", 0.05, "p.p."),
        ("X_CO2_pct", 0.05, "p.p."),
        ("H2_CO_ratio", 0.001, "-"),
        ("coke_mol_per_mol_CH4_in", 0.001, "mol/mol"),
        ("Q_total_W", 0.05, "W"),
    ]
    max_err = {name: (0.0, "") for name, _, _ in metrics}
    spot_checks = []

    with CSV_PATH.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))

    for row in rows:
        scen = row["scenario"]
        allow_coke = row["coke_model"] == "graphite_equilibrium"
        T_C = float(row["T_C"])
        P_atm = float(row["P_atm"])
        V_in = float(row["V_in_ml_min"])

        ct_res = solve_cantera(
            gas, graphite, SCENARIOS[scen], T_C, P_atm, V_in, allow_coke
        )

        tag = f"{scen} | {row['coke_model']} | {T_C:.0f} C"
        for name, _, _ in metrics:
            err = abs(float(row[name]) - ct_res[name])
            if err > max_err[name][0]:
                max_err[name] = (err, tag)

        if (scen == "screening_800" and T_C == 800.0) or (scen == "stability_700" and T_C == 700.0):
            spot_checks.append((tag, row, ct_res))

    print(f"Cantera version: {ct.__version__}")
    print(f"Total equilibrium points cross-checked: {len(rows)}\n")

    print("=== Maximum Absolute Deviation (Julia DRMInconel vs. Cantera) ===")
    all_passed = True
    for name, tol, unit in metrics:
        err, worst_case = max_err[name]
        ok = err <= tol
        all_passed = all_passed and ok
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {name:<25} max_err = {err:.6f} {unit:<7} (tol <= {tol})  @ {worst_case}")

    print("\n=== Experimental Operating Points (Julia vs. Cantera) ===")
    for tag, j_row, c_res in spot_checks:
        print(f"\n  {tag}")
        print(
            f"    Julia  : X_CH4={float(j_row['X_CH4_pct']):.3f}% | "
            f"X_CO2={float(j_row['X_CO2_pct']):.3f}% | "
            f"H2/CO={float(j_row['H2_CO_ratio']):.4f} | "
            f"Q={float(j_row['Q_total_W']):.3f} W"
        )
        print(
            f"    Cantera: X_CH4={c_res['X_CH4_pct']:.3f}% | "
            f"X_CO2={c_res['X_CO2_pct']:.3f}% | "
            f"H2/CO={c_res['H2_CO_ratio']:.4f} | "
            f"Q={c_res['Q_total_W']:.3f} W"
        )

    print("\nVERDICT:", "VALIDATED (all 186 points match Cantera)" if all_passed else "CHECK DISCREPANCIES ABOVE")


if __name__ == "__main__":
    main()