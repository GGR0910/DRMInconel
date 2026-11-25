module DRMInconel

using JuMP
using Ipopt
using YAML

include("types.jl")
include("units.jl")
include("scenarios.jl")
include("gas_properties.jl")
include("equilibrium.jl")

# Types
export FeedCondition, ReactorGeometry, EquilibriumResult
export ThermoPolynomial, SpeciesData

# Physical constants
export R_GAS, ZERO_C_K, ATM_TO_PA, BAR_TO_PA, ATM_TO_BAR

# Units & geometry helpers
export c_to_k, k_to_c
export atm_to_bar, bar_to_atm, atm_to_pa, pa_to_atm
export ml_min_flow_to_mol_s_flow, mol_s_flow_to_ml_min_flow
export tube_cross_section_cm2, superficial_velocity_cm_s

# Scenarios & feed helpers
export total_flow_ml_min, molar_fractions, feed_to_mol_s_flow
export scenario_screening_800, scenario_stability_700, scenario_benchmark_1to1, scenario_custom

# Thermodynamic properties (Phase 1 — Cantera/NASA YAML)
export BASE_SPECIES, EXTENDED_SPECIES
export load_thermo_database, default_thermo_db, get_species_data, species_list
export heat_capacity_J_mol_K, sensible_enthalpy_J_mol
export enthalpy_J_mol, entropy_J_mol_K, gibbs_energy_J_mol

# Equilibrium solver (Phase 1)
export solve_equilibrium
end