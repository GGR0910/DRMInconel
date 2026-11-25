module DRMInconel

using JuMP
using Ipopt
using YAML
using OrdinaryDiffEq

include("types.jl")
include("units.jl")
include("scenarios.jl")
include("gas_properties.jl")
include("equilibrium.jl")
include("inconel600.jl")
include("monolith_geometry.jl")
include("induction_coupling.jl")
include("termica.jl")

# Termical
export SIGMA_SB, GRAVITY_G, INCONEL_EMISSIVITY_PRESETS
export ThermalReactorConfig, ChemistryThermalLoad
export air_properties, natural_convection_h, zero_chemistry_load, build_chemistry_load
export q_gas_W, q_reac_W, surface_heat_losses, steady_state_power_0d
export simulate_ramp_0d, effective_radial_k, simulate_radial_1d

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

# Inconel 600 properties
export INCONEL600_DENSITY_20C, INCONEL600_MU_R, INCONEL600_CURIE_C, INCONEL600_MELTING_RANGE_C
export inconel_resistivity, inconel_elec_conductivity, inconel_thermal_conductivity
export inconel_cp, inconel_cte, inconel_density, inconel_permeability_rel
export inconel_thermal_diffusivity, inconel_properties

# Monolith geometry & induction coupling
export AbstractMonolith, RolledFoilMonolith, HoneycombMonolith, MetalFoamMonolith, PrintedLatticeMonolith
export experimental_monolith, bulk_volume_m3, element_thickness_m, porosity, metal_volume_m3
export specific_surface_area, surface_area_m2, hydraulic_diameter_m, monolith_mass_kg, monolith_summary
export MU_0, RHO_COPPER_COIL, skin_depth_m, skin_depth_mm
export slab_power_factor, cylinder_power_factor, coil_electrical_efficiency
export effective_resistivity, effective_skin_depth_m, coupling_regime, analyze_monolith_coupling

end