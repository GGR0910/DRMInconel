# Physical constants (CODATA / SI exact definitions)
const R_GAS = 8.31446261815324   # J / (mol * K)
const ZERO_C_K = 273.15          # K
const ATM_TO_PA = 101325.0       # Pa / atm
const BAR_TO_PA = 1.0e5          # Pa / bar
const ATM_TO_BAR = ATM_TO_PA / BAR_TO_PA  # 1.01325 bar / atm

# Temperature conversions
c_to_k(T_C::Real) = Float64(T_C) + ZERO_C_K
k_to_c(T_K::Real) = Float64(T_K) - ZERO_C_K

# Pressure conversions
atm_to_bar(P_atm::Real) = Float64(P_atm) * ATM_TO_BAR
bar_to_atm(P_bar::Real) = Float64(P_bar) / ATM_TO_BAR
atm_to_pa(P_atm::Real)  = Float64(P_atm) * ATM_TO_PA
pa_to_atm(P_pa::Real)   = Float64(P_pa)  / ATM_TO_PA

"""
    ml_min_flow_to_mol_s_flow(vol_flow_ml_min; T_ref_K=273.15, P_ref_atm=1.0) -> Float64

Convert a gas volumetric flow rate in mL/min (SCCM) at reference state
`(T_ref_K, P_ref_atm)` into molar flow rate in mol/s using the ideal gas law.
"""
function ml_min_flow_to_mol_s_flow(vol_flow_ml_min::Real; T_ref_K::Real=ZERO_C_K, P_ref_atm::Real=1.0)::Float64
    vol_flow_ml_min < 0 && throw(DomainError(vol_flow_ml_min, "Volumetric flow rate must be non-negative."))
    T_ref_K <= 0 && throw(DomainError(T_ref_K, "Reference temperature (K) must be positive."))
    P_ref_atm <= 0 && throw(DomainError(P_ref_atm, "Reference pressure (atm) must be positive."))

    V_dot_m3_s = (Float64(vol_flow_ml_min) * 1.0e-6) / 60.0
    P_pa = atm_to_pa(P_ref_atm)
    return (P_pa * V_dot_m3_s) / (R_GAS * Float64(T_ref_K))
end

"""
    mol_s_flow_to_ml_min_flow(n_mol_s; T_K=273.15, P_atm=1.0) -> Float64

Convert a molar flow rate in mol/s into actual volumetric flow rate in mL/min
at operating state `(T_K, P_atm)` using the ideal gas law.
"""
function mol_s_flow_to_ml_min_flow(n_mol_s::Real; T_K::Real=ZERO_C_K, P_atm::Real=1.0)::Float64
    n_mol_s < 0 && throw(DomainError(n_mol_s, "Molar flow rate must be non-negative."))
    T_K <= 0 && throw(DomainError(T_K, "Temperature (K) must be positive."))
    P_atm <= 0 && throw(DomainError(P_atm, "Pressure (atm) must be positive."))

    P_pa = atm_to_pa(P_atm)
    V_dot_m3_s = (Float64(n_mol_s) * R_GAS * Float64(T_K)) / P_pa
    return V_dot_m3_s * 1.0e6 * 60.0
end

"""
    tube_cross_section_cm2(diameter_cm) -> Float64

Cross-sectional area (cm²) of a circular tube.
"""
function tube_cross_section_cm2(diameter_cm::Real)::Float64
    diameter_cm <= 0 && throw(DomainError(diameter_cm, "Tube diameter must be positive."))
    return pi * (Float64(diameter_cm) / 2.0)^2
end

"""
    superficial_velocity_cm_s(vol_flow_ml_min, diameter_cm; porosity=1.0) -> Float64

Calculate gas linear velocity (cm/s) inside a tube of `diameter_cm`.
Use `porosity = 1.0` for empty-tube superficial velocity or `porosity < 1.0`
for interstitial velocity inside the monolith channels.
"""
function superficial_velocity_cm_s(
    vol_flow_ml_min::Real,
    diameter_cm::Real;
    porosity::Real=1.0,
)::Float64
    vol_flow_ml_min < 0 && throw(DomainError(vol_flow_ml_min, "Flow rate must be non-negative."))
    (porosity <= 0 || porosity > 1) && throw(DomainError(porosity, "Porosity must be in (0, 1]."))

    area_cm2 = tube_cross_section_cm2(diameter_cm) * Float64(porosity)
    flow_rate_cm3_s = Float64(vol_flow_ml_min) / 60.0
    return flow_rate_cm3_s / area_cm2
end