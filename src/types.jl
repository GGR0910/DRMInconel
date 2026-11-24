"""
    FeedCondition

Inlet gas stream definition in laboratory units.
`flow_ml_min` maps species symbols (e.g. `:CH4`, `:CO2`, `:N2`, `:O2`)
to volumetric flow rates in mL/min referenced at `(T_ref_K, P_ref_atm)`.
"""
Base.@kwdef struct FeedCondition
    name::String = "custom"
    flow_ml_min::Dict{Symbol, Float64}
    T_in_C::Float64 = 25.0
    P_atm::Float64 = 1.0
    T_ref_K::Float64 = 273.15
    P_ref_atm::Float64 = 1.0
end

"""
    ReactorGeometry

Physical dimensions of the quartz tube and Inconel 600 monolith.
Used in Phase 1 for gas velocity/residence time estimates and in Phases 2–4
for electromagnetic and thermal modeling.
"""
Base.@kwdef struct ReactorGeometry
    tube_id_cm::Float64 = 2.0
    tube_length_cm::Float64 = 40.0
    monolith_length_cm::Float64 = 3.0
    foil_thickness_um::Float64 = 50.0
    geometric_area_cm2::Float64 = 360.0
    porosity::Float64 = 0.88
end

"""
    EquilibriumResult

Output container returned by the Phase 1 Gibbs minimization solver.
"""
struct EquilibriumResult
    T_C::Float64
    P_atm::Float64
    n_in_mol_s::Dict{Symbol, Float64}
    n_out_mol_s::Dict{Symbol, Float64}
    X_CH4::Float64
    X_CO2::Float64
    H2_CO_ratio::Float64
    coke_mol_s::Float64
    Q_sensible_W::Float64
    Q_reaction_W::Float64
    Q_total_W::Float64
    V_out_ml_min::Float64
    expansion_factor::Float64
    gas_velocity_cm_s::Float64
    residence_time_s::Float64
    solver_status::String
end