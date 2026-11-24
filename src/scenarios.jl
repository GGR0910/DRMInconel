"""
    total_flow_ml_min(feed::FeedCondition) -> Float64

Return the sum of all species inlet flow rates (mL/min).
"""
total_flow_ml_min(feed::FeedCondition)::Float64 = sum(values(feed.flow_ml_min); init=0.0)

"""
    molar_fractions(feed::FeedCondition) -> Dict{Symbol, Float64}

Return inlet molar (volume) fractions of each species in `feed`.
"""
function molar_fractions(feed::FeedCondition)::Dict{Symbol, Float64}
    total = total_flow_ml_min(feed)
    total <= 0 && throw(DomainError(total, "Total feed flow must be positive."))
    return Dict(sp => q / total for (sp, q) in feed.flow_ml_min)
end

"""
    feed_to_mols(feed::FeedCondition) -> Dict{Symbol, Float64}

Convert all inlet volumetric flows in `feed` to molar flow rates (mol/s).
"""
function feed_to_mol_s_flow(feed::FeedCondition)::Dict{Symbol, Float64}
    return Dict(
        sp => ml_min_flow_to_mol_s_flow(q; T_ref_K=feed.T_ref_K, P_ref_atm=feed.P_ref_atm)
        for (sp, q) in feed.flow_ml_min
    )
end

"""
    scenario_screening_800(; T_in_C=25.0, P_atm=1.0, T_ref_K=273.15)

Experimental screening condition at 800 °C:
50 mL/min reactant mix (CH4:CO2:N2 = 45:45:10) + 12 mL/min N2 = 62 mL/min total.
"""
function scenario_screening_800(; T_in_C=25.0, P_atm=1.0, T_ref_K=ZERO_C_K)
    return FeedCondition(
        name="screening_800",
        flow_ml_min=Dict(
            :CH4 => 50.0 * 0.45,
            :CO2 => 50.0 * 0.45,
            :N2  => 50.0 * 0.10 + 12.0,
        ),
        T_in_C=T_in_C,
        P_atm=P_atm,
        T_ref_K=T_ref_K,
    )
end

"""
    scenario_stability_700(; T_in_C=25.0, P_atm=1.0, T_ref_K=273.15)

Experimental 24 h stability condition at 700 °C:
160 mL/min reactant mix (CH4:CO2:N2 = 45:45:10) + 40 mL/min N2 = 200 mL/min total.
"""
function scenario_stability_700(; T_in_C=25.0, P_atm=1.0, T_ref_K=ZERO_C_K)
    return FeedCondition(
        name="stability_700",
        flow_ml_min=Dict(
            :CH4 => 160.0 * 0.45,
            :CO2 => 160.0 * 0.45,
            :N2  => 160.0 * 0.10 + 40.0,
        ),
        T_in_C=T_in_C,
        P_atm=P_atm,
        T_ref_K=T_ref_K,
    )
end

"""
    scenario_benchmark_1to1(; total_ml_min=100.0, T_in_C=25.0, P_atm=1.0, T_ref_K=273.15)

Classical equimolar 1:1 CH4:CO2 literature benchmark (no inert gas).
"""
function scenario_benchmark_1to1(; total_ml_min=100.0, T_in_C=25.0, P_atm=1.0, T_ref_K=ZERO_C_K)
    half = Float64(total_ml_min) / 2.0
    return FeedCondition(
        name="benchmark_1to1",
        flow_ml_min=Dict(:CH4 => half, :CO2 => half),
        T_in_C=T_in_C,
        P_atm=P_atm,
        T_ref_K=T_ref_K,
    )
end

"""
    scenario_custom(; CH4=50.0, CO2=50.0, N2=0.0, O2=0.0, name="custom", T_in_C=25.0, P_atm=1.0, T_ref_K=273.15)

Helper to build arbitrary feed mixtures and flow rates for parametric sweeps.
"""
function scenario_custom(;
    CH4::Real=50.0,
    CO2::Real=50.0,
    N2::Real=0.0,
    O2::Real=0.0,
    name::String="custom",
    T_in_C::Real=25.0,
    P_atm::Real=1.0,
    T_ref_K::Real=ZERO_C_K,
)
    flows = Dict{Symbol, Float64}()
    CH4 > 0 && (flows[:CH4] = Float64(CH4))
    CO2 > 0 && (flows[:CO2] = Float64(CO2))
    N2  > 0 && (flows[:N2]  = Float64(N2))
    O2  > 0 && (flows[:O2]  = Float64(O2))
    isempty(flows) && throw(ArgumentError("At least one species must have positive flow."))

    return FeedCondition(
        name=name,
        flow_ml_min=flows,
        T_in_C=Float64(T_in_C),
        P_atm=Float64(P_atm),
        T_ref_K=Float64(T_ref_K),
    )
end