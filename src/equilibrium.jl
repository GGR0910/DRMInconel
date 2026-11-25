const _EPS_GAS = 1.0e-12

"""
    solve_equilibrium(feed::FeedCondition; T_C, P_atm=feed.P_atm, mode=:base,
                      allow_coke=true, dG_carbon_J_mol=0.0,
                      geometry=ReactorGeometry(), db=default_thermo_db()) -> EquilibriumResult

Compute 0D isothermal chemical equilibrium at `(T_C, P_atm)` via Gibbs free energy
minimization (`JuMP` + `Ipopt`), returning molar outlet flows, conversions, H2/CO ratio,
carbon deposition, minimum thermal power demand (Watts), and outlet gas velocity.
"""
function solve_equilibrium(
    feed::FeedCondition;
    T_C::Real,
    P_atm::Real=feed.P_atm,
    mode::Symbol=:base,
    allow_coke::Bool=true,
    dG_carbon_J_mol::Real=0.0,
    geometry::ReactorGeometry=ReactorGeometry(),
    db::Dict{Symbol, SpeciesData}=default_thermo_db(),
)::EquilibriumResult
    T_K = c_to_k(T_C)
    T_in_K = c_to_k(feed.T_in_C)
    P_bar = atm_to_bar(P_atm)
    P_bar <= 0 && throw(DomainError(P_atm, "Reactor pressure must be positive."))

    # 1. Inlet molar flows (mol/s) and dimensionless 1-mol basis for solver scaling
    F_in = feed_to_mol_s_flow(feed)
    F_tot_in = sum(values(F_in))
    F_tot_in <= 0 && throw(DomainError(F_tot_in, "Total inlet molar flow must be positive."))

    # Elemental balance vector b_e on the normalized 1-mol inlet basis
    b_elem = Dict(:C => 0.0, :H => 0.0, :O => 0.0, :N => 0.0)
    for (sp, F_val) in F_in
        atoms = get_species_data(sp; db=db).atoms
        n_norm = F_val / F_tot_in
        b_elem[:C] += atoms.C * n_norm
        b_elem[:H] += atoms.H * n_norm
        b_elem[:O] += atoms.O * n_norm
        b_elem[:N] += atoms.N * n_norm
    end

    # 2. Select active species compatible with elements present in the feed
    all_species = species_list(mode; db=db)
    active_species = Symbol[]
    for sp in all_species
        if sp === :C_s && !allow_coke
            continue
        end
        atoms = get_species_data(sp; db=db).atoms
        present = (atoms.C == 0 || b_elem[:C] > 0) &&
                  (atoms.H == 0 || b_elem[:H] > 0) &&
                  (atoms.O == 0 || b_elem[:O] > 0) &&
                  (atoms.N == 0 || b_elem[:N] > 0)
        present && push!(active_species, sp)
    end

    gas_active = [sp for sp in active_species if get_species_data(sp; db=db).phase === :gas]
    has_solid_c = (:C_s in active_species)
    isempty(gas_active) && error("No active gas species available for the provided feed.")

    # 3. Dimensionless standard Gibbs free energies G°(T) / (R * T) at P° = 1 bar
    RT = R_GAS * T_K
    g_dimless = Dict(sp => gibbs_energy_J_mol(sp, T_K; db=db) / RT for sp in active_species)
    ln_P_ratio = log(P_bar)  # P° = 1.0 bar for NASA-9

    # 4. Build and solve nonlinear Gibbs minimization model
    model = Model(Ipopt.Optimizer)
    set_silent(model)
    set_attribute(model, "tol", 1.0e-10)
    set_attribute(model, "constr_viol_tol", 1.0e-12)
    set_attribute(model, "max_iter", 500)

    n_init_gas = 1.0 / length(gas_active)
    @variable(model, n[sp in active_species])
    for sp in gas_active
        set_lower_bound(n[sp], _EPS_GAS)
        set_start_value(n[sp], n_init_gas)
    end
    if has_solid_c
        set_lower_bound(n[:C_s], 0.0)
        set_start_value(n[:C_s], 0.05)
    end

    @variable(model, n_gas >= _EPS_GAS, start=1.0)
    @constraint(model, n_gas == sum(n[sp] for sp in gas_active))

    # Elemental conservation constraints: sum_i (a_{e,i} * n_i) == b_e
    for elem in (:C, :H, :O, :N)
        b_elem[elem] > 0 || continue
        @constraint(
            model,
            sum(getproperty(get_species_data(sp; db=db).atoms, elem) * n[sp] for sp in active_species) == b_elem[elem]
        )
    end

    # Objective: dimensionless total Gibbs free energy G_total / (R * T)
    gas_obj = sum(
        n[sp] * (g_dimless[sp] + ln_P_ratio + log(n[sp]) - log(n_gas))
        for sp in gas_active
    )
    if has_solid_c
        g_coke = g_dimless[:C_s] + Float64(dG_carbon_J_mol) / RT
        @objective(model, Min, gas_obj + n[:C_s] * g_coke)
    else
        @objective(model, Min, gas_obj)
    end

    optimize!(model)
    status = string(termination_status(model))
    if !(status in ("LOCALLY_SOLVED", "ALMOST_LOCALLY_SOLVED"))
        error("Gibbs minimization failed at T=$(T_C) °C, P=$(P_atm) atm with status: $status")
    end

    # 5. Rescale back to physical molar flow rates (mol/s)
    F_out = Dict{Symbol, Float64}()
    for sp in all_species
        if sp in active_species
            val = value(n[sp]) * F_tot_in
            F_out[sp] = val < 0.0 ? 0.0 : val
        else
            F_out[sp] = 0.0
        end
    end

    # 6. Performance metrics
    F_CH4_in = get(F_in, :CH4, 0.0)
    F_CO2_in = get(F_in, :CO2, 0.0)
    X_CH4 = F_CH4_in > 0 ? (F_CH4_in - get(F_out, :CH4, 0.0)) / F_CH4_in : NaN
    X_CO2 = F_CO2_in > 0 ? (F_CO2_in - get(F_out, :CO2, 0.0)) / F_CO2_in : NaN

    F_H2_out = get(F_out, :H2, 0.0)
    F_CO_out = get(F_out, :CO, 0.0)
    H2_CO_ratio = F_CO_out > 1.0e-15 ? (F_H2_out / F_CO_out) : NaN
    coke_mol_s = get(F_out, :C_s, 0.0)

    # 7. Minimum thermal power balance (W = J/s)
    H_in_at_Tin = sum(F * enthalpy_J_mol(sp, T_in_K; db=db) for (sp, F) in F_in)
    H_in_at_T   = sum(F * enthalpy_J_mol(sp, T_K; db=db) for (sp, F) in F_in)
    H_out_at_T  = sum(F * enthalpy_J_mol(sp, T_K; db=db) for (sp, F) in F_out if F > 0.0)

    Q_sensible_W = H_in_at_T - H_in_at_Tin
    Q_reaction_W = H_out_at_T - H_in_at_T
    Q_total_W    = Q_sensible_W + Q_reaction_W

    # 8. Outlet gas expansion, interstitial velocity, and residence time
    F_gas_out = sum(F for (sp, F) in F_out if get_species_data(sp; db=db).phase === :gas)
    V_out_ml_min = mol_s_flow_to_ml_min_flow(F_gas_out; T_K=T_K, P_atm=P_atm)
    expansion_factor = V_out_ml_min / total_flow_ml_min(feed)
    gas_velocity_cm_s = superficial_velocity_cm_s(
        V_out_ml_min,
        geometry.tube_id_cm;
        porosity=geometry.porosity,
    )
    residence_time_s = gas_velocity_cm_s > 0 ? (geometry.monolith_length_cm / gas_velocity_cm_s) : Inf

    return EquilibriumResult(
        Float64(T_C),
        Float64(P_atm),
        F_in,
        F_out,
        X_CH4,
        X_CO2,
        H2_CO_ratio,
        coke_mol_s,
        Q_sensible_W,
        Q_reaction_W,
        Q_total_W,
        V_out_ml_min,
        expansion_factor,
        gas_velocity_cm_s,
        residence_time_s,
        status,
    )
end