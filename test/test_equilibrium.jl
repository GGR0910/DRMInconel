function _elemental_totals(flows::Dict{Symbol, Float64}, db)
    totals = Dict(:C => 0.0, :H => 0.0, :O => 0.0, :N => 0.0)
    for (sp, F) in flows
        atoms = get_species_data(sp; db=db).atoms
        totals[:C] += atoms.C * F
        totals[:H] += atoms.H * F
        totals[:O] += atoms.O * F
        totals[:N] += atoms.N * F
    end
    return totals
end

@testset "Gibbs Equilibrium Solver (Phase 1)" begin
    db = default_thermo_db()

    @testset "Atomic conservation & convergence across sweep of scenarios" begin
        feeds = [
            scenario_benchmark_1to1(),
            scenario_screening_800(),
            scenario_stability_700(),
            scenario_custom(CH4=40.0, CO2=80.0, N2=20.0),
            scenario_custom(CH4=50.0, CO2=50.0, N2=10.0, O2=5.0),
        ]
        for feed in feeds, T_C in (600.0, 700.0, 800.0, 900.0), allow_coke in (true, false)
            res = solve_equilibrium(feed; T_C=T_C, allow_coke=allow_coke, db=db)
            @test res.solver_status in ("LOCALLY_SOLVED", "ALMOST_LOCALLY_SOLVED")
            @test 0.0 <= res.X_CH4 <= 1.0
            @test 0.0 <= res.X_CO2 <= 1.0
            @test res.coke_mol_s >= 0.0

            atoms_in = _elemental_totals(res.n_in_mol_s, db)
            atoms_out = _elemental_totals(res.n_out_mol_s, db)
            for elem in (:C, :H, :O, :N)
                @test atoms_out[elem] ≈ atoms_in[elem] atol=1e-12 rtol=1e-9
            end
        end
    end

    @testset "Classical 1:1 DRM benchmark physics at 1 atm" begin
        bench = scenario_benchmark_1to1()

        # Without solid carbon (pure gas phase DRM + RWGS): X_CO2 > X_CH4 and H2/CO < 1
        for T_C in (600.0, 700.0, 800.0, 900.0)
            res = solve_equilibrium(bench; T_C=T_C, allow_coke=false, db=db)
            @test res.coke_mol_s == 0.0
            @test res.X_CO2 > res.X_CH4
            @test res.H2_CO_ratio < 1.0
        end

        # With ideal graphite in undiluted 1:1 at 1 atm: coke drops sharply from 600 to 900 °C,
        # leaving only a small residual (~2% of total C) at 900 °C
        r600 = solve_equilibrium(bench; T_C=600.0, allow_coke=true, db=db)
        r750 = solve_equilibrium(bench; T_C=750.0, allow_coke=true, db=db)
        r900 = solve_equilibrium(bench; T_C=900.0, allow_coke=true, db=db)

        F_C_total = r900.n_in_mol_s[:CH4] + r900.n_in_mol_s[:CO2]
        @test r600.coke_mol_s > r750.coke_mol_s > r900.coke_mol_s > 0.0
        @test (r900.coke_mol_s / F_C_total) < 0.05
        @test r900.X_CH4 > 0.95

        # With excess CO2 (1.5:1) at 900 °C: ideal graphite vanishes completely
        co2_rich = scenario_custom(CH4=40.0, CO2=60.0)
        r900_rich = solve_equilibrium(co2_rich; T_C=900.0, allow_coke=true, db=db)
        @test r900_rich.coke_mol_s < 1e-11
        @test r900_rich.H2_CO_ratio < 1.0

        # With filamentous carbon penalty (+10 kJ/mol) on 1:1 at 900 °C: coke vanishes
        # and gas-phase RWGS dominance (X_CO2 > X_CH4, H2/CO < 1) is restored
        r900_fil = solve_equilibrium(bench; T_C=900.0, allow_coke=true, dG_carbon_J_mol=10_000.0, db=db)
        @test r900_fil.coke_mol_s < 1e-11
        @test r900_fil.X_CO2 > r900_fil.X_CH4
        @test r900_fil.H2_CO_ratio < 1.0
    end
    
    @testset "Energy balance, flow-rate scaling & gas expansion" begin
        f1 = scenario_screening_800()
        f2 = scenario_custom(
            CH4 = 2.0 * f1.flow_ml_min[:CH4],
            CO2 = 2.0 * f1.flow_ml_min[:CO2],
            N2  = 2.0 * f1.flow_ml_min[:N2],
        )

        r1 = solve_equilibrium(f1; T_C=800.0, db=db)
        r2 = solve_equilibrium(f2; T_C=800.0, db=db)

        @test r1.Q_sensible_W > 0.0
        @test r1.Q_reaction_W > 0.0
        @test r1.Q_total_W ≈ (r1.Q_sensible_W + r1.Q_reaction_W) rtol=1e-12

        # 0D equilibrium conversion is invariant to flow rate; power and velocity scale linearly by 2x
        @test r2.X_CH4 ≈ r1.X_CH4 rtol=1e-6
        @test r2.X_CO2 ≈ r1.X_CO2 rtol=1e-6
        @test r2.Q_total_W ≈ 2.0 * r1.Q_total_W rtol=1e-6
        @test r2.gas_velocity_cm_s ≈ 2.0 * r1.gas_velocity_cm_s rtol=1e-6
        @test r2.residence_time_s ≈ 0.5 * r1.residence_time_s rtol=1e-6
        @test r1.expansion_factor > 4.0  # thermal + chemical expansion at 800 °C
    end

    @testset "Filamentous carbon penalty & O2 exploratory co-feeding" begin
        f700 = scenario_stability_700()
        r_graphite = solve_equilibrium(f700; T_C=700.0, allow_coke=true, dG_carbon_J_mol=0.0, db=db)
        r_filament = solve_equilibrium(f700; T_C=700.0, allow_coke=true, dG_carbon_J_mol=10_000.0, db=db)
        @test r_filament.coke_mol_s < r_graphite.coke_mol_s

        f_no_o2 = scenario_custom(CH4=50.0, CO2=50.0, N2=20.0)
        f_with_o2 = scenario_custom(CH4=50.0, CO2=50.0, N2=20.0, O2=10.0)
        r_no_o2 = solve_equilibrium(f_no_o2; T_C=750.0, allow_coke=true, db=db)
        r_with_o2 = solve_equilibrium(f_with_o2; T_C=750.0, allow_coke=true, db=db)
        @test r_with_o2.Q_reaction_W < r_no_o2.Q_reaction_W
        @test r_with_o2.coke_mol_s < r_no_o2.coke_mol_s
    end

    if all(haskey(db, sp) for sp in EXTENDED_SPECIES)
        @testset "Mode consistency: :base vs :extended (12 species)" begin
            f800 = scenario_screening_800()
            r_base = solve_equilibrium(f800; T_C=800.0, mode=:base, db=db)
            r_ext  = solve_equilibrium(f800; T_C=800.0, mode=:extended, db=db)

            @test r_ext.X_CH4 ≈ r_base.X_CH4 atol=1e-4
            @test r_ext.X_CO2 ≈ r_base.X_CO2 atol=1e-4
            @test r_ext.H2_CO_ratio ≈ r_base.H2_CO_ratio atol=1e-4

            F_tot = sum(values(r_ext.n_out_mol_s))
            for minor in (:C2H6, :C2H4, :CH3OH, :CH3OCH3)
                @test r_ext.n_out_mol_s[minor] / F_tot < 1e-4
            end
        end
    end
end

@testset "External Oracle Validation (Cantera 3.2.0, P_ref = 1 bar)" begin
    db = default_thermo_db()

    # screening_800 @ 800 °C — fase gasosa
    r_scr_gas = solve_equilibrium(scenario_screening_800(); T_C=800.0, allow_coke=false, db=db)
    @test isapprox(r_scr_gas.X_CH4,       0.91803; atol=1e-4)
    @test isapprox(r_scr_gas.X_CO2,       0.95653; atol=1e-4)
    @test isapprox(r_scr_gas.H2_CO_ratio, 0.9589;  atol=5e-4)
    @test isapprox(r_scr_gas.Q_total_W,   5.675;   atol=5e-3)

    # screening_800 @ 800 °C — com equilíbrio de grafite
    r_scr_coke = solve_equilibrium(scenario_screening_800(); T_C=800.0, allow_coke=true, db=db)
    @test isapprox(r_scr_coke.X_CH4,       0.96534; atol=1e-4)
    @test isapprox(r_scr_coke.X_CO2,       0.90696; atol=1e-4)
    @test isapprox(r_scr_coke.H2_CO_ratio, 1.0678;  atol=5e-4)
    @test isapprox(r_scr_coke.Q_total_W,   5.484;   atol=5e-3)

    # stability_700 @ 700 °C — fase gasosa
    r_stab_gas = solve_equilibrium(scenario_stability_700(); T_C=700.0, allow_coke=false, db=db)
    @test isapprox(r_stab_gas.X_CH4,       0.75843; atol=1e-4)
    @test isapprox(r_stab_gas.X_CO2,       0.84454; atol=1e-4)
    @test isapprox(r_stab_gas.H2_CO_ratio, 0.8926;  atol=5e-4)
    @test isapprox(r_stab_gas.Q_total_W,   15.241;  atol=5e-3)

    # stability_700 @ 700 °C — com equilíbrio de grafite
    r_stab_coke = solve_equilibrium(scenario_stability_700(); T_C=700.0, allow_coke=true, db=db)
    @test isapprox(r_stab_coke.X_CH4,       0.92179; atol=1e-4)
    @test isapprox(r_stab_coke.X_CO2,       0.69309; atol=1e-4)
    @test isapprox(r_stab_coke.H2_CO_ratio, 1.4092;  atol=5e-4)
    @test isapprox(r_stab_coke.Q_total_W,   13.307;  atol=5e-3)
end