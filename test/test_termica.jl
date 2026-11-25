@testset "Phase 3: Thermal Dynamics (0D) & Radial Profile (1D)" begin
    m = experimental_monolith()
    cfg_clean = ThermalReactorConfig(surface_state=:clean)
    cfg_ox    = ThermalReactorConfig(surface_state=:oxidized)
    cfg_coked = ThermalReactorConfig(surface_state=:coked)
    cfg_ins   = ThermalReactorConfig(surface_state=:oxidized, insulation_thickness_m=0.015)

    @testset "0D Steady-State Balance Closure (< 1e-6) & Emissivity/Insulation Physics" begin
        chem = build_chemistry_load(scenario_screening_800(); allow_coke=false)

        # Consistência com Fase 1: Q_gas + Q_reac == Q_total_W (~5.675 W a 800 °C)
        @test isapprox(q_gas_W(chem, 800.0) + q_reac_W(chem, 800.0), 5.675; atol=1e-2)
        @test q_gas_W(chem, 800.0) > 0.0
        @test q_reac_W(chem, 800.0) > q_gas_W(chem, 800.0)

        # Rampa 0D até regime permanente com potência exata para 800 °C
        bd800 = steady_state_power_0d(m, cfg_ox, 800.0; chem=chem)
        @test isapprox(bd800.pct_rad + bd800.pct_conv + bd800.pct_gas + bd800.pct_reac, 100.0; atol=1e-10)

        ramp = simulate_ramp_0d(m, cfg_ox; P_ind_W=bd800.P_ind_W, t_end_s=250.0, chem=chem)
        @test isapprox(ramp.T_final_C, 800.0; atol=1e-4)
        @test ramp.balance_rel_err < 1.0e-6
        @test all(diff(ramp.T_C) .>= -1e-9)

        # Efeito da emissividade: limpo < oxidado < coqueado na potência necessária a 800 °C
        p_clean = steady_state_power_0d(m, cfg_clean, 800.0).P_ind_W
        p_ox    = steady_state_power_0d(m, cfg_ox,    800.0).P_ind_W
        p_coked = steady_state_power_0d(m, cfg_coked, 800.0).P_ind_W
        @test p_clean < p_ox < p_coked

        # Isolamento de 15 mm reduz drasticamente as perdas térmicas
        p_ins = steady_state_power_0d(m, cfg_ins, 800.0).P_ind_W
        @test p_ins < 0.55 * p_ox
    end

    @testset "1D Radial Flux Conservation (< 1e-6) & Axis Symmetry (dT/dr = 0 at r = 0)" begin
        chem = build_chemistry_load(scenario_screening_800(); allow_coke=false)
        P_ref = steady_state_power_0d(m, cfg_ox, 800.0; chem=chem).P_ind_W

        # 1. Caso completo (com química + perdas axiais e laterais) a 200 kHz
        r200 = simulate_radial_1d(m, cfg_ox; P_ind_W=P_ref, f_hz=200e3, N=25, t_end_s=250.0, chem=chem)
        @test r200.flux_balance_rel_err < 1.0e-6
        @test 750.0 < r200.T_rim_C < 820.0
        @test r200.T_avg_C > r200.T_rim_C
        @test 750.0 < r200.T_avg_C < 920.0


        # 2. Simetria em r = 0: gradiente junto ao eixo é muito menor que na borda
        @test abs(r200.dT_dr_center_C_mm) < 0.25 * abs(r200.dT_dr_rim_C_mm)

        # 3. Refinamento de malha (N=20 -> N=40): dT/dr no centro converge para 0
        r_n20 = simulate_radial_1d(m, cfg_ox; P_ind_W=P_ref, f_hz=200e3, N=20, t_end_s=250.0, include_ends=false)
        r_n40 = simulate_radial_1d(m, cfg_ox; P_ind_W=P_ref, f_hz=200e3, N=40, t_end_s=250.0, include_ends=false)
        @test r_n20.flux_balance_rel_err < 1.0e-6
        @test r_n40.flux_balance_rel_err < 1.0e-6
        @test abs(r_n40.dT_dr_center_C_mm) < abs(r_n20.dT_dr_center_C_mm)

        # 4. Frequência maior (400 kHz) concentra mais calor junto à borda onde a radiação escapa,
        #    encurtando o caminho radial de condução e reduzindo o ΔT estacionário frente a 50 kHz
        r50  = simulate_radial_1d(m, cfg_ox; P_ind_W=P_ref, f_hz=50e3,  N=25, t_end_s=250.0, chem=chem)
        r400 = simulate_radial_1d(m, cfg_ox; P_ind_W=P_ref, f_hz=400e3, N=25, t_end_s=250.0, chem=chem)
        @test r400.q_ind_MW_m3[end] > r50.q_ind_MW_m3[end]
        @test r400.delta_T_max_min_C < r50.delta_T_max_min_C
    end
end