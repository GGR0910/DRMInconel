@testset "Units & Physical Conversions" begin
    @test c_to_k(0.0) == 273.15
    @test c_to_k(800.0) == 1073.15
    @test k_to_c(1073.15) ≈ 800.0

    @test atm_to_bar(1.0) ≈ 1.01325
    @test bar_to_atm(1.01325) ≈ 1.0
    @test atm_to_pa(1.0) == 101325.0
    @test pa_to_atm(101325.0) == 1.0

    # Ideal gas molar volume at STP (273.15 K, 1 atm) = 22.413969545 L/mol
    # 22413.969545 mL/min must equal exactly 1/60 mol/s
    @test ml_min_flow_to_mol_s_flow(22413.969545; T_ref_K=273.15, P_ref_atm=1.0) ≈ (1.0 / 60.0) atol=1e-9

    # Exact round-trip conversion
    @test mol_s_flow_to_ml_min_flow(ml_min_flow_to_mol_s_flow(62.0)) ≈ 62.0 rtol=1e-12
    @test mol_s_flow_to_ml_min_flow(ml_min_flow_to_mol_s_flow(200.0; T_ref_K=298.15); T_K=298.15) ≈ 200.0 rtol=1e-12

    # Pure thermal expansion from STP (273.15 K) to 800 °C (1073.15 K)
    n_62 = ml_min_flow_to_mol_s_flow(62.0; T_ref_K=273.15, P_ref_atm=1.0)
    Q_hot = mol_s_flow_to_ml_min_flow(n_62; T_K=c_to_k(800.0), P_atm=1.0)
    @test Q_hot ≈ 62.0 * (1073.15 / 273.15) rtol=1e-12

    # Tube cross-section (2.0 cm ID) and gas velocity
    @test tube_cross_section_cm2(2.0) ≈ pi rtol=1e-12
    @test superficial_velocity_cm_s(60.0, 2.0) ≈ (1.0 / pi) rtol=1e-12
    @test superficial_velocity_cm_s(60.0, 2.0; porosity=0.5) ≈ (2.0 / pi) rtol=1e-12

    # Domain guards
    @test_throws DomainError ml_min_flow_to_mol_s_flow(-1.0)
    @test_throws DomainError tube_cross_section_cm2(0.0)
    @test_throws DomainError superficial_velocity_cm_s(10.0, 2.0; porosity=1.5)
end

@testset "Experimental Scenarios" begin
    s800 = scenario_screening_800()
    @test total_flow_ml_min(s800) ≈ 62.0
    frac800 = molar_fractions(s800)
    @test frac800[:CH4] ≈ 22.5 / 62.0
    @test frac800[:CO2] ≈ 22.5 / 62.0
    @test frac800[:N2]  ≈ 17.0 / 62.0
    @test sum(values(frac800)) ≈ 1.0

    s700 = scenario_stability_700()
    @test total_flow_ml_min(s700) ≈ 200.0
    frac700 = molar_fractions(s700)
    @test frac700[:CH4] ≈ 0.36
    @test frac700[:CO2] ≈ 0.36
    @test frac700[:N2]  ≈ 0.28

    custom = scenario_custom(CH4=30.0, CO2=60.0, N2=10.0, O2=5.0)
    @test total_flow_ml_min(custom) ≈ 105.0
    @test haskey(feed_to_mol_s_flow(custom), :O2)
end