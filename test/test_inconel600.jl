@testset "Inconel 600 Material Properties (Special Metals SMC-027)" begin
    # 1. Constantes físicas de referência
    @test inconel_density(20.0) == 8470.0
    @test isapprox(inconel_permeability_rel(800.0), 1.010; atol=1e-6)
    @test INCONEL600_CURIE_C == -124.0

    # 2. Pontos exatos da Tabela 3 (Special Metals)
    ref_rho = [
        (20.0,  1.03e-6), (100.0, 1.04e-6), (200.0, 1.05e-6),
        (500.0, 1.12e-6), (700.0, 1.13e-6), (800.0, 1.13e-6), (900.0, 1.15e-6),
    ]
    for (T, expected) in ref_rho
        @test isapprox(inconel_resistivity(T), expected; rtol=1e-12)
    end

    ref_cp = [
        (20.0, 444.0), (100.0, 465.0), (300.0, 502.0),
        (500.0, 536.0), (700.0, 595.0), (800.0, 611.0), (900.0, 628.0),
    ]
    for (T, expected) in ref_cp
        @test isapprox(inconel_cp(T), expected; rtol=1e-12)
    end

    ref_k = [
        (20.0, 14.9), (100.0, 15.9), (400.0, 20.5),
        (600.0, 23.9), (700.0, 25.7), (800.0, 27.5),
    ]
    for (T, expected) in ref_k
        @test isapprox(inconel_thermal_conductivity(T), expected; rtol=1e-12)
    end

    # Extrapolação linear suave de k(T) a 900 °C (inclinação 1.8 W/(m·K) por 100 °C)
    @test isapprox(inconel_thermal_conductivity(900.0), 29.3; atol=1e-9)

    # 3. Monotonicidade física e consistência entre funções
    Ts = collect(25.0:25.0:900.0)
    @test all(diff(inconel_thermal_conductivity.(Ts)) .> 0.0)
    @test all(diff(inconel_cp.(Ts)) .> 0.0)
    @test all(diff(inconel_density.(Ts)) .< 0.0)      # densidade diminui ao dilatar
    @test all(diff(inconel_resistivity.(Ts)) .>= 0.0) # resistividade não-decrescente

    # 4. Verificação da massa real do monólito (2 folhas 3x30 cm x 50 um a 20 °C)
    V_metal_m3 = 2.0 * (0.03 * 0.30 * 50.0e-6)
    mass_g = inconel_density(20.0) * V_metal_m3 * 1000.0
    @test isapprox(mass_g, 7.623; atol=0.01)

    # 5. Consistência do agregador e limites de domínio
    p800 = inconel_properties(800.0)
    @test isapprox(p800.sigma_S_m * p800.rho_e_ohm_m, 1.0; atol=1e-12)
    @test isapprox(p800.alpha_m2_s, inconel_thermal_diffusivity(800.0); rtol=1e-12)
    @test_throws DomainError inconel_resistivity(0.0)
    @test_throws DomainError inconel_cp(1200.0)
end