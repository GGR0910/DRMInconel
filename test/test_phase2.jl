@testset "Phase 2: Monolith Geometry & Induction Coupling" begin
    @testset "Experimental Rolled Foil Monolith Audit" begin
        m = experimental_monolith()
        s = monolith_summary(m; T_C=20.0)

        # Bate exatamente com os dados de laboratório (Goodfellow 2x 3x30 cm, 50 um, 7.62 g)
        @test isapprox(s.mass_g, 7.623; atol=0.01)
        @test isapprox(s.surface_area_cm2, 360.0; atol=1e-9)
        @test isapprox(s.porosity_pct, 90.4507; atol=1e-3)
        @test isapprox(s.Sv_m2_m3, 3819.72; atol=0.05)
        @test isapprox(s.hydraulic_diam_mm, 0.9472; atol=1e-3)
    end

    @testset "Alternative Monolith Geometries Consistency" begin
        hc = HoneycombMonolith(cpsi=400.0, wall_thickness_m=150e-6)
        pitch = 0.0254 / sqrt(400.0)
        # Em dutos quadrados, dh == largura livre do canal (pitch - wall)
        @test isapprox(hydraulic_diameter_m(hc), pitch - 150e-6; rtol=1e-12)
        @test 0.5 < porosity(hc) < 0.9

        foam = MetalFoamMonolith(ppi=30.0, porosity=0.90)
        @test isapprox(porosity(foam), 0.90; atol=1e-12)
        @test 0.05e-3 < element_thickness_m(foam) < 0.5e-3
        @test 1000.0 < specific_surface_area(foam) < 8000.0

        lat = PrintedLatticeMonolith(cell_size_m=2.0e-3, strut_diameter_m=0.4e-3)
        @test 0.85 < porosity(lat) < 0.95
        @test isapprox(hydraulic_diameter_m(lat), 4.0 * porosity(lat) / specific_surface_area(lat); rtol=1e-12)
    end

    @testset "Skin Depth Scaling Laws" begin
        # Inconel 600 a 800 °C (rho = 1.13e-6 Ohm*m, mu_r = 1.010)
        d_100 = skin_depth_mm(100e3, 800.0)
        d_400 = skin_depth_mm(400e3, 800.0)
        d_10  = skin_depth_mm(10e3,  800.0)

        @test isapprox(d_100, 1.6834; atol=1e-3)
        # Escala exata com 1 / sqrt(f): quadruplicar f corta delta exatamente pela metade
        @test isapprox(d_100 / d_400, 2.0; rtol=1e-12)
        @test isapprox(d_10 / d_100, sqrt(10.0); rtol=1e-12)
        # Resistividade sobe com T -> skin depth sobe de 25 °C para 900 °C
        @test skin_depth_mm(100e3, 900.0) > skin_depth_mm(100e3, 25.0)
    end

    @testset "Analytical Inductive Coupling Limits (Slab & Cylinder)" begin
        # 1. Placa fina: F_slab(xi) -> xi^3 / 6 para xi << 1
        for xi in (1e-4, 1e-2, 0.05)
            @test isapprox(slab_power_factor(xi), xi^3 / 6.0; rtol=1e-3)
        end
        # Placa espessa: F_slab -> 1.0; overshoot clássico em xi = pi (~1.0903)
        @test isapprox(slab_power_factor(50.0), 1.0; atol=1e-12)
        @test isapprox(slab_power_factor(pi), 1.09033; atol=1e-4)

        # 2. Cilindro fino: F_cyl(zeta) -> zeta^3 / 32 para zeta << 1
        for zeta in (1e-4, 1e-2, 0.1)
            @test isapprox(cylinder_power_factor(zeta), zeta^3 / 32.0; rtol=2e-3)
        end
        # Cilindro espesso: F_cyl(zeta) ~ 1 - 1/zeta
        @test isapprox(cylinder_power_factor(100.0), 0.99; atol=2e-4)
        @test isapprox(cylinder_power_factor(1000.0), 0.999; atol=2e-5)

        # 3. Continuidade suave na fronteira de algoritmo (|z| = 35 <=> zeta ≈ 49.5)
        @test abs(cylinder_power_factor(49.4) - cylinder_power_factor(49.6)) < 2e-4
    end

    @testset "50 um Foil Paradox (Isolated Foil vs. Macroscopic Rolled Monolith)" begin
        m = experimental_monolith(contact_factor=0.5)
        res = analyze_monolith_coupling(m, 200e3, 800.0)

        # Folha isolada de 50 um é praticamente transparente a 200 kHz
        @test res.wall_to_delta < 0.05
        @test res.F_wall_isolated < 1e-4
        @test res.eta_wall_pct < 0.1

        # Enrolado com contato elétrico entre voltas, o cilindro de 20 mm acopla forte
        @test res.D_to_delta_eff > 3.5
        @test res.F_macro_cyl > 0.65
        @test res.eta_macro_pct > 85.0
        @test res.macro_regime === :optimal
    end
end