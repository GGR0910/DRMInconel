@testset "Thermodynamic Properties (Cantera / NASA-9 YAML)" begin
    db = default_thermo_db()

    @testset "Database loading & stoichiometry" begin
        for sp in BASE_SPECIES
            @test haskey(db, sp)
        end
        @test db[:CH4].atoms == (C=1, H=4, O=0, N=0)
        @test db[:CO2].atoms == (C=1, H=0, O=2, N=0)
        @test db[:H2O].atoms == (C=0, H=2, O=1, N=0)
        @test db[:C_s].atoms == (C=1, H=0, O=0, N=0)
        @test db[:C_s].phase === :solid
        @test db[:CH4].phase === :gas
    end

    @testset "Exact thermodynamic calculus identities (dH/dT = Cp, T*dS/dT = Cp, dG/dT = -S)" begin
        # Independent mathematical proof that the NASA-9 / NASA-7 formulas have no
        # sign, exponent, or unit errors — without relying on hardcoded reference tables.
        h = 1.0e-3
        test_temps = [350.0, 700.0, 950.0, 1073.15, 1300.0]

        for sp in keys(db)
            for T in test_temps
                cp = heat_capacity_J_mol_K(sp, T; db=db)
                s  = entropy_J_mol_K(sp, T; db=db)

                dH_dT = (enthalpy_J_mol(sp, T + h; db=db) - enthalpy_J_mol(sp, T - h; db=db)) / (2.0 * h)
                dS_dT = (entropy_J_mol_K(sp, T + h; db=db) - entropy_J_mol_K(sp, T - h; db=db)) / (2.0 * h)
                dG_dT = (gibbs_energy_J_mol(sp, T + h; db=db) - gibbs_energy_J_mol(sp, T - h; db=db)) / (2.0 * h)

                @test dH_dT ≈ cp rtol=1e-6
                @test T * dS_dT ≈ cp rtol=1e-6
                @test dG_dT ≈ -s rtol=1e-6
            end
        end
    end

    @testset "Piecewise continuity across YAML temperature boundaries" begin
        eps_T = 1.0e-3
        for sp in keys(db)
            ranges = db[sp].ranges
            for i in 1:(length(ranges) - 1)
                T_split = ranges[i].T_max_K
                @test heat_capacity_J_mol_K(sp, T_split - eps_T; db=db) ≈
                      heat_capacity_J_mol_K(sp, T_split + eps_T; db=db) atol=0.1
                @test enthalpy_J_mol(sp, T_split - eps_T; db=db) ≈
                      enthalpy_J_mol(sp, T_split + eps_T; db=db) atol=50.0
                @test entropy_J_mol_K(sp, T_split - eps_T; db=db) ≈
                      entropy_J_mol_K(sp, T_split + eps_T; db=db) atol=0.1
            end
        end
    end

    @testset "Reference state elements & reaction enthalpies vs. qualification tables (298.15 K)" begin
        T298 = 298.15
        H(sp) = enthalpy_J_mol(sp, T298; db=db) / 1000.0  # kJ/mol

        # Pure elements in their reference state have ΔfH°(298.15 K) = 0
        for elem in [:H2, :N2, :O2, :C_s]
            @test H(elem) ≈ 0.0 atol=0.05
        end
        @test sensible_enthalpy_J_mol(:CH4, T298; db=db) == 0.0

        # Cross-check against the qualification exam reaction tables (kJ/mol)
        @test (2H(:CO) + 2H(:H2) - H(:CH4) - H(:CO2)) ≈  247.0  atol=0.5  # DRM (+247)
        @test (H(:CO) + H(:H2O) - H(:CO2) - H(:H2))   ≈   41.0  atol=0.5  # RWGS (+41)
        @test (H(:C_s) + H(:CO2) - 2H(:CO))           ≈ -172.4  atol=0.5  # Boudouard (-172.4)
        @test (H(:C_s) + 2H(:H2) - H(:CH4))           ≈   74.9  atol=0.5  # CH4 cracking (+74.9)
        @test (H(:C_s) + 2H(:H2O) - H(:CO2) - 2H(:H2))≈  -90.0  atol=0.5  # CO2 hydrogenation to C(s) (-90)
        @test (H(:C_s) + H(:H2O) - H(:CO) - H(:H2))   ≈ -131.3  atol=0.5  # CO hydrogenation to C(s) (-131.3)
        @test (H(:CH4) + 2H(:H2O) - H(:CO2) - 4H(:H2))≈ -165.0  atol=0.5  # Sabatier (-165)
        @test (H(:CH4) + H(:H2O) - H(:CO) - 3H(:H2))  ≈ -206.2  atol=0.5  # CO methanation (-206.2)

        # Extended species (if present in nasa_gas.yaml)
        if all(haskey(db, sp) for sp in EXTENDED_SPECIES)
            @test (H(:C2H6) + H(:CO) + H(:H2O) - 2H(:CH4) - H(:CO2)) ≈ 106.0 atol=1.5
            @test (H(:C2H4) + H(:H2) - H(:C2H6))                     ≈ 136.0 atol=1.5
            @test (H(:CH3OH) - H(:CO) - 2H(:H2))                     ≈ -90.6 atol=1.0
            @test (H(:CH3OCH3) + H(:H2O) - 2H(:CH3OH))               ≈ -37.0 atol=15.0 # gas vs liquid ref check
        end
    end

    @testset "Domain & lookup guards" begin
        @test_throws KeyError get_species_data(:UNKNOWN; db=db)
        @test_throws DomainError enthalpy_J_mol(:CH4, 50.0; db=db)
        @test_throws DomainError gibbs_energy_J_mol(:CO2, 10000.0; db=db)
        @test_throws ArgumentError species_list(:invalid_mode; db=db)
    end
end