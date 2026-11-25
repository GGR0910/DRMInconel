# =============================================================================
# Caracterização Geométrica de Monólitos Metálicos (Inconel 600)
# Calcula porosidade (ε), área específica (Sv), massa metálica e diâmetro hidráulico (dh)
# =============================================================================

abstract type AbstractMonolith end

"""
    RolledFoilMonolith

Monólito enrolado em espiral composto por uma folha lisa e uma folha corrugada
(configuração experimental real do projeto: 2 folhas 3x30 cm, 50 µm, D=2 cm, L=3 cm).
"""
struct RolledFoilMonolith <: AbstractMonolith
    name::String
    diameter_m::Float64
    length_m::Float64
    foil_thickness_m::Float64
    flat_strip_length_m::Float64
    corrugated_strip_length_m::Float64
    contact_factor::Float64

    function RolledFoilMonolith(;
        name::String = "Rolled Foil (Exp. 50 um)",
        diameter_m::Real = 0.02,
        length_m::Real = 0.03,
        foil_thickness_m::Real = 50.0e-6,
        flat_strip_length_m::Real = 0.30,
        corrugated_strip_length_m::Real = 0.30,
        contact_factor::Real = 0.5,
    )
        D, L, w = Float64(diameter_m), Float64(length_m), Float64(foil_thickness_m)
        Lf, Lc = Float64(flat_strip_length_m), Float64(corrugated_strip_length_m)
        cf = Float64(contact_factor)
        (D > 0 && L > 0 && w > 0 && Lf > 0 && Lc > 0) ||
            throw(ArgumentError("Dimensões do monólito devem ser positivas."))
        (0.0 < cf <= 1.0) || throw(ArgumentError("contact_factor deve estar em (0, 1]."))
        V_bulk = pi * (D / 2)^2 * L
        V_metal = (Lf + Lc) * L * w
        V_metal < V_bulk || throw(ArgumentError("Volume metálico excede o volume do cilindro."))
        return new(name, D, L, w, Lf, Lc, cf)
    end
end

"""
    HoneycombMonolith

Monólito colmeia de canais retos quadrados definido por CPSI (células/pol²) e espessura de parede.
"""
struct HoneycombMonolith <: AbstractMonolith
    name::String
    diameter_m::Float64
    length_m::Float64
    cpsi::Float64
    wall_thickness_m::Float64
    contact_factor::Float64

    function HoneycombMonolith(;
        name::String = "Honeycomb 400 CPSI",
        diameter_m::Real = 0.02,
        length_m::Real = 0.03,
        cpsi::Real = 400.0,
        wall_thickness_m::Real = 150.0e-6,
        contact_factor::Real = 0.75,
    )
        D, L = Float64(diameter_m), Float64(length_m)
        cpsi_f, w, cf = Float64(cpsi), Float64(wall_thickness_m), Float64(contact_factor)
        (D > 0 && L > 0 && cpsi_f > 0 && w > 0) ||
            throw(ArgumentError("Parâmetros geométricos devem ser positivos."))
        pitch = 0.0254 / sqrt(cpsi_f)
        w < pitch || throw(ArgumentError("Espessura de parede >= passo da célula."))
        (0.0 < cf <= 1.0) || throw(ArgumentError("contact_factor deve estar em (0, 1]."))
        return new(name, D, L, cpsi_f, w, cf)
    end
end

"""
    MetalFoamMonolith

Espuma metálica de poros abertos caracterizada por PPI (poros/polegada) e porosidade ε.
Usa relação geométrica padrão de filamentos (Calmidi-Mahajan / Richardson).
"""
struct MetalFoamMonolith <: AbstractMonolith
    name::String
    diameter_m::Float64
    length_m::Float64
    ppi::Float64
    porosity_val::Float64
    contact_factor::Float64

    function MetalFoamMonolith(;
        name::String = "Metal Foam 30 PPI",
        diameter_m::Real = 0.02,
        length_m::Real = 0.03,
        ppi::Real = 30.0,
        porosity::Real = 0.90,
        contact_factor::Real = 0.65,
    )
        D, L = Float64(diameter_m), Float64(length_m)
        ppi_f, eps, cf = Float64(ppi), Float64(porosity), Float64(contact_factor)
        (D > 0 && L > 0 && ppi_f > 0) || throw(ArgumentError("Dimensões e PPI devem ser positivos."))
        (0.5 < eps < 0.99) || throw(ArgumentError("Porosidade da espuma deve estar em (0.5, 0.99)."))
        (0.0 < cf <= 1.0) || throw(ArgumentError("contact_factor deve estar em (0, 1]."))
        return new(name, D, L, ppi_f, eps, cf)
    end
end

"""
    PrintedLatticeMonolith

Estrutura periódica cúbica impressa em 3D (SLM/LPBF) com passo `cell_size_m` e
diâmetro de haste `strut_diameter_m` (com correção exata da interseção triaxial nos nós).
"""
struct PrintedLatticeMonolith <: AbstractMonolith
    name::String
    diameter_m::Float64
    length_m::Float64
    cell_size_m::Float64
    strut_diameter_m::Float64
    contact_factor::Float64

    function PrintedLatticeMonolith(;
        name::String = "3D Printed Lattice (0.4 mm)",
        diameter_m::Real = 0.02,
        length_m::Real = 0.03,
        cell_size_m::Real = 2.0e-3,
        strut_diameter_m::Real = 0.4e-3,
        contact_factor::Real = 0.85,
    )
        D, L = Float64(diameter_m), Float64(length_m)
        a, ds, cf = Float64(cell_size_m), Float64(strut_diameter_m), Float64(contact_factor)
        (D > 0 && L > 0 && a > 0 && ds > 0) || throw(ArgumentError("Dimensões devem ser positivas."))
        ds < 0.8 * a || throw(ArgumentError("strut_diameter_m muito grande para o passo da célula."))
        (0.0 < cf <= 1.0) || throw(ArgumentError("contact_factor deve estar em (0, 1]."))
        return new(name, D, L, a, ds, cf)
    end
end

"""
    experimental_monolith(; contact_factor=0.5) -> RolledFoilMonolith

Retorna o monólito experimental exato (Goodfellow Inconel 600, 2 folhas 3x30 cm, 50 µm, 7.62 g).
"""
experimental_monolith(; contact_factor::Real=0.5) = RolledFoilMonolith(; contact_factor=contact_factor)

# --- Funções comuns de geometria ---
bulk_volume_m3(m::AbstractMonolith) = pi * (m.diameter_m / 2)^2 * m.length_m
element_thickness_m(m::RolledFoilMonolith) = m.foil_thickness_m
element_thickness_m(m::HoneycombMonolith) = m.wall_thickness_m
element_thickness_m(m::PrintedLatticeMonolith) = m.strut_diameter_m

function element_thickness_m(m::MetalFoamMonolith)
    d_p = 0.0254 / m.ppi
    solid_frac = 1.0 - m.porosity_val
    return d_p * 1.18 * sqrt(solid_frac / (3.0 * pi)) / (1.0 - exp(-solid_frac / 0.04))
end

# Porosidade (ε)
function porosity(m::RolledFoilMonolith)
    V_metal = (m.flat_strip_length_m + m.corrugated_strip_length_m) * m.length_m * m.foil_thickness_m
    return 1.0 - V_metal / bulk_volume_m3(m)
end

function porosity(m::HoneycombMonolith)
    pitch = 0.0254 / sqrt(m.cpsi)
    open_w = pitch - m.wall_thickness_m
    return (open_w / pitch)^2
end

porosity(m::MetalFoamMonolith) = m.porosity_val

function porosity(m::PrintedLatticeMonolith)
    r = m.strut_diameter_m / m.cell_size_m
    # União exata de 3 cilindros ortogonais (sólido de Steinmetz triaxial)
    solid_frac = 3.0 * (pi / 4.0) * r^2 - sqrt(2.0) * r^3
    return 1.0 - solid_frac
end

metal_volume_m3(m::AbstractMonolith) = (1.0 - porosity(m)) * bulk_volume_m3(m)

# Área superficial específica Sv [m² / m³]
function specific_surface_area(m::RolledFoilMonolith)
    # Ambas as faces das duas folhas (lisa + corrugada)
    A_geo = 2.0 * (m.flat_strip_length_m + m.corrugated_strip_length_m) * m.length_m
    return A_geo / bulk_volume_m3(m)
end

function specific_surface_area(m::HoneycombMonolith)
    pitch = 0.0254 / sqrt(m.cpsi)
    open_w = pitch - m.wall_thickness_m
    return 4.0 * open_w / (pitch^2)
end

function specific_surface_area(m::MetalFoamMonolith)
    d_s = element_thickness_m(m)
    return 4.0 * (1.0 - m.porosity_val) / d_s
end

function specific_surface_area(m::PrintedLatticeMonolith)
    a, ds = m.cell_size_m, m.strut_diameter_m
    # Derivada exata dV/dr da união de 3 cilindros ortogonais
    A_cell = 3.0 * pi * ds * a - 6.0 * sqrt(2.0) * ds^2
    return A_cell / (a^3)
end

surface_area_m2(m::AbstractMonolith) = specific_surface_area(m) * bulk_volume_m3(m)
hydraulic_diameter_m(m::AbstractMonolith) = 4.0 * porosity(m) / specific_surface_area(m)
monolith_mass_kg(m::AbstractMonolith; T_C::Real=20.0) = inconel_density(T_C) * metal_volume_m3(m)

"""
    monolith_summary(m::AbstractMonolith; T_C=20.0) -> NamedTuple

Resumo quantitativo completo da geometria do monólito.
"""
function monolith_summary(m::AbstractMonolith; T_C::Real=20.0)
    eps = porosity(m)
    Sv = specific_surface_area(m)
    dh = hydraulic_diameter_m(m)
    A_m2 = surface_area_m2(m)
    mass_g = monolith_mass_kg(m; T_C=T_C) * 1000.0
    return (
        name = m.name,
        diameter_cm = m.diameter_m * 100.0,
        length_cm = m.length_m * 100.0,
        porosity_pct = 100.0 * eps,
        solid_frac_pct = 100.0 * (1.0 - eps),
        wall_or_strut_um = element_thickness_m(m) * 1.0e6,
        surface_area_cm2 = A_m2 * 1.0e4,
        Sv_m2_m3 = Sv,
        hydraulic_diam_mm = dh * 1000.0,
        mass_g = mass_g,
        area_per_gram_cm2_g = (A_m2 * 1.0e4) / mass_g,
    )
end