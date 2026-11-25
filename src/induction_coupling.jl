# =============================================================================
# Acoplamento Eletromagnético Indutivo e Efeito Pelicular (Skin Depth)
# Referências analíticas: Davies (1990) "Conduction and Induction Heating";
# Rudnev et al. (2017) "Handbook of Induction Heating".
# =============================================================================

const MU_0 = 4.0 * pi * 1.0e-7          # Permeabilidade do vácuo [H / m]
const RHO_COPPER_COIL = 1.85e-8         # Resistividade do cobre refrigerado (~40 °C) [Ω·m]

"""
    skin_depth_m(f_hz::Real, T_C::Real) -> Float64

Profundidade de penetração eletromagnética (skin depth) `δ(f, T)` no Inconel 600 em metros:
    δ = sqrt(ρ(T) / (π * f * μ₀ * μᵣ))
"""
function skin_depth_m(f_hz::Real, T_C::Real)
    f = Float64(f_hz)
    f > 0.0 || throw(DomainError(f, "Frequência f_hz deve ser estritamente positiva."))
    rho_e = inconel_resistivity(T_C)
    mu_r = inconel_permeability_rel(T_C)
    return sqrt(rho_e / (pi * f * MU_0 * mu_r))
end

skin_depth_mm(f_hz::Real, T_C::Real) = 1000.0 * skin_depth_m(f_hz, T_C)

"""
    slab_power_factor(xi::Real) -> Float64

Fator analítico adimensional de dissipação indutiva em placa plana de espessura `w`,
com `ξ = w / δ`:
    F_slab(ξ) = (sinh(ξ) - sin(ξ)) / (cosh(ξ) + cos(ξ))

Limites físicos:
- Placa fina (ξ << 1): F_slab ≈ ξ³ / 6  (peça transparente ao campo)
- Placa espessa (ξ >> 1): F_slab -> 1.0 (limite de efeito pelicular superficial)
"""
function slab_power_factor(xi::Real)
    x = Float64(xi)
    x >= 0.0 || throw(DomainError(x, "Razão w/δ deve ser não-negativa."))
    if x < 1.0e-3
        # Série de Taylor exata: evita cancelamento numérico em sinh(x) - sin(x)
        return (x^3 / 6.0) * (1.0 - (17.0 / 420.0) * x^4)
    elseif x > 40.0
        return 1.0
    end
    return (sinh(x) - sin(x)) / (cosh(x) + cos(x))
end

# Razão I₁(z) / I₀(z) para z = (1 + i)*u (u > 0) em Julia puro com precisão de máquina
function _bessel_I1_over_I0(z::ComplexF64)
    az = abs(z)
    if az > 35.0
        # Expansão assintótica de Hankel: num = I1/factor, den = I0/factor
        num = 1.0 + 0.0im
        den = 1.0 + 0.0im
        term_num = 1.0 + 0.0im
        term_den = 1.0 + 0.0im
        inv_z = 1.0 / z
        for k in 1:14
            odd_sq = (2.0 * k - 1.0)^2
            term_num *= -(4.0 - odd_sq) * inv_z / (8.0 * k)
            term_den *= -(-odd_sq) * inv_z / (8.0 * k)
            num += term_num
            den += term_den
        end
        return num / den
    else
        # Fração contínua exata de Gauss: I₁(z)/I₀(z) = 1 / (2/z + 1 / (4/z + ...))
        M = Int(ceil(2.5 * az)) + 45
        r = 0.0 + 0.0im
        inv_z = 1.0 / z
        for k in M:-1:1
            r = 1.0 / ((2.0 * k) * inv_z + r)
        end
        return r
    end
end

"""
    cylinder_power_factor(zeta::Real) -> Float64

Fator analítico adimensional de dissipação indutiva em cilindro condutor de diâmetro `D`
sob campo axial alternado, com `ζ = D / δ`:
    F_cyl(ζ) = Re[(1 + i) * I₁((1 + i)ζ/2) / I₀((1 + i)ζ/2)]

Limites físicos:
- Cilindro fino (ζ << 1): F_cyl ≈ ζ³ / 32 (transparente ao campo)
- Cilindro espesso (ζ >> 1): F_cyl ≈ 1 - 1/ζ -> 1.0
"""
function cylinder_power_factor(zeta::Real)
    z_ratio = Float64(zeta)
    z_ratio >= 0.0 || throw(DomainError(z_ratio, "Razão D/δ deve ser não-negativa."))
    if z_ratio < 1.0e-3
        return (z_ratio^3) / 32.0
    end
    u = 0.5 * z_ratio
    z = complex(u, u)
    return real((1.0 + 1.0im) * _bessel_I1_over_I0(z))
end

"""
    coil_electrical_efficiency(power_factor, rho_workpiece; mu_r=INCONEL600_MU_R, D_coil_ratio=1.5, rho_coil=RHO_COPPER_COIL) -> Float64

Eficiência elétrica analítica (solenoide longo) de transferência de potência da bobina
de cobre para a peça metálica:
    η_el = 1 / (1 + (D_coil / D_w) * sqrt(ρ_coil / (μ_r * ρ_w)) / F)
"""
function coil_electrical_efficiency(
    power_factor::Real,
    rho_workpiece::Real;
    mu_r::Real = INCONEL600_MU_R,
    D_coil_ratio::Real = 1.5,
    rho_coil::Real = RHO_COPPER_COIL,
)
    F = Float64(power_factor)
    F <= 0.0 && return 0.0
    loss_ratio = Float64(D_coil_ratio) * sqrt(Float64(rho_coil) / (Float64(mu_r) * Float64(rho_workpiece))) / F
    return 1.0 / (1.0 + loss_ratio)
end

"""
    effective_resistivity(m::AbstractMonolith, T_C::Real; contact_factor=m.contact_factor) -> Float64

Resistividade elétrica efetiva macroscópica do monólito poroso [Ω·m]:
    ρ_eff = ρ_metal(T) / (contact_factor * (1 - ε))
"""
function effective_resistivity(m::AbstractMonolith, T_C::Real; contact_factor::Real=m.contact_factor)
    cf = Float64(contact_factor)
    (0.0 < cf <= 1.0) || throw(ArgumentError("contact_factor deve estar em (0, 1]."))
    solid_frac = 1.0 - porosity(m)
    return inconel_resistivity(T_C) / (cf * solid_frac)
end

"""
    effective_skin_depth_m(m::AbstractMonolith, f_hz::Real, T_C::Real; contact_factor=m.contact_factor) -> Float64

Skin depth efetivo macroscópico `δ_eff` do monólito poroso [m].
"""
function effective_skin_depth_m(m::AbstractMonolith, f_hz::Real, T_C::Real; contact_factor::Real=m.contact_factor)
    f = Float64(f_hz)
    f > 0.0 || throw(DomainError(f, "Frequência deve ser positiva."))
    rho_eff = effective_resistivity(m, T_C; contact_factor=contact_factor)
    mu_r = inconel_permeability_rel(T_C)
    return sqrt(rho_eff / (pi * f * MU_0 * mu_r))
end

"""
    coupling_regime(ratio::Real; geometry::Symbol=:cylinder) -> Symbol

Classifica o regime eletromagnético de operação:
- `:transparent` (acoplamento fraco, peça fina demais para a frequência)
- `:optimal` (boa penetração e alta eficiência de acoplamento)
- `:surface_concentrated` (correntes concentradas na casca externa, risco de gradiente radial)
"""
function coupling_regime(ratio::Real; geometry::Symbol=:cylinder)
    r = Float64(ratio)
    if geometry === :slab
        r < 1.5 && return :transparent
        r <= 6.0 && return :optimal
        return :surface_concentrated
    elseif geometry === :cylinder
        r < 3.5 && return :transparent
        r <= 12.0 && return :optimal
        return :surface_concentrated
    else
        throw(ArgumentError("geometry deve ser :cylinder ou :slab."))
    end
end

"""
    analyze_monolith_coupling(m::AbstractMonolith, f_hz::Real, T_C::Real; contact_factor=m.contact_factor, D_coil_ratio=1.5) -> NamedTuple

Avalia simultaneamente o limite local (parede/folha isolada) e o acoplamento macroscópico
do cilindro poroso para o monólito `m` na frequência `f_hz` e temperatura `T_C`.
"""
function analyze_monolith_coupling(
    m::AbstractMonolith,
    f_hz::Real,
    T_C::Real;
    contact_factor::Real = m.contact_factor,
    D_coil_ratio::Real = 1.5,
)
    rho_m = inconel_resistivity(T_C)
    mu_r = inconel_permeability_rel(T_C)
    delta_metal = skin_depth_m(f_hz, T_C)

    # 1. Limite de folha/parede isolada (sem laços fechados de corrente entre voltas)
    w = element_thickness_m(m)
    xi_wall = w / delta_metal
    F_wall = slab_power_factor(xi_wall)
    eta_wall = coil_electrical_efficiency(F_wall, rho_m; mu_r=mu_r, D_coil_ratio=D_coil_ratio)

    # 2. Acoplamento macroscópico do cilindro (correntes fechadas pela circunferência)
    rho_eff = effective_resistivity(m, T_C; contact_factor=contact_factor)
    delta_eff = effective_skin_depth_m(m, f_hz, T_C; contact_factor=contact_factor)
    zeta_macro = m.diameter_m / delta_eff
    F_macro = cylinder_power_factor(zeta_macro)
    eta_macro = coil_electrical_efficiency(F_macro, rho_eff; mu_r=mu_r, D_coil_ratio=D_coil_ratio)

    return (
        name = m.name,
        f_kHz = Float64(f_hz) / 1000.0,
        T_C = Float64(T_C),
        delta_metal_mm = delta_metal * 1000.0,
        delta_eff_mm = delta_eff * 1000.0,
        wall_to_delta = xi_wall,
        F_wall_isolated = F_wall,
        eta_wall_pct = 100.0 * eta_wall,
        D_to_delta_eff = zeta_macro,
        F_macro_cyl = F_macro,
        eta_macro_pct = 100.0 * eta_macro,
        macro_regime = coupling_regime(zeta_macro; geometry=:cylinder),
    )
end