# =============================================================================
# Propriedades Termofísicas e Elétricas do Inconel 600 (UNS N06600)
# Fonte: Special Metals Corporation — INCONEL alloy 600 Technical Bulletin (SMC-027)
# =============================================================================

const INCONEL600_DENSITY_20C = 8470.0           # kg / m^3 (8.47 g/cm^3 a 20 °C)
const INCONEL600_MU_R = 1.010                   # permeabilidade relativa a 200 Oe (não-ferromagnético)
const INCONEL600_CURIE_C = -124.0               # °C (Curie criogênico)
const INCONEL600_MELTING_RANGE_C = (1354.0, 1413.0) # °C (solidus, liquidus)

# Tabela 3 (Special Metals SMC-027) — Temperatura (°C)
const _INC_T_C = Float64[
    20.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0
]

# Resistividade elétrica rho(T) [Ohm * m] (tabela em micro-ohm * m)
const _INC_RHO_E_OHM_M = Float64[
    1.03, 1.04, 1.05, 1.07, 1.09, 1.12, 1.13, 1.13, 1.13, 1.15
] .* 1.0e-6

# Calor específico cp(T) [J / (kg * K)]
const _INC_CP_J_KG_K = Float64[
    444.0, 465.0, 486.0, 502.0, 519.0, 536.0, 578.0, 595.0, 611.0, 628.0
]

# Condutividade térmica k(T) [W / (m * K)] — tabelada até 800 °C na SMC-027
const _INC_K_T_C = Float64[
    20.0, 100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0
]
const _INC_K_W_M_K = Float64[
    14.9, 15.9, 17.3, 19.0, 20.5, 22.1, 23.9, 25.7, 27.5
]

# Coeficiente médio de expansão térmica linear (21 °C até T) [1 / K]
# (O ponto de 20 °C na tabela original é medido a partir de -196 °C, por isso inicia-se em 100 °C)
const _INC_CTE_T_C = Float64[
    100.0, 200.0, 300.0, 400.0, 500.0, 600.0, 700.0, 800.0, 900.0
]
const _INC_CTE_1_K = Float64[
    13.3, 13.8, 14.2, 14.5, 14.9, 15.3, 15.8, 16.1, 16.4
] .* 1.0e-6

function _check_inconel_temp(T_C::Real)
    T = Float64(T_C)
    if T < 20.0 || T > 1000.0
        throw(DomainError(T, "Temperatura fora da faixa validada para Inconel 600 [20 °C, 1000 °C]."))
    end
    return T
end

function _interp_linear(xs::Vector{Float64}, ys::Vector{Float64}, x::Float64)
    if x <= xs[1]
        slope = (ys[2] - ys[1]) / (xs[2] - xs[1])
        return ys[1] + (x - xs[1]) * slope
    elseif x >= xs[end]
        slope = (ys[end] - ys[end - 1]) / (xs[end] - xs[end - 1])
        return ys[end] + (x - xs[end]) * slope
    end
    i = searchsortedlast(xs, x)
    t = (x - xs[i]) / (xs[i + 1] - xs[i])
    return ys[i] + t * (ys[i + 1] - ys[i])
end

"""
    inconel_resistivity(T_C::Real) -> Float64

Resistividade elétrica do Inconel 600 em `Ω·m` na temperatura `T_C` (°C).
"""
function inconel_resistivity(T_C::Real)
    T = _check_inconel_temp(T_C)
    return _interp_linear(_INC_T_C, _INC_RHO_E_OHM_M, T)
end

"""
    inconel_elec_conductivity(T_C::Real) -> Float64

Condutividade elétrica `σ(T) = 1 / ρ(T)` do Inconel 600 em `S/m`.
"""
inconel_elec_conductivity(T_C::Real) = 1.0 / inconel_resistivity(T_C)

"""
    inconel_thermal_conductivity(T_C::Real) -> Float64

Condutividade térmica `k(T)` do Inconel 600 em `W/(m·K)`.
Interpolação linear até 800 °C e extrapolação linear suave até 900–1000 °C.
"""
function inconel_thermal_conductivity(T_C::Real)
    T = _check_inconel_temp(T_C)
    return _interp_linear(_INC_K_T_C, _INC_K_W_M_K, T)
end

"""
    inconel_cp(T_C::Real) -> Float64

Calor específico `cp(T)` do Inconel 600 em `J/(kg·K)`.
"""
function inconel_cp(T_C::Real)
    T = _check_inconel_temp(T_C)
    return _interp_linear(_INC_T_C, _INC_CP_J_KG_K, T)
end

"""
    inconel_cte(T_C::Real) -> Float64

Coeficiente médio de expansão térmica linear (de 20 °C até `T_C`) em `1/K`.
"""
function inconel_cte(T_C::Real)
    T = _check_inconel_temp(T_C)
    if T <= _INC_CTE_T_C[1]
        return _INC_CTE_1_K[1]
    end
    return _interp_linear(_INC_CTE_T_C, _INC_CTE_1_K, T)
end

"""
    inconel_density(T_C::Real=20.0) -> Float64

Densidade `ρ_m(T)` do Inconel 600 em `kg/m³`, corrigida pela expansão volumétrica isotrópica.
Retorna exatamente `8470.0 kg/m³` a 20 °C.
"""
function inconel_density(T_C::Real=20.0)
    T = _check_inconel_temp(T_C)
    dT = T - 20.0
    linear_factor = 1.0 + inconel_cte(T) * dT
    return INCONEL600_DENSITY_20C / (linear_factor^3)
end

"""
    inconel_permeability_rel(T_C::Real=20.0) -> Float64

Permeabilidade magnética relativa `μ_r` do Inconel 600.
Acima da temperatura de Curie (-124 °C), a liga é paramagnética (μ_r ≈ 1.010).
"""
function inconel_permeability_rel(T_C::Real=20.0)
    _check_inconel_temp(T_C)
    return INCONEL600_MU_R
end

"""
    inconel_thermal_diffusivity(T_C::Real) -> Float64

Difusividade térmica `α(T) = k(T) / (ρ_m(T) * cp(T))` em `m²/s`.
"""
function inconel_thermal_diffusivity(T_C::Real)
    k = inconel_thermal_conductivity(T_C)
    rho_m = inconel_density(T_C)
    cp = inconel_cp(T_C)
    return k / (rho_m * cp)
end

"""
    inconel_properties(T_C::Real) -> NamedTuple

Retorna todas as propriedades termofísicas e eletromagnéticas do Inconel 600 em `T_C` (°C).
"""
function inconel_properties(T_C::Real)
    T = _check_inconel_temp(T_C)
    rho_e = inconel_resistivity(T)
    k = inconel_thermal_conductivity(T)
    cp = inconel_cp(T)
    rho_m = inconel_density(T)
    mu_r = inconel_permeability_rel(T)
    return (
        T_C = T,
        rho_e_ohm_m = rho_e,
        sigma_S_m = 1.0 / rho_e,
        k_W_m_K = k,
        cp_J_kg_K = cp,
        density_kg_m3 = rho_m,
        mu_r = mu_r,
        alpha_m2_s = k / (rho_m * cp),
    )
end