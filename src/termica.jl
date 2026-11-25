# =============================================================================
# Fase 3 — Dinâmica de Aquecimento (0D) e Perfil Térmico Radial (1D) no Quartzo
# Stack: OrdinaryDiffEq.jl (SciML / DifferentialEquations) + Módulos Fases 1 e 2
# =============================================================================

using OrdinaryDiffEq

const SIGMA_SB = 5.670374419e-8   # Constante de Stefan-Boltzmann [W / (m^2 * K^4)]
const GRAVITY_G = 9.80665         # Aceleração da gravidade [m / s^2]

# Emissividades típicas do Inconel 600 conforme estado superficial
const INCONEL_EMISSIVITY_PRESETS = Dict{Symbol, Float64}(
    :clean    => 0.30,  # Metal polido / não tratado
    :oxidized => 0.75,  # Ativado em ar a 800 °C (condição experimental padrão)
    :coked    => 0.90,  # Superfície escurecida por deposição de carbono
)

# Polinômios Cp/R = A + B*T + C*T^2 + D*T^-2 (Smith-Van Ness-Abbott / NIST, 298–1300 K)
# Usados para separar Q_gas (calor sensível) de Q_reac (entalpia da reação)
const _CP_COEFFS = Dict{Symbol, NTuple{4, Float64}}(
    :CH4 => (1.702, 9.081e-3, -2.164e-6,  0.0),
    :CO2 => (5.457, 1.045e-3,  0.0,      -1.157e5),
    :N2  => (3.280, 0.593e-3,  0.0,       0.040e5),
    :O2  => (3.639, 0.506e-3,  0.0,      -0.227e5),
    :H2  => (3.249, 0.422e-3,  0.0,       0.083e5),
    :CO  => (3.376, 0.557e-3,  0.0,      -0.031e5),
    :H2O => (3.470, 1.450e-3,  0.0,       0.121e5),
)

"""
    ThermalReactorConfig

Parâmetros térmicos do monólito de Inconel 600 dentro do tubo de quartzo.
"""
struct ThermalReactorConfig
    quartz_id_m::Float64
    quartz_od_m::Float64
    quartz_hot_length_m::Float64
    emissivity_monolith::Float64
    emissivity_quartz::Float64
    quartz_ir_transmittance::Float64
    end_view_factor::Float64
    gap_gas_k_W_m_K::Float64
    ambient_T_C::Float64
    inlet_gas_T_C::Float64
    insulation_thickness_m::Float64
    insulation_k_W_m_K::Float64
    insulation_emissivity::Float64
    radial_k_contact_factor::Float64

    function ThermalReactorConfig(;
        quartz_id_m::Real = 0.022,
        quartz_od_m::Real = 0.025,
        quartz_hot_length_m::Real = 0.050,
        surface_state::Symbol = :oxidized,
        emissivity_monolith::Union{Nothing, Real} = nothing,
        emissivity_quartz::Real = 0.85,
        quartz_ir_transmittance::Real = 0.35,
        end_view_factor::Real = 0.65,
        gap_gas_k_W_m_K::Real = 0.075,
        ambient_T_C::Real = 25.0,
        inlet_gas_T_C::Real = 25.0,
        insulation_thickness_m::Real = 0.0,
        insulation_k_W_m_K::Real = 0.12,
        insulation_emissivity::Real = 0.85,
        radial_k_contact_factor::Real = 0.35,
    )
        eps_m = if isnothing(emissivity_monolith)
            haskey(INCONEL_EMISSIVITY_PRESETS, surface_state) ||
                throw(ArgumentError("surface_state inválido: use :clean, :oxidized ou :coked."))
            INCONEL_EMISSIVITY_PRESETS[surface_state]
        else
            Float64(emissivity_monolith)
        end

        qid, qod, qL = Float64(quartz_id_m), Float64(quartz_od_m), Float64(quartz_hot_length_m)
        qod > qid > 0.0 || throw(ArgumentError("Exige-se quartz_od_m > quartz_id_m > 0."))
        qL > 0.0 || throw(ArgumentError("quartz_hot_length_m deve ser positivo."))
        (0.0 < eps_m <= 1.0 && 0.0 < emissivity_quartz <= 1.0) ||
            throw(ArgumentError("Emissividades devem estar em (0, 1]."))
        (0.0 <= quartz_ir_transmittance < 1.0) ||
            throw(ArgumentError("quartz_ir_transmittance deve estar em [0, 1)."))
        (0.0 <= end_view_factor <= 1.0) ||
            throw(ArgumentError("end_view_factor deve estar em [0, 1]."))
        ins_t = Float64(insulation_thickness_m)
        ins_t >= 0.0 || throw(ArgumentError("insulation_thickness_m não pode ser negativo."))
        (0.0 < radial_k_contact_factor <= 1.0) ||
            throw(ArgumentError("radial_k_contact_factor deve estar em (0, 1]."))

        return new(
            qid, qod, qL, eps_m, Float64(emissivity_quartz),
            Float64(quartz_ir_transmittance), Float64(end_view_factor),
            Float64(gap_gas_k_W_m_K), Float64(ambient_T_C), Float64(inlet_gas_T_C),
            ins_t, Float64(insulation_k_W_m_K), Float64(insulation_emissivity),
            Float64(radial_k_contact_factor),
        )
    end
end

"""
    air_properties(T_K::Real) -> NamedTuple

Propriedades termofísicas do ar atmosférico (1 atm) em `T_K` [K].
"""
function air_properties(T_K::Real)
    T = Float64(T_K)
    T > 150.0 || throw(DomainError(T, "Temperatura do ar deve ser > 150 K."))
    rho = 101325.0 / (287.05 * T)
    mu = 1.458e-6 * (T^1.5) / (T + 110.4)
    nu = mu / rho
    k = 0.02624 * (T / 300.0)^0.82
    Pr = 0.70
    return (rho=rho, mu=mu, nu=nu, k=k, Pr=Pr, alpha=nu / Pr, beta=1.0 / T)
end

"""
    natural_convection_h(T_surf_K::Real, T_amb_K::Real, D_outer_m::Real) -> Float64

Coeficiente de convecção natural `h` [W/(m²·K)] em cilindro externo (Churchill–Chu).
"""
function natural_convection_h(T_surf_K::Real, T_amb_K::Real, D_outer_m::Real)
    dT = abs(Float64(T_surf_K) - Float64(T_amb_K))
    dT < 1e-9 && return 0.0
    T_film = 0.5 * (Float64(T_surf_K) + Float64(T_amb_K))
    air = air_properties(T_film)
    D = Float64(D_outer_m)
    Ra = (GRAVITY_G * air.beta * dT * D^3) / (air.nu * air.alpha)
    denom = (1.0 + (0.559 / air.Pr)^(9.0 / 16.0))^(8.0 / 27.0)
    Nu = (0.60 + 0.387 * (Ra^(1.0 / 6.0)) / denom)^2
    return Nu * air.k / D
end

# Integral analítica de Cp(T) de T_in até T_out [J/mol]
function _sensible_delta_h_J_mol(sp::Symbol, T_in_K::Float64, T_out_K::Float64)
    A, B, C, D = get(_CP_COEFFS, sp, _CP_COEFFS[:N2])
    term = A * (T_out_K - T_in_K) +
           0.5 * B * (T_out_K^2 - T_in_K^2) +
           (C / 3.0) * (T_out_K^3 - T_in_K^3) -
           D * (1.0 / T_out_K - 1.0 / T_in_K)
    return R_GAS * term
end

"""
    ChemistryThermalLoad

Curvas interpoladas de `Q_gas(T)` (aquecimento sensível do gás) e `Q_reac(T)`
(calor endotérmico da reação química) pré-computadas a partir do equilíbrio da Fase 1.
"""
struct ChemistryThermalLoad
    T_grid_C::Vector{Float64}
    Q_gas_W::Vector{Float64}
    Q_reac_W::Vector{Float64}
end

"""
    zero_chemistry_load() -> ChemistryThermalLoad

Carga química nula (aquecimento em atmosfera inerte / sem fluxo reativo).
"""
zero_chemistry_load() = ChemistryThermalLoad([20.0, 1000.0], [0.0, 0.0], [0.0, 0.0])

"""
    build_chemistry_load(feed; T_in_C=25.0, allow_coke=false, db=default_thermo_db()) -> ChemistryThermalLoad

Pré-calcula `Q_gas(T)` e `Q_reac(T)` de 25 °C até 950 °C acoplando o solver da Fase 1,
garantindo `Q_gas(T) + Q_reac(T) == Q_total_W(T)` em toda a faixa reativa.
"""
function build_chemistry_load(
    feed;
    T_in_C::Real = 25.0,
    allow_coke::Bool = false,
    db = default_thermo_db(),
)
    # Obtém vazões molares de entrada via uma chamada de referência
    ref = solve_equilibrium(feed; T_C=600.0, allow_coke=allow_coke, db=db)
    n_in = ref.n_in_mol_s
    T_in_K = c_to_k(T_in_C)

    T_nodes = vcat([20.0, 25.0, 100.0, 200.0], collect(300.0:25.0:950.0))
    Q_gas_vec = zeros(Float64, length(T_nodes))
    Q_reac_vec = zeros(Float64, length(T_nodes))

    for (idx, T_C) in enumerate(T_nodes)
        T_K = c_to_k(T_C)
        q_sens = sum(F * _sensible_delta_h_J_mol(sp, T_in_K, T_K) for (sp, F) in n_in)
        Q_gas_vec[idx] = max(0.0, q_sens)

        if T_C >= 300.0
            eq = solve_equilibrium(feed; T_C=T_C, allow_coke=allow_coke, db=db)
            Q_reac_vec[idx] = max(0.0, eq.Q_total_W - Q_gas_vec[idx])
        else
            Q_reac_vec[idx] = 0.0
        end
    end
    return ChemistryThermalLoad(T_nodes, Q_gas_vec, Q_reac_vec)
end

function q_gas_W(load::ChemistryThermalLoad, T_C::Real)
    T = clamp(Float64(T_C), load.T_grid_C[1], load.T_grid_C[end])
    return _interp_linear(load.T_grid_C, load.Q_gas_W, T)
end

function q_reac_W(load::ChemistryThermalLoad, T_C::Real)
    T = clamp(Float64(T_C), load.T_grid_C[1], load.T_grid_C[end])
    return _interp_linear(load.T_grid_C, load.Q_reac_W, T)
end

# -----------------------------------------------------------------------------
# Perdas Térmicas pelo Tubo de Quartzo (Radiação + Convecção Natural)
# -----------------------------------------------------------------------------
function _quartz_outer_loss(cfg::ThermalReactorConfig, T_q_K::Float64)
    T_amb_K = c_to_k(cfg.ambient_T_C)
    if T_q_K <= T_amb_K + 1e-12
        return (Q_rad=0.0, Q_conv=0.0, Q_total=0.0, T_outer_K=T_q_K, h_conv=0.0)
    end

    if cfg.insulation_thickness_m <= 0.0
        D_out = cfg.quartz_od_m
        A_out = pi * D_out * cfg.quartz_hot_length_m
        h_conv = natural_convection_h(T_q_K, T_amb_K, D_out)
        Q_rad = cfg.emissivity_quartz * SIGMA_SB * A_out * (T_q_K^4 - T_amb_K^4)
        Q_conv = h_conv * A_out * (T_q_K - T_amb_K)
        return (Q_rad=Q_rad, Q_conv=Q_conv, Q_total=Q_rad + Q_conv, T_outer_K=T_q_K, h_conv=h_conv)
    else
        D_out = cfg.quartz_od_m + 2.0 * cfg.insulation_thickness_m
        A_out = pi * D_out * cfg.quartz_hot_length_m
        R_ins = log(D_out / cfg.quartz_od_m) / (2.0 * pi * cfg.insulation_k_W_m_K * cfg.quartz_hot_length_m)

        lo, hi = T_amb_K, T_q_K
        for _ in 1:60
            mid = 0.5 * (lo + hi)
            q_cond = (T_q_K - mid) / R_ins
            h_mid = natural_convection_h(mid, T_amb_K, D_out)
            q_surf = cfg.insulation_emissivity * SIGMA_SB * A_out * (mid^4 - T_amb_K^4) +
                     h_mid * A_out * (mid - T_amb_K)
            if q_cond > q_surf
                lo = mid
            else
                hi = mid
            end
        end
        T_s = 0.5 * (lo + hi)
        h_conv = natural_convection_h(T_s, T_amb_K, D_out)
        Q_rad = cfg.insulation_emissivity * SIGMA_SB * A_out * (T_s^4 - T_amb_K^4)
        Q_conv = h_conv * A_out * (T_s - T_amb_K)
        return (Q_rad=Q_rad, Q_conv=Q_conv, Q_total=Q_rad + Q_conv, T_outer_K=T_s, h_conv=h_conv)
    end
end

"""
    surface_heat_losses(m::AbstractMonolith, cfg::ThermalReactorConfig, T_surf_C::Real; include_ends=true) -> NamedTuple

Calcula `Q_rad_W` (radiação direta + faces + reemissão do quartzo/manta) e
`Q_conv_W` (convecção natural externa) para uma superfície de monólito em `T_surf_C` (°C).
"""
function surface_heat_losses(
    m::AbstractMonolith,
    cfg::ThermalReactorConfig,
    T_surf_C::Real;
    include_ends::Bool = true,
)
    T_m_K = c_to_k(T_surf_C)
    T_amb_K = c_to_k(cfg.ambient_T_C)
    if T_m_K <= T_amb_K + 1e-12
        return (
            T_q_C = cfg.ambient_T_C, T_outer_C = cfg.ambient_T_C,
            Q_rad_direct_W = 0.0, Q_rad_ends_W = 0.0,
            Q_rad_outer_W = 0.0, Q_rad_W = 0.0, Q_conv_W = 0.0,
            Q_lat_total_W = 0.0, Q_loss_total_W = 0.0,
        )
    end

    A_lat = pi * m.diameter_m * m.length_m
    A_ends = 2.0 * pi * (m.diameter_m / 2.0)^2
    A_q_in = pi * cfg.quartz_id_m * cfg.quartz_hot_length_m

    tau_eff = cfg.insulation_thickness_m > 0.0 ? 0.0 : cfg.quartz_ir_transmittance
    rad_res = (1.0 / cfg.emissivity_monolith) +
              ((1.0 - cfg.emissivity_quartz) / cfg.emissivity_quartz) * (A_lat / A_q_in)
    ratio_gap = max(cfg.quartz_id_m / m.diameter_m, 1.01)
    G_gap = (2.0 * pi * cfg.gap_gas_k_W_m_K * m.length_m) / log(ratio_gap)

    # Temperatura de equilíbrio quasi-estacionário da parede do quartzo T_q ∈ [T_amb, T_m]
    lo, hi = T_amb_K, T_m_K
    for _ in 1:65
        mid = 0.5 * (lo + hi)
        q_in = (1.0 - tau_eff) * SIGMA_SB * A_lat * (T_m_K^4 - mid^4) / rad_res +
               G_gap * (T_m_K - mid)
        q_out = _quartz_outer_loss(cfg, mid).Q_total
        if q_in > q_out
            lo = mid
        else
            hi = mid
        end
    end
    T_q_K = 0.5 * (lo + hi)
    ext = _quartz_outer_loss(cfg, T_q_K)

    Q_rad_direct = tau_eff * cfg.emissivity_monolith * SIGMA_SB * A_lat * (T_m_K^4 - T_amb_K^4)
    Q_rad_ends = include_ends ?
        (cfg.end_view_factor * cfg.emissivity_monolith * SIGMA_SB * A_ends * (T_m_K^4 - T_amb_K^4)) : 0.0

    Q_rad_total = Q_rad_direct + Q_rad_ends + ext.Q_rad
    Q_conv_total = ext.Q_conv
    Q_lat_total = Q_rad_direct + ext.Q_total
    Q_loss_total = Q_rad_total + Q_conv_total

    return (
        T_q_C = k_to_c(T_q_K),
        T_outer_C = k_to_c(ext.T_outer_K),
        Q_rad_direct_W = Q_rad_direct,
        Q_rad_ends_W = Q_rad_ends,
        Q_rad_outer_W = ext.Q_rad,
        Q_rad_W = Q_rad_total,
        Q_conv_W = Q_conv_total,
        Q_lat_total_W = Q_lat_total,
        Q_loss_total_W = Q_loss_total,
    )
end

# -----------------------------------------------------------------------------
# 1. Modelo Transiente Concentrado (0D) via OrdinaryDiffEq.jl
#    m * cp(T) * dT/dt = P_ind - Q_rad(T) - Q_conv(T) - Q_gas(T) - Q_reac(T)
# -----------------------------------------------------------------------------

"""
    steady_state_power_0d(m, cfg, T_C; chem=zero_chemistry_load()) -> NamedTuple

Decomposição exata das parcelas de potência no estado estacionário 0D em `T_C` (°C).
"""
function steady_state_power_0d(
    m::AbstractMonolith,
    cfg::ThermalReactorConfig,
    T_C::Real;
    chem::ChemistryThermalLoad = zero_chemistry_load(),
)
    losses = surface_heat_losses(m, cfg, T_C; include_ends=true)
    qg = q_gas_W(chem, T_C)
    qr = q_reac_W(chem, T_C)
    P_ind = losses.Q_rad_W + losses.Q_conv_W + qg + qr
    return (
        T_C = Float64(T_C),
        T_q_C = losses.T_q_C,
        T_outer_C = losses.T_outer_C,
        P_ind_W = P_ind,
        Q_rad_W = losses.Q_rad_W,
        Q_conv_W = losses.Q_conv_W,
        Q_gas_W = qg,
        Q_reac_W = qr,
        pct_rad  = P_ind > 0 ? 100.0 * losses.Q_rad_W  / P_ind : 0.0,
        pct_conv = P_ind > 0 ? 100.0 * losses.Q_conv_W / P_ind : 0.0,
        pct_gas  = P_ind > 0 ? 100.0 * qg              / P_ind : 0.0,
        pct_reac = P_ind > 0 ? 100.0 * qr              / P_ind : 0.0,
    )
end

"""
    simulate_ramp_0d(m, cfg; P_ind_W, t_end_s=200.0, T0_C=cfg.ambient_T_C, chem=zero_chemistry_load()) -> NamedTuple

Resolve a ODE concentrada 0D com `OrdinaryDiffEq.Tsit5()` desde `T0_C` até `t_end_s`.
Retorna série temporal, temperatura final, tempo para atingir 600/700/800 °C e resíduo do balanço.
"""
function simulate_ramp_0d(
    m::AbstractMonolith,
    cfg::ThermalReactorConfig;
    P_ind_W::Real,
    t_end_s::Real = 200.0,
    T0_C::Real = cfg.ambient_T_C,
    chem::ChemistryThermalLoad = zero_chemistry_load(),
)
    P_in = Float64(P_ind_W)
    P_in >= 0.0 || throw(ArgumentError("P_ind_W deve ser não-negativa."))
    mass_kg = monolith_mass_kg(m; T_C=20.0)

    function rhs_0d(T_C, p, t)
        T_clamped = clamp(T_C, 20.0, 1000.0)
        C_th = mass_kg * inconel_cp(T_clamped)
        losses = surface_heat_losses(m, cfg, T_C; include_ends=true)
        qg = q_gas_W(chem, T_C)
        qr = q_reac_W(chem, T_C)
        return (P_in - losses.Q_rad_W - losses.Q_conv_W - qg - qr) / C_th
    end

    prob = ODEProblem(rhs_0d, Float64(T0_C), (0.0, Float64(t_end_s)))
    sol = solve(prob, Tsit5(); reltol=1e-9, abstol=1e-9, maxiters=10^6)

    t_grid = collect(range(0.0, Float64(t_end_s), length=400))
    T_grid = [sol(t) for t in t_grid]
    T_final = sol.u[end]

    bd = steady_state_power_0d(m, cfg, T_final; chem=chem)
    balance_rel_err = P_in > 0.0 ? abs(P_in - bd.P_ind_W) / P_in : abs(bd.P_ind_W)

    function time_to_temp(target::Float64)
        idx = findfirst(x -> x >= target, T_grid)
        return isnothing(idx) ? NaN : t_grid[idx]
    end

    return (
        t_s = t_grid,
        T_C = T_grid,
        T_final_C = T_final,
        P_ind_W = P_in,
        t_600_s = time_to_temp(600.0),
        t_700_s = time_to_temp(700.0),
        t_800_s = time_to_temp(800.0),
        breakdown_final = bd,
        balance_rel_err = balance_rel_err,
    )
end

# -----------------------------------------------------------------------------
# 2. Modelo Radial 1D em Coordenadas Cilíndricas (Método das Linhas + OrdinaryDiffEq)
#    ρ_eff * cp(T) * ∂T/∂t = (1/r) ∂/∂r (r * k_eff(T) * ∂T/∂r) + q_ind(r) - q_chem(T)
# -----------------------------------------------------------------------------

# Perfil adimensional |I₁((1+i)*r/δ)|² / |I₁((1+i)*R/δ)|² (correntes parasitas exatas em cilindro)
function _radial_eddy_heating_weight(r_m::Float64, R_m::Float64, delta_eff_m::Float64)
    r_m <= 0.0 && return 0.0
    u_r = r_m / delta_eff_m
    u_R = R_m / delta_eff_m
    if u_R > 14.0
        # Forma assintótica exata de Hankel (evita overflow em I₁)
        return (R_m / r_m) * exp(-2.0 * (u_R - u_r))
    end
    function bessel_I1_series(z::ComplexF64)
        half_z = 0.5 * z
        term = half_z
        acc = term
        z2_4 = half_z * half_z
        for k in 1:30
            term *= z2_4 / (k * (k + 1))
            acc += term
            abs(term) < 1e-16 * abs(acc) && break
        end
        return acc
    end
    I_r = bessel_I1_series(complex(u_r, u_r))
    I_R = bessel_I1_series(complex(u_R, u_R))
    return abs2(I_r) / abs2(I_R)
end

"""
    effective_radial_k(m::AbstractMonolith, cfg::ThermalReactorConfig, T_C::Real) -> Float64

Condutividade térmica radial efetiva `k_eff(T)` [W/(m·K)] do monólito poroso.
"""
function effective_radial_k(m::AbstractMonolith, cfg::ThermalReactorConfig, T_C::Real)
    T_clamped = clamp(Float64(T_C), 20.0, 1000.0)
    eps = porosity(m)
    k_metal = inconel_thermal_conductivity(T_clamped)
    return cfg.radial_k_contact_factor * (1.0 - eps) * k_metal + eps * cfg.gap_gas_k_W_m_K
end

"""
    simulate_radial_1d(m, cfg; P_ind_W, f_hz=200e3, N=25, t_end_s=250.0, T0_C=cfg.ambient_T_C, chem=zero_chemistry_load(), include_ends=true) -> NamedTuple

Resolve a EDP cilíndrica radial 1D por Volumes Finitos (Método das Linhas) via `OrdinaryDiffEq.Tsit5()`.
Impõe simetria exata `∂T/∂r = 0` em `r = 0` e balanço de radiação/convecção na parede `r = R`.
"""
function simulate_radial_1d(
    m::AbstractMonolith,
    cfg::ThermalReactorConfig;
    P_ind_W::Real,
    f_hz::Real = 200.0e3,
    N::Int = 25,
    t_end_s::Real = 250.0,
    T0_C::Real = cfg.ambient_T_C,
    chem::ChemistryThermalLoad = zero_chemistry_load(),
    include_ends::Bool = true,
)
    N >= 5 || throw(ArgumentError("Número de células radiais N deve ser >= 5."))
    P_in = Float64(P_ind_W)
    R = 0.5 * m.diameter_m
    L = m.length_m
    dr = R / N

    # Malha de volumes finitos concêntricos
    r_faces = [i * dr for i in 0:N]
    r_centers = [(i - 0.5) * dr for i in 1:N]
    A_ends_ring = [2.0 * pi * (r_faces[i + 1]^2 - r_faces[i]^2) for i in 1:N]
    V_cell = [pi * (r_faces[i + 1]^2 - r_faces[i]^2) * L for i in 1:N]
    V_total = sum(V_cell)
    m_cell = [monolith_mass_kg(m; T_C=20.0) * (V_cell[i] / V_total) for i in 1:N]

    # Perfil radial de potência induzida normalizado para somar exatamente P_ind_W
    delta_eff = effective_skin_depth_m(m, f_hz, 800.0)
    raw_weights = [_radial_eddy_heating_weight(r_centers[i], R, delta_eff) * V_cell[i] for i in 1:N]
    w_sum = sum(raw_weights)
    P_cell = [P_in * (raw_weights[i] / w_sum) for i in 1:N]
    q_ind_vol = [P_cell[i] / V_cell[i] for i in 1:N]

    T_amb_K = c_to_k(cfg.ambient_T_C)

    function rhs_1d!(dT, T, p, t)
        # Fluxos condutivos nas faces internas (i = 1..N-1)
        # Na face i = 0 (r = 0), simetria exata => Q_inner = 0.0
        Q_prev = 0.0
        for i in 1:N
            if i < N
                k_face = 0.5 * (effective_radial_k(m, cfg, T[i]) + effective_radial_k(m, cfg, T[i + 1]))
                A_face = 2.0 * pi * r_faces[i + 1] * L
                Q_next = -k_face * A_face * (T[i + 1] - T[i]) / dr
            else
                # Fronteira externa r = R: perda lateral total (radiação + convecção pelo quartzo/manta)
                Q_next = surface_heat_losses(m, cfg, T[N]; include_ends=false).Q_lat_total_W
            end

            # Perda radiativa pelas duas faces axiais do anel i (opcional)
            T_i_K = c_to_k(T[i])
            Q_end_i = include_ends ?
                (cfg.end_view_factor * cfg.emissivity_monolith * SIGMA_SB * A_ends_ring[i] *
                 max(0.0, T_i_K^4 - T_amb_K^4)) : 0.0

            # Demanda química proporcional ao volume do anel i avaliada na temperatura local T[i]
            Q_chem_i = (V_cell[i] / V_total) * (q_gas_W(chem, T[i]) + q_reac_W(chem, T[i]))

            C_i = m_cell[i] * inconel_cp(clamp(T[i], 20.0, 1000.0))
            dT[i] = (P_cell[i] + Q_prev - Q_next - Q_end_i - Q_chem_i) / C_i
            Q_prev = Q_next
        end
    end

    u0 = fill(Float64(T0_C), N)
    prob = ODEProblem(rhs_1d!, u0, (0.0, Float64(t_end_s)))
    sol = solve(prob, Tsit5(); reltol=1e-9, abstol=1e-9, maxiters=10^6)

    T_final = sol.u[end]
    # Auditoria de conservação de fluxo térmico no estado final
    Q_lat_out = surface_heat_losses(m, cfg, T_final[N]; include_ends=false).Q_lat_total_W
    Q_ends_out = include_ends ? sum(
        cfg.end_view_factor * cfg.emissivity_monolith * SIGMA_SB * A_ends_ring[i] *
        max(0.0, c_to_k(T_final[i])^4 - T_amb_K^4) for i in 1:N
    ) : 0.0
    Q_chem_out = sum((V_cell[i] / V_total) * (q_gas_W(chem, T_final[i]) + q_reac_W(chem, T_final[i])) for i in 1:N)
    Q_out_total = Q_lat_out + Q_ends_out + Q_chem_out
    flux_balance_rel_err = P_in > 0.0 ? abs(P_in - Q_out_total) / P_in : abs(Q_out_total)

    # Gradiente numérico junto ao eixo central vs. borda externa [°C / mm]
    dT_dr_center_C_mm = (T_final[2] - T_final[1]) / (dr * 1000.0)
    dT_dr_rim_C_mm    = (T_final[N] - T_final[N - 1]) / (dr * 1000.0)
    T_avg_C = sum(T_final[i] * V_cell[i] for i in 1:N) / V_total

    return (
        r_mm = r_centers .* 1000.0,
        T_C = T_final,
        q_ind_MW_m3 = q_ind_vol .* 1.0e-6,
        T_center_C = T_final[1],
        T_rim_C = T_final[N],
        T_max_C = maximum(T_final),
        T_min_C = minimum(T_final),
        T_avg_C = T_avg_C,
        delta_T_rim_minus_center_C = T_final[N] - T_final[1],
        delta_T_max_min_C = maximum(T_final) - minimum(T_final),
        delta_eff_mm = delta_eff * 1000.0,
        P_ind_W = P_in,
        Q_lat_out_W = Q_lat_out,
        Q_ends_out_W = Q_ends_out,
        Q_chem_out_W = Q_chem_out,
        flux_balance_rel_err = flux_balance_rel_err,
        dT_dr_center_C_mm = dT_dr_center_C_mm,
        dT_dr_rim_C_mm = dT_dr_rim_C_mm,
    )
end