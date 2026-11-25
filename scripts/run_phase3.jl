const T_START = time()

using DRMInconel
using DataFrames
using CSV
using Plots

default(
    fontfamily = "Computer Modern",
    linewidth = 2.2,
    framestyle = :box,
    grid = true,
    gridalpha = 0.25,
    dpi = 300,
    margin = 5Plots.mm,
)

const OUTPUT_DIR = normpath(joinpath(@__DIR__, "..", "output"))
mkpath(OUTPUT_DIR)

const M_EXP = experimental_monolith(contact_factor=0.5)
const CFG_CLEAN = ThermalReactorConfig(surface_state=:clean)
const CFG_OX    = ThermalReactorConfig(surface_state=:oxidized)
const CFG_COKED = ThermalReactorConfig(surface_state=:coked)
const CFG_INS   = ThermalReactorConfig(surface_state=:oxidized, insulation_thickness_m=0.015)

println("[1/4] Pre-computando curvas de Q_gas(T) e Q_reac(T) a partir da Fase 1...")
const CHEM_SCR  = build_chemistry_load(scenario_screening_800(); allow_coke=false)
const CHEM_STAB = build_chemistry_load(scenario_stability_700(); allow_coke=false)

# -----------------------------------------------------------------------------
# Gráfico 1 — Curvas de Aquecimento T(t) (Potência, Isolamento e Emissividade)
# -----------------------------------------------------------------------------
println("[2/4] Gerando Grafico 1 (Curvas de Aquecimento T(t) via OrdinaryDiffEq)...")

p1_left = plot(
    xlabel = "Time (s)",
    ylabel = "Monolith Temperature (°C)",
    title = "Heating Ramps: Power & Insulation (Oxidized, 62 mL/min)",
    ylims = (0, 960),
    legend = :bottomright,
)
for (P_W, col, ls, cfg, lbl) in (
    (60.0,  :steelblue,   :solid, CFG_OX,  "60 W (Bare Quartz)"),
    (115.0, :darkorange,  :solid, CFG_OX,  "115 W (Bare Quartz)"),
    (180.0, :crimson,     :solid, CFG_OX,  "180 W (Bare Quartz)"),
    (60.0,  :forestgreen, :dash,  CFG_INS, "60 W (15 mm Insulation)"),
    (115.0, :purple,      :dash,  CFG_INS, "115 W (15 mm Insulation)"),
)
    sim = simulate_ramp_0d(M_EXP, cfg; P_ind_W=P_W, t_end_s=150.0, chem=CHEM_SCR)
    plot!(p1_left, sim.t_s, sim.T_C, label = lbl, color = col, linestyle = ls)
end
hline!(p1_left, [700.0, 800.0], label = "DRM Targets (700 & 800°C)", color = :gray, linestyle = :dot)

p1_right = plot(
    xlabel = "Time (s)",
    ylabel = "Monolith Temperature (°C)",
    title = "Surface Emissivity Impact at Constant P_ind = 115 W",
    ylims = (0, 960),
    legend = :bottomright,
)
for (cfg, col, lbl) in (
    (CFG_CLEAN, :steelblue,  "Clean Inconel (eps = 0.30)"),
    (CFG_OX,    :darkorange, "Oxidized / Activated (eps = 0.75)"),
    (CFG_COKED, :crimson,    "Coked Surface (eps = 0.90)"),
)
    sim = simulate_ramp_0d(M_EXP, cfg; P_ind_W=115.0, t_end_s=150.0, chem=CHEM_SCR)
    plot!(p1_right, sim.t_s, sim.T_C, label = lbl, color = col)
end
hline!(p1_right, [800.0], label = "800°C Reference", color = :gray, linestyle = :dot)

fig1 = plot(p1_left, p1_right, layout = (1, 2), size = (1420, 490))
fig1_path = joinpath(OUTPUT_DIR, "fig_p3_1_heating_ramps.png")
savefig(fig1, fig1_path)
println("      Salvo: ", fig1_path)

# -----------------------------------------------------------------------------
# Gráfico 2 — Raio-X das Perdas Energéticas (Decomposição % e Watts)
# -----------------------------------------------------------------------------
println("[3/4] Gerando Grafico 2 (Raio-X das Perdas Energeticas)...")

T_sweep = collect(600.0:25.0:900.0)
rows = NamedTuple[]
for (scen, chem) in (("screening_62", CHEM_SCR), ("stability_200", CHEM_STAB)),
    (cfg_lbl, cfg) in (("bare_oxidized", CFG_OX), ("insulated_15mm", CFG_INS)),
    T_C in T_sweep
    bd = steady_state_power_0d(M_EXP, cfg, T_C; chem=chem)
    push!(rows, (
        scenario = scen, config = cfg_lbl, T_C = T_C,
        P_ind_W = round(bd.P_ind_W, digits=2),
        Q_rad_W = round(bd.Q_rad_W, digits=2),
        Q_conv_W = round(bd.Q_conv_W, digits=2),
        Q_gas_W = round(bd.Q_gas_W, digits=3),
        Q_reac_W = round(bd.Q_reac_W, digits=3),
        pct_rad = round(bd.pct_rad, digits=2),
        pct_conv = round(bd.pct_conv, digits=2),
        pct_gas = round(bd.pct_gas, digits=2),
        pct_reac = round(bd.pct_reac, digits=2),
        T_q_C = round(bd.T_q_C, digits=1),
        T_outer_C = round(bd.T_outer_C, digits=1),
    ))
end
df_p3 = DataFrame(rows)
csv_p3 = joinpath(OUTPUT_DIR, "phase3_energy_breakdown.csv")
CSV.write(csv_p3, df_p3)

println("\n=== Raio-X Energetico Estacionario (700 °C e 800 °C) ===")
for (scen, T_mark) in (("screening_62", 800.0), ("stability_200", 700.0))
    for cfg_lbl in ("bare_oxidized", "insulated_15mm")
        r = only(filter(x -> x.scenario == scen && x.config == cfg_lbl && x.T_C == T_mark, df_p3))
        println(
            "  ", rpad(scen, 14), " | ", rpad(cfg_lbl, 15), " | ", Int(T_mark), "°C => ",
            "P_ind=", r.P_ind_W, " W | ",
            "Rad=", r.pct_rad, "% | Conv=", r.pct_conv, "% | ",
            "Gas=", r.pct_gas, "% | Reac=", r.pct_reac, "% | ",
            "T_out=", r.T_outer_C, "°C"
        )
    end
end

# Painel esquerdo: barras empilhadas via somas acumuladas (fundo -> frente no GR)
case_labels = ["62 mL/min\n800°C Bare", "62 mL/min\n800°C Insul.", "200 mL/min\n700°C Bare", "200 mL/min\n700°C Insul."]
case_rows = [
    only(filter(x -> x.scenario == "screening_62"  && x.config == "bare_oxidized"  && x.T_C == 800.0, df_p3)),
    only(filter(x -> x.scenario == "screening_62"  && x.config == "insulated_15mm" && x.T_C == 800.0, df_p3)),
    only(filter(x -> x.scenario == "stability_200" && x.config == "bare_oxidized"  && x.T_C == 700.0, df_p3)),
    only(filter(x -> x.scenario == "stability_200" && x.config == "insulated_15mm" && x.T_C == 700.0, df_p3)),
]
pct_reac_v = [r.pct_reac for r in case_rows]
pct_gas_v  = [r.pct_gas  for r in case_rows]
pct_conv_v = [r.pct_conv for r in case_rows]
pct_rad_v  = [r.pct_rad  for r in case_rows]

cum_reac = pct_reac_v
cum_gas  = cum_reac .+ pct_gas_v
cum_conv = cum_gas  .+ pct_conv_v
cum_rad  = cum_conv .+ pct_rad_v

x_idx = 1:4
p2_left = bar(x_idx, cum_rad, label = "Quartz Radiation", color = :crimson,
    xticks = (x_idx, case_labels), ylims = (0, 108),
    ylabel = "Share of Induced Power (%)", title = "Energy X-Ray: Where Each Watt Goes",
    legend = :outertopright)
bar!(p2_left, x_idx, cum_conv, label = "Natural Convection", color = :gold)
bar!(p2_left, x_idx, cum_gas,  label = "Gas Heating (Q_gas)", color = :steelblue)
bar!(p2_left, x_idx, cum_reac, label = "DRM Reaction (Q_reac)", color = :forestgreen)

# Painel direito: potência total necessária vs. T (600-900 °C)
df_scr_bare = filter(x -> x.scenario == "screening_62" && x.config == "bare_oxidized", df_p3)
df_scr_ins  = filter(x -> x.scenario == "screening_62" && x.config == "insulated_15mm", df_p3)
df_stb_bare = filter(x -> x.scenario == "stability_200" && x.config == "bare_oxidized", df_p3)

p2_right = plot(
    df_scr_bare.T_C, df_scr_bare.P_ind_W,
    label = "Total P_ind - Bare (62 mL/min)", color = :crimson,
    xlabel = "Monolith Temperature (°C)",
    ylabel = "Required Induced Power (W)",
    title = "Induced Power vs. Temperature (T^4 Radiation Effect)",
    legend = :topleft,
)
plot!(p2_right, df_stb_bare.T_C, df_stb_bare.P_ind_W,
    label = "Total P_ind - Bare (200 mL/min)", color = :darkorange, linestyle = :dash)
plot!(p2_right, df_scr_ins.T_C, df_scr_ins.P_ind_W,
    label = "Total P_ind - 15 mm Insulation (62 mL/min)", color = :forestgreen)
plot!(p2_right, df_scr_bare.T_C, df_scr_bare.Q_gas_W .+ df_scr_bare.Q_reac_W,
    label = "Chemistry Only (62 mL/min)", color = :steelblue, linestyle = :dot, linewidth = 2.6)

fig2 = plot(p2_left, p2_right, layout = (1, 2), size = (1460, 490))
fig2_path = joinpath(OUTPUT_DIR, "fig_p3_2_energy_xray.png")
savefig(fig2, fig2_path)
println("      Salvo: ", csv_p3, " e ", fig2_path)

# -----------------------------------------------------------------------------
# Gráfico 3 — Perfil Térmico Radial 1D T(r) e Geração Indutiva q_ind(r)
# -----------------------------------------------------------------------------
println("[4/4] Gerando Grafico 3 (Perfil Termico Radial 1D via Metodo das Linhas)...")

P_800_bare = steady_state_power_0d(M_EXP, CFG_OX,  800.0; chem=CHEM_SCR).P_ind_W
P_800_ins  = steady_state_power_0d(M_EXP, CFG_INS, 800.0; chem=CHEM_SCR).P_ind_W

rad_50  = simulate_radial_1d(M_EXP, CFG_OX;  P_ind_W=P_800_bare, f_hz=50e3,  N=30, t_end_s=250.0, chem=CHEM_SCR)
rad_200 = simulate_radial_1d(M_EXP, CFG_OX;  P_ind_W=P_800_bare, f_hz=200e3, N=30, t_end_s=250.0, chem=CHEM_SCR)
rad_400 = simulate_radial_1d(M_EXP, CFG_OX;  P_ind_W=P_800_bare, f_hz=400e3, N=30, t_end_s=250.0, chem=CHEM_SCR)
rad_ins = simulate_radial_1d(M_EXP, CFG_INS; P_ind_W=P_800_ins,  f_hz=200e3, N=30, t_end_s=800.0, chem=CHEM_SCR)

println("\n=== Perfil Radial 1D em ~800 °C (Centro r=0 vs. Borda r=10 mm) ===")
for (lbl, r) in (
    ("Bare @ 50 kHz ", rad_50),
    ("Bare @ 200 kHz", rad_200),
    ("Bare @ 400 kHz", rad_400),
    ("Insul @ 200kHz", rad_ins),
)
    println(
        "  ", lbl,
        " | T_center=", round(r.T_center_C, digits=1), "°C",
        " | T_max=", round(r.T_max_C, digits=1), "°C",
        " | T_rim=", round(r.T_rim_C, digits=1), "°C",
        " | DeltaT=", round(r.delta_T_max_min_C, digits=1), "°C",
        " | err_flux=", round(r.flux_balance_rel_err, sigdigits=2)
    )
end

p3_left = plot(
    rad_50.r_mm, rad_50.T_C,
    label = "Bare Quartz @ 50 kHz (deep skin)", color = :steelblue,
    xlabel = "Radial Position, r (mm)  [0 = Center, 10 = Outer Rim]",
    ylabel = "Monolith Temperature, T(r) (°C)",
    title = "1D Radial Temperature Profile T(r) (~800°C Target)",
    legend = :left,
)
plot!(p3_left, rad_200.r_mm, rad_200.T_C, label = "Bare Quartz @ 200 kHz", color = :darkorange)
plot!(p3_left, rad_400.r_mm, rad_400.T_C, label = "Bare Quartz @ 400 kHz (shallow skin)", color = :crimson)
plot!(p3_left, rad_ins.r_mm, rad_ins.T_C, label = "15 mm Insulation @ 200 kHz", color = :forestgreen, linestyle = :dash)

p3_right = plot(
    rad_50.r_mm, rad_50.q_ind_MW_m3,
    label = "50 kHz (delta_eff = " * string(round(rad_50.delta_eff_mm, digits=1)) * " mm)",
    color = :steelblue,
    xlabel = "Radial Position, r (mm)",
    ylabel = "Volumetric Inductive Heating (MW / m3)",
    title = "Spatial Distribution of Induced Heat q_ind(r)",
    legend = :topleft,
)
plot!(p3_right, rad_200.r_mm, rad_200.q_ind_MW_m3,
    label = "200 kHz (delta_eff = " * string(round(rad_200.delta_eff_mm, digits=1)) * " mm)",
    color = :darkorange)
plot!(p3_right, rad_400.r_mm, rad_400.q_ind_MW_m3,
    label = "400 kHz (delta_eff = " * string(round(rad_400.delta_eff_mm, digits=1)) * " mm)",
    color = :crimson)

fig3 = plot(p3_left, p3_right, layout = (1, 2), size = (1420, 490))
fig3_path = joinpath(OUTPUT_DIR, "fig_p3_3_radial_profile.png")
savefig(fig3, fig3_path)
println("\n      Salvo: ", fig3_path)

elapsed = round(time() - T_START, digits=2)
println("Concluido em ", elapsed, " s!")