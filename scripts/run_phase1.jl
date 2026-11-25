const T_START = time()

using DRMInconel
using DataFrames
using CSV
using Plots

# Configuração visual padrão para os gráficos
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

const MOLAR_MASS_C_G_MOL = 12.011
const DB = default_thermo_db()
const GEOM = ReactorGeometry()

# -----------------------------------------------------------------------------
# 1. Varredura de temperatura (600–900 °C, passo 10 °C) e exportação para CSV
# -----------------------------------------------------------------------------
println("[1/7] Rodando varreduras de temperatura e exportando CSV de referência...")

T_sweep = collect(600.0:10.0:900.0)
scenarios = [
    scenario_screening_800(),
    scenario_stability_700(),
    scenario_benchmark_1to1(),
]

rows = NamedTuple[]
for feed in scenarios, allow_coke in (false, true), T_C in T_sweep
    res = solve_equilibrium(feed; T_C=T_C, allow_coke=allow_coke, geometry=GEOM, db=DB)

    F_CH4_in = get(res.n_in_mol_s, :CH4, 0.0)
    coke_yield = F_CH4_in > 0 ? (res.coke_mol_s / F_CH4_in) : 0.0
    coke_mg_h = res.coke_mol_s * MOLAR_MASS_C_G_MOL * 1000.0 * 3600.0

    F_gas_out = sum(
        F for (sp, F) in res.n_out_mol_s
        if get_species_data(sp; db=DB).phase === :gas
    )
    y_pct(sp) = 100.0 * get(res.n_out_mol_s, sp, 0.0) / F_gas_out

    push!(rows, (
        scenario = feed.name,
        coke_model = allow_coke ? "graphite_equilibrium" : "gas_phase_only",
        T_C = res.T_C,
        P_atm = res.P_atm,
        V_in_ml_min = total_flow_ml_min(feed),
        X_CH4_pct = 100.0 * res.X_CH4,
        X_CO2_pct = 100.0 * res.X_CO2,
        H2_CO_ratio = res.H2_CO_ratio,
        coke_mol_per_mol_CH4_in = coke_yield,
        coke_mg_h = coke_mg_h,
        Q_sensible_W = res.Q_sensible_W,
        Q_reaction_W = res.Q_reaction_W,
        Q_total_W = res.Q_total_W,
        V_out_ml_min = res.V_out_ml_min,
        expansion_factor = res.expansion_factor,
        gas_velocity_cm_s = res.gas_velocity_cm_s,
        residence_time_s = res.residence_time_s,
        y_H2_pct = y_pct(:H2),
        y_CO_pct = y_pct(:CO),
        y_CH4_pct = y_pct(:CH4),
        y_CO2_pct = y_pct(:CO2),
        y_H2O_pct = y_pct(:H2O),
        y_N2_pct = y_pct(:N2),
    ))
end

df_ref = DataFrame(rows)
csv_path = joinpath(OUTPUT_DIR, "equilibrium_reference.csv")
CSV.write(csv_path, df_ref)
println("      Salvo: $csv_path ($(nrow(df_ref)) linhas)")

# Resumo rápido dos pontos operacionais principais no terminal
println("\n=== Resumo dos pontos principais (Fase 1) ===")
key_points = filter(
    r -> (r.scenario == "screening_800" && r.T_C == 800.0) ||
         (r.scenario == "stability_700" && r.T_C == 700.0),
    df_ref,
)
for r in eachrow(key_points)
    println(
        "  ", rpad(r.scenario, 14),
        " | ", rpad(r.coke_model, 20),
        " | ", Int(r.T_C), " °C => ",
        "X_CH4=", round(r.X_CH4_pct, digits=2), "%",
        " | X_CO2=", round(r.X_CO2_pct, digits=2), "%",
        " | H2/CO=", round(r.H2_CO_ratio, digits=3),
        " | Q=", round(r.Q_total_W, digits=2), " W",
        " | coke=", round(r.coke_mg_h, digits=1), " mg/h"
    )
end
println()

slice_df(scen, model) = filter(r -> r.scenario == scen && r.coke_model == model, df_ref)

# -----------------------------------------------------------------------------
# 2. Figura 1 — Curvas de equilíbrio (600–900 °C) vs. pontos experimentais
# -----------------------------------------------------------------------------
println("[2/7] Gerando Figura 1 (Curvas de equilíbrio vs. experimental)...")

df_scr_gas  = slice_df("screening_800", "gas_phase_only")
df_scr_coke = slice_df("screening_800", "graphite_equilibrium")

p1_conv = plot(
    df_scr_gas.T_C, df_scr_gas.X_CO2_pct,
    label = "CO2 eq. (gas phase)",
    color = :steelblue,
    xlabel = "Temperature (°C)",
    ylabel = "Conversion (%)",
    title = "Equilibrium Ceiling vs. Experimental Points",
    ylims = (25, 102),
    legend = :bottomright,
)
plot!(p1_conv, df_scr_gas.T_C, df_scr_gas.X_CH4_pct,
    label = "CH4 eq. (gas phase)", color = :darkorange)
plot!(p1_conv, df_scr_coke.T_C, df_scr_coke.X_CO2_pct,
    label = "CO2 eq. (w/ graphite)", color = :steelblue, linestyle = :dash, alpha = 0.65)
plot!(p1_conv, df_scr_coke.T_C, df_scr_coke.X_CH4_pct,
    label = "CH4 eq. (w/ graphite)", color = :darkorange, linestyle = :dash, alpha = 0.65)

# Pontos experimentais de referência (800 °C triagem e 700 °C estabilidade 24h)
scatter!(p1_conv, [800.0], [95.0], label = "Exp. MON-Ind 800°C (>95% CO2, >90% CH4)",
    color = :forestgreen, markershape = :star5, markersize = 8)
scatter!(p1_conv, [800.0], [90.0], label = false,
    color = :forestgreen, markershape = :star5, markersize = 8)

scatter!(p1_conv, [800.0], [90.0], label = "Exp. MON-Conv 800°C (~90% CO2, ~80% CH4)",
    color = :royalblue, markershape = :diamond, markersize = 6)
scatter!(p1_conv, [800.0], [80.0], label = false,
    color = :royalblue, markershape = :diamond, markersize = 6)

scatter!(p1_conv, [700.0, 700.0], [80.0, 55.0], label = "Exp. 700°C 24h CO2 (80% -> 55%)",
    color = :steelblue, markershape = :circle, markersize = 5, markerstrokewidth = 1.5)
scatter!(p1_conv, [700.0, 700.0], [60.0, 39.0], label = "Exp. 700°C 24h CH4 (60% -> 39%)",
    color = :darkorange, markershape = :square, markersize = 5, markerstrokewidth = 1.5)

p1_ratio = plot(
    df_scr_gas.T_C, df_scr_gas.H2_CO_ratio,
    label = "H2/CO eq. (gas phase)",
    color = :purple,
    xlabel = "Temperature (°C)",
    ylabel = "H2 / CO Molar Ratio",
    title = "Syngas H2/CO Ratio (1 atm, 27.4% N2)",
    ylims = (0.6, 1.6),
    legend = :topright,
)
plot!(p1_ratio, df_scr_coke.T_C, df_scr_coke.H2_CO_ratio,
    label = "H2/CO eq. (w/ graphite)", color = :purple, linestyle = :dash, alpha = 0.65)
hline!(p1_ratio, [1.0], label = "Stoichiometric DRM (1.0)", color = :gray, linestyle = :dot)
scatter!(p1_ratio, [700.0], [0.8], label = "Exp. 700°C 24h (~0.8)",
    color = :crimson, markershape = :circle, markersize = 7)

fig1 = plot(p1_conv, p1_ratio, layout = (1, 2), size = (1300, 520))
fig1_path = joinpath(OUTPUT_DIR, "fig1_equilibrium_curves.png")
savefig(fig1, fig1_path)
println("      Salvo: $fig1_path")

# -----------------------------------------------------------------------------
# 3. Figura 2 — Mapa de deposição de carbono (T x razão CO2/CH4)
# -----------------------------------------------------------------------------
println("[3/7] Gerando Figura 2 (Mapa da zona livre de coque)...")

T_grid = collect(600.0:15.0:900.0)
ratio_grid = collect(0.5:0.1:2.5)  # razão molar CO2 / CH4

function compute_coke_map(T_vals, ratio_vals; N2_frac=0.28, O2_to_CH4=0.0)
    Z = zeros(Float64, length(ratio_vals), length(T_vals))
    for (j, ratio) in enumerate(ratio_vals)
        CH4_flow = 50.0
        CO2_flow = ratio * CH4_flow
        O2_flow  = O2_to_CH4 * CH4_flow
        reactive_flow = CH4_flow + CO2_flow + O2_flow
        N2_flow = N2_frac > 0 ? reactive_flow * (N2_frac / (1.0 - N2_frac)) : 0.0

        feed = scenario_custom(CH4=CH4_flow, CO2=CO2_flow, N2=N2_flow, O2=O2_flow)
        for (i, T_C) in enumerate(T_vals)
            res = solve_equilibrium(feed; T_C=T_C, allow_coke=true, geometry=GEOM, db=DB)
            Z[j, i] = res.coke_mol_s / res.n_in_mol_s[:CH4]
        end
    end
    return Z
end

Z_no_O2   = compute_coke_map(T_grid, ratio_grid; N2_frac=0.28, O2_to_CH4=0.0)
Z_with_O2 = compute_coke_map(T_grid, ratio_grid; N2_frac=0.28, O2_to_CH4=0.2)
c_max = max(maximum(Z_no_O2), maximum(Z_with_O2))

p2_base = heatmap(
    T_grid, ratio_grid, Z_no_O2,
    xlabel = "Temperature (°C)",
    ylabel = "Inlet CO2 / CH4 Molar Ratio",
    title = "DRM Feed (28% N2, No O2)",
    colorbar_title = "mol C(s) / mol CH4 in",
    c = :inferno,
    clims = (0.0, c_max),
)
contour!(p2_base, T_grid, ratio_grid, Z_no_O2,
    levels = [0.005], color = :cyan, linestyle = :dash, linewidth = 2.2)
scatter!(p2_base, [700.0, 800.0], [1.0, 1.0],
    label = "Exp. points (700 & 800°C)",
    color = :lime, markershape = :circle, markersize = 6, legend = :topright)

p2_o2 = heatmap(
    T_grid, ratio_grid, Z_with_O2,
    xlabel = "Temperature (°C)",
    ylabel = "Inlet CO2 / CH4 Molar Ratio",
    title = "Exploratory Co-Feed (28% N2, O2/CH4 = 0.2)",
    colorbar_title = "mol C(s) / mol CH4 in",
    c = :inferno,
    clims = (0.0, c_max),
)
contour!(p2_o2, T_grid, ratio_grid, Z_with_O2,
    levels = [0.005], color = :cyan, linestyle = :dash, linewidth = 2.2)
scatter!(p2_o2, [700.0, 800.0], [1.0, 1.0],
    label = "Ref. (700 & 800°C, 1:1)",
    color = :lime, markershape = :circle, markersize = 6, legend = :topright)

fig2 = plot(p2_base, p2_o2, layout = (1, 2), size = (1380, 500))
fig2_path = joinpath(OUTPUT_DIR, "fig2_carbon_boundary.png")
savefig(fig2, fig2_path)
println("      Salvo: $fig2_path")

# -----------------------------------------------------------------------------
# 4. Figura 3 — Demanda mínima de potência térmica (Watts)
# -----------------------------------------------------------------------------
println("[4/7] Gerando Figura 3 (Demanda mínima de potência térmica)...")

df_stab_gas = slice_df("stability_700", "gas_phase_only")

p3_scr = plot(
    df_scr_gas.T_C, df_scr_gas.Q_total_W,
    label = "Total Power (Q_total)", color = :black,
    xlabel = "Temperature (°C)",
    ylabel = "Minimum Thermal Power (W)",
    title = "Screening Feed (62 mL/min, 27.4% N2)",
    legend = :topleft,
)
plot!(p3_scr, df_scr_gas.T_C, df_scr_gas.Q_reaction_W,
    label = "Reaction Enthalpy (Q_reaction)", color = :crimson, linestyle = :dash)
plot!(p3_scr, df_scr_gas.T_C, df_scr_gas.Q_sensible_W,
    label = "Gas Heating 25°C -> T (Q_sensible)", color = :steelblue, linestyle = :dot)
vline!(p3_scr, [800.0], label = "Screening T (800°C)", color = :forestgreen, linestyle = :dashdot)

p3_stab = plot(
    df_stab_gas.T_C, df_stab_gas.Q_total_W,
    label = "Total Power (Q_total)", color = :black,
    xlabel = "Temperature (°C)",
    ylabel = "Minimum Thermal Power (W)",
    title = "Stability Feed (200 mL/min, 28.0% N2)",
    legend = :topleft,
)
plot!(p3_stab, df_stab_gas.T_C, df_stab_gas.Q_reaction_W,
    label = "Reaction Enthalpy (Q_reaction)", color = :crimson, linestyle = :dash)
plot!(p3_stab, df_stab_gas.T_C, df_stab_gas.Q_sensible_W,
    label = "Gas Heating 25°C -> T (Q_sensible)", color = :steelblue, linestyle = :dot)
vline!(p3_stab, [700.0], label = "Stability T (700°C)", color = :darkorange, linestyle = :dashdot)

fig3 = plot(p3_scr, p3_stab, layout = (1, 2), size = (1300, 500))
fig3_path = joinpath(OUTPUT_DIR, "fig3_thermal_power_demand.png")
savefig(fig3, fig3_path)
println("      Salvo: $fig3_path")

# -----------------------------------------------------------------------------
# 5. Figura 4 — Sensibilidade à razão molar CO2 / CH4 (0.5 a 2.5, com 28% N2)
# -----------------------------------------------------------------------------
println("[5/7] Rodando sensibilidade à razão CO2/CH4 (Figura 4 + CSV)...")

ratios = collect(0.5:0.1:2.5)
rows_ratio = NamedTuple[]

for T_C in (700.0, 800.0), ratio in ratios
    V_tot = 100.0
    N2_flow = 0.28 * V_tot
    reactive = V_tot - N2_flow
    CH4_flow = reactive / (1.0 + ratio)
    CO2_flow = ratio * CH4_flow

    feed = scenario_custom(CH4=CH4_flow, CO2=CO2_flow, N2=N2_flow)
    r_gas  = solve_equilibrium(feed; T_C=T_C, allow_coke=false, geometry=GEOM, db=DB)
    r_coke = solve_equilibrium(feed; T_C=T_C, allow_coke=true,  geometry=GEOM, db=DB)

    push!(rows_ratio, (
        T_C = T_C,
        CO2_CH4_ratio = ratio,
        CH4_ml_min = CH4_flow,
        CO2_ml_min = CO2_flow,
        N2_ml_min = N2_flow,
        X_CH4_gas_pct = 100.0 * r_gas.X_CH4,
        X_CO2_gas_pct = 100.0 * r_gas.X_CO2,
        H2_CO_gas = r_gas.H2_CO_ratio,
        coke_mol_per_mol_CH4 = r_coke.coke_mol_s / r_coke.n_in_mol_s[:CH4],
        Q_total_gas_W = r_gas.Q_total_W,
    ))
end

df_ratio = DataFrame(rows_ratio)
csv_ratio_path = joinpath(OUTPUT_DIR, "sensitivity_co2_ch4.csv")
CSV.write(csv_ratio_path, df_ratio)

df_r700 = filter(r -> r.T_C == 700.0, df_ratio)
df_r800 = filter(r -> r.T_C == 800.0, df_ratio)

pa1 = plot(df_r800.CO2_CH4_ratio, df_r800.X_CH4_gas_pct,
    label="CH4 (800°C)", color=:darkorange,
    xlabel="Inlet CO2 / CH4 Molar Ratio", ylabel="Gas-Phase Conversion (%)",
    title="Conversion vs. CO2/CH4", ylims=(20, 102), legend=:bottomright)
plot!(pa1, df_r800.CO2_CH4_ratio, df_r800.X_CO2_gas_pct,
    label="CO2 (800°C)", color=:steelblue)
plot!(pa1, df_r700.CO2_CH4_ratio, df_r700.X_CH4_gas_pct,
    label="CH4 (700°C)", color=:darkorange, linestyle=:dash)
plot!(pa1, df_r700.CO2_CH4_ratio, df_r700.X_CO2_gas_pct,
    label="CO2 (700°C)", color=:steelblue, linestyle=:dash)
vline!(pa1, [1.0], label="1:1 Ref.", color=:gray, linestyle=:dot)

pa2 = plot(df_r800.CO2_CH4_ratio, df_r800.H2_CO_gas,
    label="H2/CO (800°C)", color=:purple,
    xlabel="Inlet CO2 / CH4 Molar Ratio", ylabel="H2 / CO Molar Ratio",
    title="Syngas Quality vs. CO2/CH4", legend=:topright)
plot!(pa2, df_r700.CO2_CH4_ratio, df_r700.H2_CO_gas,
    label="H2/CO (700°C)", color=:purple, linestyle=:dash)
hline!(pa2, [1.0], label="Stoich. (1.0)", color=:gray, linestyle=:dot)
vline!(pa2, [1.0], label=false, color=:gray, linestyle=:dot)

pa3 = plot(df_r800.CO2_CH4_ratio, df_r800.coke_mol_per_mol_CH4,
    label="Coke (800°C)", color=:crimson,
    xlabel="Inlet CO2 / CH4 Molar Ratio", ylabel="mol C(s) / mol CH4 in",
    title="Thermodynamic Coke vs. CO2/CH4", legend=:topright)
plot!(pa3, df_r700.CO2_CH4_ratio, df_r700.coke_mol_per_mol_CH4,
    label="Coke (700°C)", color=:crimson, linestyle=:dash)
vline!(pa3, [1.0], label="1:1 Ref.", color=:gray, linestyle=:dot)

fig4 = plot(pa1, pa2, pa3, layout=(1, 3), size=(1550, 460))
fig4_path = joinpath(OUTPUT_DIR, "fig4_sensitivity_co2_ch4.png")
savefig(fig4, fig4_path)
println("      Salvo: $csv_ratio_path e$fig4_path")

# -----------------------------------------------------------------------------
# 6. Figura 5 — Sensibilidade à diluição com N2 (0% a 60%, 1:1, 100 mL/min)
# -----------------------------------------------------------------------------
println("[6/7] Rodando sensibilidade à diluição de N2 (Figura 5 + CSV)...")

n2_fracs = collect(0.0:0.05:0.60)
rows_n2 = NamedTuple[]

for T_C in (700.0, 800.0), f_n2 in n2_fracs
    V_tot = 100.0
    N2_flow = f_n2 * V_tot
    CH4_flow = 0.5 * (V_tot - N2_flow)
    CO2_flow = CH4_flow

    feed = scenario_custom(CH4=CH4_flow, CO2=CO2_flow, N2=N2_flow)
    r_gas  = solve_equilibrium(feed; T_C=T_C, allow_coke=false, geometry=GEOM, db=DB)
    r_coke = solve_equilibrium(feed; T_C=T_C, allow_coke=true,  geometry=GEOM, db=DB)

    push!(rows_n2, (
        T_C = T_C,
        N2_pct = 100.0 * f_n2,
        X_CH4_gas_pct = 100.0 * r_gas.X_CH4,
        X_CO2_gas_pct = 100.0 * r_gas.X_CO2,
        H2_CO_gas = r_gas.H2_CO_ratio,
        coke_mol_per_mol_CH4 = r_coke.coke_mol_s / r_coke.n_in_mol_s[:CH4],
        Q_total_gas_W = r_gas.Q_total_W,
        residence_time_s = r_gas.residence_time_s,
    ))
end

df_n2 = DataFrame(rows_n2)
csv_n2_path = joinpath(OUTPUT_DIR, "sensitivity_n2_dilution.csv")
CSV.write(csv_n2_path, df_n2)

df_n700 = filter(r -> r.T_C == 700.0, df_n2)
df_n800 = filter(r -> r.T_C == 800.0, df_n2)

elapsed = round(time() - T_START, digits=2)
println("\nConcluído em ", elapsed, " s! Todos os 4 CSVs e 6 figuras foram salvos em: ", OUTPUT_DIR)