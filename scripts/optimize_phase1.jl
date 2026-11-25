const T_START = time()

using DRMInconel
using DataFrames
using CSV
using Plots

default(
    fontfamily = "Computer Modern",
    linewidth = 2.0,
    framestyle = :box,
    grid = true,
    gridalpha = 0.25,
    dpi = 300,
    margin = 5Plots.mm,
)

const OUTPUT_DIR = normpath(joinpath(@__DIR__, "..", "output"))
mkpath(OUTPUT_DIR)

const DB = default_thermo_db()
const GEOM = ReactorGeometry()
const MOLAR_MASS_C_G_MOL = 12.011
const TOTAL_FLOW_ML_MIN = 100.0  # Base fixa (100 mL/min) para comparacao justa

# -----------------------------------------------------------------------------
# 1. Gerador de composicoes (simplex 3 gases sem O2 + 4 gases com O2)
# -----------------------------------------------------------------------------
function build_compositions()
    comps = NamedTuple{(:mode, :CH4, :CO2, :N2, :O2), Tuple{String, Float64, Float64, Float64, Float64}}[]
    seen = Set{NTuple{4, Int}}()

    function add_comp!(mode::String, ch4::Float64, co2::Float64, n2::Float64, o2::Float64)
        s = ch4 + co2 + n2 + o2
        c1, c2, c3, c4 = 100.0 * ch4 / s, 100.0 * co2 / s, 100.0 * n2 / s, 100.0 * o2 / s
        key = (round(Int, c1 * 10), round(Int, c2 * 10), round(Int, c3 * 10), round(Int, c4 * 10))
        if !(key in seen)
            push!(seen, key)
            push!(comps, (mode=mode, CH4=c1, CO2=c2, N2=c3, O2=c4))
        end
    end

    # --- A) Sem O2 (3 gases: CH4, CO2, N2 somando 100%) ---
    # Grade simplex de passo 10% (cada gas >= 10%)
    for ch4 in 10.0:10.0:80.0, co2 in 10.0:10.0:(90.0 - ch4)
        n2 = 100.0 - ch4 - co2
        if n2 >= 10.0 - 1e-9
            add_comp!("no_O2", ch4, co2, n2, 0.0)
        end
    end
    # Extremos (98/1/1), igualdade (33.3/33.3/33.3), bordas com 1% N2 e ref. experimental
    special_no_o2 = [
        (98.0,  1.0,  1.0),
        ( 1.0, 98.0,  1.0),
        ( 1.0,  1.0, 98.0),
        (100.0/3, 100.0/3, 100.0/3),
        (49.5, 49.5,  1.0),
        (49.5,  1.0, 49.5),
        ( 1.0, 49.5, 49.5),
        (45.0, 54.0,  1.0),
        (40.0, 59.0,  1.0),
        (36.0, 36.0, 28.0),
    ]
    for (ch4, co2, n2) in special_no_o2
        add_comp!("no_O2", ch4, co2, n2, 0.0)
    end

    # --- B) Com O2 (4 gases: CH4, CO2, N2, O2 somando 100%) ---
    # Fracoes relativas para distribuir (100 - O2) entre CH4, CO2 e N2
    rel_templates = [
        (0.98, 0.01, 0.01),
        (0.01, 0.98, 0.01),
        (0.01, 0.01, 0.98),
        (1.0/3, 1.0/3, 1.0/3),
        (0.495, 0.495, 0.01),
        (0.60, 0.39, 0.01),
        (0.70, 0.29, 0.01),
        (0.39, 0.60, 0.01),
        (0.36, 0.36, 0.28),
        (0.48, 0.24, 0.28),
        (0.58, 0.14, 0.28),
        (0.24, 0.48, 0.28),
    ]
    for o2 in (5.0, 10.0, 15.0, 20.0, 30.0)
        rem = 100.0 - o2
        for (f1, f2, f3) in rel_templates
            add_comp!("with_O2", rem * f1, rem * f2, rem * f3, o2)
        end
    end
    # Extremos 4D: 97% O2 + 1% do resto, e igualdade exata dos 4 gases
    add_comp!("with_O2", 1.0, 1.0, 1.0, 97.0)
    add_comp!("with_O2", 25.0, 25.0, 25.0, 25.0)

    return comps
end

# -----------------------------------------------------------------------------
# 2. Varredura global (600 a 1200 °C)
#    Dica: troque 50.0 por 100.0 abaixo se quiser rodar na metade do tempo
# -----------------------------------------------------------------------------
T_GRID = collect(600.0:50.0:1200.0)
COMPOSITIONS = build_compositions()
total_runs = length(COMPOSITIONS) * length(T_GRID)

println("[1/4] Iniciando varredura global: ", length(COMPOSITIONS), " composicoes x ",
        length(T_GRID), " temperaturas = ", total_runs, " simulacoes...")

rows = NamedTuple[]
n_failed = 0

for comp in COMPOSITIONS, T_C in T_GRID
    feed = scenario_custom(
        CH4 = TOTAL_FLOW_ML_MIN * comp.CH4 / 100.0,
        CO2 = TOTAL_FLOW_ML_MIN * comp.CO2 / 100.0,
        N2  = TOTAL_FLOW_ML_MIN * comp.N2  / 100.0,
        O2  = TOTAL_FLOW_ML_MIN * comp.O2  / 100.0,
    )

    try
        # allow_coke=true detecta formacao de carbono; onde coke==0 coincide com fase gasosa
        res = solve_equilibrium(feed; T_C=T_C, allow_coke=true, geometry=GEOM, db=DB)

        F_in_tot = sum(values(res.n_in_mol_s))
        F_C_in = get(res.n_in_mol_s, :CH4, 0.0) + get(res.n_in_mol_s, :CO2, 0.0)
        F_H2 = get(res.n_out_mol_s, :H2, 0.0)
        F_CO = get(res.n_out_mol_s, :CO, 0.0)
        F_syngas = F_H2 + F_CO

        syngas_per_mol_in = F_syngas / F_in_tot
        coke_pct_C = F_C_in > 0 ? 100.0 * res.coke_mol_s / F_C_in : 0.0
        coke_mg_h = res.coke_mol_s * MOLAR_MASS_C_G_MOL * 1000.0 * 3600.0
        h2_co = F_CO > 1e-10 ? (F_H2 / F_CO) : NaN
        kJ_per_mol_syn = F_syngas > 1e-9 ? (res.Q_total_W / F_syngas / 1000.0) : NaN

        F_gas_out = sum(
            F for (sp, F) in res.n_out_mol_s
            if get_species_data(sp; db=DB).phase === :gas
        )
        y_syngas_pct = 100.0 * F_syngas / F_gas_out

        push!(rows, (
            mode = comp.mode,
            T_C = T_C,
            CH4_pct = round(comp.CH4, digits=2),
            CO2_pct = round(comp.CO2, digits=2),
            N2_pct  = round(comp.N2,  digits=2),
            O2_pct  = round(comp.O2,  digits=2),
            X_CH4_pct = 100.0 * res.X_CH4,
            X_CO2_pct = 100.0 * res.X_CO2,
            syngas_per_mol_in = syngas_per_mol_in,
            y_syngas_out_pct = y_syngas_pct,
            H2_CO_ratio = h2_co,
            coke_pct_C = coke_pct_C,
            coke_mg_h = coke_mg_h,
            coke_free = coke_pct_C <= 0.1,
            Q_total_W = res.Q_total_W,
            kJ_per_mol_syngas = kJ_per_mol_syn,
            residence_time_s = res.residence_time_s,
        ))
    catch
        global n_failed += 1
    end
end

df_all = DataFrame(rows)
csv_all = joinpath(OUTPUT_DIR, "optimization_all_points.csv")
CSV.write(csv_all, df_all)
println("      Concluido: ", nrow(df_all), " convergidos, ", n_failed, " falhas em extremos.")
println("      Salvo: ", csv_all)

# -----------------------------------------------------------------------------
# 3. Extracao da Fronteira de Pareto (livre de coque + syngas de boa qualidade)
#    Objetivos simultaneos:
#      1) MAXIMIZAR syngas_per_mol_in
#      2) MINIMIZAR T_C
#      3) MINIMIZAR kJ_per_mol_syngas
# -----------------------------------------------------------------------------
println("[2/4] Calculando Fronteira de Pareto (livre de coque)...")

function pareto_front(df::DataFrame)
    n = nrow(df)
    keep = trues(n)
    syn = df.syngas_per_mol_in
    t   = df.T_C
    en  = df.kJ_per_mol_syngas

    for i in 1:n
        if !isfinite(en[i])
            keep[i] = false
            continue
        end
        for j in 1:n
            i == j && continue
            if !isfinite(en[j])
                continue
            end
            better_or_eq = (syn[j] >= syn[i] - 1e-9) &&
                           (t[j]   <= t[i]   + 1e-9) &&
                           (en[j]  <= en[i]  + 1e-9)
            strictly_better = (syn[j] > syn[i] + 1e-6) ||
                              (t[j]   < t[i]   - 1e-6) ||
                              (en[j]  < en[i]  - 1e-4)
            if better_or_eq && strictly_better
                keep[i] = false
                break
            end
        end
    end
    return sort(df[keep, :], [:mode, :T_C, :syngas_per_mol_in], rev=[false, false, true])
end

# Filtro de viabilidade: zero coque (<=0.1% C), H2/CO util e producao significativa
df_feas_no_o2 = filter(
    r -> r.mode == "no_O2" && r.coke_free &&
         isfinite(r.H2_CO_ratio) && 0.85 <= r.H2_CO_ratio <= 1.15 &&
         r.syngas_per_mol_in >= 0.3,
    df_all,
)
df_feas_with_o2 = filter(
    r -> r.mode == "with_O2" && r.coke_free &&
         isfinite(r.H2_CO_ratio) && 0.85 <= r.H2_CO_ratio <= 2.05 &&
         r.syngas_per_mol_in >= 0.3,
    df_all,
)

df_pareto_no_o2   = pareto_front(df_feas_no_o2)
df_pareto_with_o2 = pareto_front(df_feas_with_o2)
df_pareto = vcat(df_pareto_no_o2, df_pareto_with_o2)

csv_pareto = joinpath(OUTPUT_DIR, "optimization_pareto_front.csv")
CSV.write(csv_pareto, df_pareto)
println("      Pontos livres de coque (viaveis): Sem O2 = ", nrow(df_feas_no_o2),
        " | Com O2 = ", nrow(df_feas_with_o2))
println("      Pontos na Fronteira de Pareto:    Sem O2 = ", nrow(df_pareto_no_o2),
        " | Com O2 = ", nrow(df_pareto_with_o2))
println("      Salvo: ", csv_pareto)

function print_top(title::String, df::DataFrame; limit::Int=8)
    println("\n=== ", title, " ===")
    if nrow(df) == 0
        println("  (nenhum ponto satisfez os criterios)")
        return
    end
    top = first(sort(df, :syngas_per_mol_in, rev=true), min(limit, nrow(df)))
    for r in eachrow(top)
        println(
            "  T=", lpad(Int(r.T_C), 4), "°C",
            " | CH4/CO2/N2/O2=",
            lpad(round(r.CH4_pct, digits=1), 4), "/",
            lpad(round(r.CO2_pct, digits=1), 4), "/",
            lpad(round(r.N2_pct,  digits=1), 4), "/",
            lpad(round(r.O2_pct,  digits=1), 4), "%",
            " | Syngas=", round(r.syngas_per_mol_in, digits=3), " mol/mol",
            " | H2/CO=", round(r.H2_CO_ratio, digits=2),
            " | E=", round(r.kJ_per_mol_syngas, digits=1), " kJ/mol",
            " | Q=", round(r.Q_total_W, digits=2), " W"
        )
    end
end

print_top("TOP 8 SEM O2 (600-1200°C, livre de coque, 0.85<=H2/CO<=1.15)", df_pareto_no_o2)
print_top("TOP 8 COM O2 (600-1200°C, livre de coque, 0.85<=H2/CO<=2.05)", df_pareto_with_o2)
print_top("TOP 8 ZONA SEGURA INCONEL/QUARTZO (T <= 900°C, livre de coque)",
          filter(r -> r.T_C <= 900.0, df_pareto))

# -----------------------------------------------------------------------------
# 4. Figura 7 — Fronteira de Pareto (Produtividade vs. Temperatura e Energia)
# -----------------------------------------------------------------------------
println("\n[3/4] Gerando Figura 7 (Fronteira de Pareto)...")

df_coke_bad = filter(r -> !r.coke_free, df_all)

p7_left = scatter(
    df_coke_bad.T_C, df_coke_bad.syngas_per_mol_in,
    label = "Forms Coke (discarded)",
    color = :lightgray, alpha = 0.45, markersize = 3, markerstrokewidth = 0,
    xlabel = "Temperature (°C)",
    ylabel = "Syngas Yield (mol H2+CO / mol feed)",
    title = "Yield vs. Temperature (600-1200°C)",
    legend = :topleft,
)
scatter!(p7_left, df_feas_no_o2.T_C, df_feas_no_o2.syngas_per_mol_in,
    label = "Coke-free (No O2)", color = :steelblue, alpha = 0.45, markersize = 4, markerstrokewidth = 0)
scatter!(p7_left, df_feas_with_o2.T_C, df_feas_with_o2.syngas_per_mol_in,
    label = "Coke-free (With O2)", color = :darkorange, alpha = 0.45, markersize = 4, markerstrokewidth = 0)
scatter!(p7_left, df_pareto_no_o2.T_C, df_pareto_no_o2.syngas_per_mol_in,
    label = "Pareto Optimal (No O2)", color = :navy, markershape = :star5, markersize = 7)
scatter!(p7_left, df_pareto_with_o2.T_C, df_pareto_with_o2.syngas_per_mol_in,
    label = "Pareto Optimal (With O2)", color = :crimson, markershape = :diamond, markersize = 6)
vline!(p7_left, [900.0], label = "Quartz/Inconel practical limit (~900-1000°C)",
    color = :black, linestyle = :dash)

p7_right = scatter(
    df_feas_no_o2.kJ_per_mol_syngas, df_feas_no_o2.syngas_per_mol_in,
    label = "Coke-free (No O2)", color = :steelblue, alpha = 0.45, markersize = 4, markerstrokewidth = 0,
    xlabel = "Thermal Energy Cost (kJ / mol syngas)",
    ylabel = "Syngas Yield (mol H2+CO / mol feed)",
    title = "Yield vs. Specific Energy (Coke-Free Only)",
    legend = :bottomright,
)
scatter!(p7_right, df_feas_with_o2.kJ_per_mol_syngas, df_feas_with_o2.syngas_per_mol_in,
    label = "Coke-free (With O2)", color = :darkorange, alpha = 0.45, markersize = 4, markerstrokewidth = 0)
scatter!(p7_right, df_pareto_no_o2.kJ_per_mol_syngas, df_pareto_no_o2.syngas_per_mol_in,
    label = "Pareto Optimal (No O2)", color = :navy, markershape = :star5, markersize = 7)
scatter!(p7_right, df_pareto_with_o2.kJ_per_mol_syngas, df_pareto_with_o2.syngas_per_mol_in,
    label = "Pareto Optimal (With O2)", color = :crimson, markershape = :diamond, markersize = 6)
vline!(p7_right, [0.0], label = "Autothermal (0 kJ/mol)", color = :forestgreen, linestyle = :dot)

fig7 = plot(p7_left, p7_right, layout = (1, 2), size = (1380, 520))
fig7_path = joinpath(OUTPUT_DIR, "fig7_optimization_pareto.png")
savefig(fig7, fig7_path)
println("      Salvo: ", fig7_path)

# -----------------------------------------------------------------------------
# 5. Figura 8 — Mapa de composicoes sem O2 (CH4% vs CO2%, restante = N2%)
# -----------------------------------------------------------------------------
println("[4/4] Gerando Figura 8 (Mapa de composicoes a 800°C vs 1000°C)...")

function plot_composition_slice(T_target::Float64)
    sub = filter(r -> r.mode == "no_O2" && r.T_C == T_target, df_all)
    ok  = filter(r -> r.coke_free, sub)
    bad = filter(r -> !r.coke_free, sub)

    p = plot(
        [0.0, 100.0], [100.0, 0.0],
        label = "0% N2 boundary", color = :gray, linestyle = :dot,
        xlabel = "Inlet CH4 (%)",
        ylabel = "Inlet CO2 (%)",
        title = string("No-O2 Composition Map at ", Int(T_target), "°C (N2 = 100-CH4-CO2)"),
        xlims = (-2, 102), ylims = (-2, 102),
        legend = :topright,
    )
    scatter!(p, ok.CH4_pct, ok.CO2_pct,
        marker_z = ok.syngas_per_mol_in, c = :viridis, clims = (0.0, 2.0),
        markersize = 7, markerstrokewidth = 0.5,
        label = "Coke-free (color = syngas yield)",
        colorbar_title = "mol syngas / mol feed")
    scatter!(p, bad.CH4_pct, bad.CO2_pct,
        color = :crimson, markershape = :xcross, markersize = 6,
        label = "Forms Coke")
    scatter!(p, [36.0], [36.0],
        color = :lime, markershape = :star5, markersize = 9,
        label = "Exp. (36/36/28)")
    return p
end

fig8 = plot(plot_composition_slice(800.0), plot_composition_slice(1000.0),
            layout = (1, 2), size = (1380, 520))
fig8_path = joinpath(OUTPUT_DIR, "fig8_optimization_composition_map.png")
savefig(fig8, fig8_path)
println("      Salvo: ", fig8_path)

elapsed = round(time() - T_START, digits=2)
println("\nConcluido em ", elapsed, " s! Resultados salvos em: ", OUTPUT_DIR)