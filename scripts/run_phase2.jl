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

# -----------------------------------------------------------------------------
# 1. Gráfico 1 — Curvas de Skin Depth vs. Temperatura (10 kHz a 500 kHz)
# -----------------------------------------------------------------------------
println("[1/3] Gerando Grafico 1 (Skin Depth vs. Temperatura)...")

T_vals = collect(25.0:25.0:900.0)
freqs_kHz = [10.0, 30.0, 50.0, 100.0, 200.0, 300.0, 500.0]
colors = palette(:viridis, length(freqs_kHz))

p1_left = plot(
    xlabel = "Temperature (°C)",
    ylabel = "Skin Depth, delta (mm)",
    title = "Inconel 600 Skin Depth vs. Temperature",
    yscale = :log10,
    ylims = (0.03, 10.0),
    legend = :outertopright,
)
for (idx, fk) in enumerate(freqs_kHz)
    deltas = [skin_depth_mm(fk * 1000.0, T) for T in T_vals]
    plot!(p1_left, T_vals, deltas, label = string(Int(fk), " kHz"), color = colors[idx])
end
hline!(p1_left, [0.05], label = "Exp. foil (0.05 mm)", color = :crimson, linestyle = :dash, linewidth = 2.4)

# Painel direito: curvas de propriedades do Inconel 600 (resistividade e condutividade)
rho_uohm_m = [inconel_resistivity(T) * 1.0e6 for T in T_vals]
k_vals = [inconel_thermal_conductivity(T) for T in T_vals]
cp_vals = [inconel_cp(T) for T in T_vals]

p1_right = plot(
    T_vals, rho_uohm_m,
    label = "Resistivity rho (uOhm*m)",
    color = :navy,
    xlabel = "Temperature (°C)",
    ylabel = "Electrical Resistivity (uOhm*m)",
    title = "Inconel 600 Electrical Resistivity (Special Metals)",
    ylims = (1.00, 1.18),
    legend = :bottomright,
)
vline!(p1_right, [700.0, 800.0], label = "DRM Exp. (700 & 800°C)", color = :forestgreen, linestyle = :dot)

fig1 = plot(p1_left, p1_right, layout = (1, 2), size = (1380, 480))
fig1_path = joinpath(OUTPUT_DIR, "fig_p2_1_skin_depth.png")
savefig(fig1, fig1_path)
println("      Salvo: ", fig1_path)

# -----------------------------------------------------------------------------
# 2. Gráfico 2 — Mapa de Acoplamento (Frequência x Espessura/Diâmetro a 800 °C)
# -----------------------------------------------------------------------------
println("[2/3] Gerando Grafico 2 (Mapa de Acoplamento Indutivo a 800 °C)...")

T_ref = 800.0
f_grid_kHz = 10.0 .^ range(log10(10.0), log10(500.0), length=90)

# Mapa A: Placa/Parede isolada (espessura w de 0.03 mm a 10 mm)
w_grid_mm = 10.0 .^ range(log10(0.03), log10(10.0), length=90)
Z_slab = [slab_power_factor(w / skin_depth_mm(f * 1000.0, T_ref)) for w in w_grid_mm, f in f_grid_kHz]

p2_slab = heatmap(
    f_grid_kHz, w_grid_mm, Z_slab,
    xscale = :log10, yscale = :log10,
    xlabel = "Induction Frequency (kHz)",
    ylabel = "Isolated Wall / Foil Thickness, w (mm)",
    title = "Isolated Foil Coupling Factor F_slab (800°C)",
    colorbar_title = "F_slab (-)",
    c = :inferno, clims = (0.0, 1.1),
)
contour!(p2_slab, f_grid_kHz, w_grid_mm, Z_slab,
    levels = [0.5], color = :cyan, linestyle = :dash, linewidth = 2.0)
hline!(p2_slab, [0.05], label = "Exp. 50 um foil (transparent)",
    color = :lime, linestyle = :solid, linewidth = 2.4, legend = :bottomleft)

# Mapa B: Cilindro poroso macroscópico (D de 2 mm a 50 mm, eps=90.45%, contact=0.5)
D_grid_mm = 10.0 .^ range(log10(2.0), log10(50.0), length=90)
m_exp = experimental_monolith(contact_factor=0.5)
Z_cyl = zeros(Float64, length(D_grid_mm), length(f_grid_kHz))
for (j, D_mm) in enumerate(D_grid_mm), (i, f_kHz) in enumerate(f_grid_kHz)
    delta_eff_mm = effective_skin_depth_m(m_exp, f_kHz * 1000.0, T_ref) * 1000.0
    Z_cyl[j, i] = cylinder_power_factor(D_mm / delta_eff_mm)
end

p2_cyl = heatmap(
    f_grid_kHz, D_grid_mm, Z_cyl,
    xscale = :log10, yscale = :log10,
    xlabel = "Induction Frequency (kHz)",
    ylabel = "Monolith Outer Diameter, D (mm)",
    title = "Macroscopic Monolith Coupling F_cyl (800°C, 90.5% eps)",
    colorbar_title = "F_cyl (-)",
    c = :inferno, clims = (0.0, 1.0),
)
contour!(p2_cyl, f_grid_kHz, D_grid_mm, Z_cyl,
    levels = [0.5, 0.75], color = :cyan, linestyle = :dash, linewidth = 2.0)
hline!(p2_cyl, [20.0], label = "Exp. monolith (D = 20 mm)",
    color = :lime, linestyle = :solid, linewidth = 2.4, legend = :bottomleft)

fig2 = plot(p2_slab, p2_cyl, layout = (1, 2), size = (1420, 500))
fig2_path = joinpath(OUTPUT_DIR, "fig_p2_2_coupling_map.png")
savefig(fig2, fig2_path)
println("      Salvo: ", fig2_path)

# -----------------------------------------------------------------------------
# 3. Tabela e Gráfico 3 — Comparação de Geometrias (mesmo envelope 2 cm x 3 cm)
# -----------------------------------------------------------------------------
println("[3/3] Comparando geometrias de monolito (CSV + Grafico 3)...")

geometries = AbstractMonolith[
    experimental_monolith(contact_factor=0.5),
    RolledFoilMonolith(name="Rolled Foil (100 um, ideal)", foil_thickness_m=100e-6, contact_factor=0.8),
    HoneycombMonolith(name="Honeycomb 400 CPSI (150 um)", cpsi=400.0, wall_thickness_m=150e-6),
    HoneycombMonolith(name="Honeycomb 200 CPSI (250 um)", cpsi=200.0, wall_thickness_m=250e-6),
    MetalFoamMonolith(name="Metal Foam 30 PPI (90% eps)", ppi=30.0, porosity=0.90),
    PrintedLatticeMonolith(name="3D Lattice (2 mm / 0.4 mm)", cell_size_m=2.0e-3, strut_diameter_m=0.4e-3),
]

rows = NamedTuple[]
for g in geometries
    s = monolith_summary(g; T_C=20.0)
    c50  = analyze_monolith_coupling(g, 50e3,  800.0)
    c200 = analyze_monolith_coupling(g, 200e3, 800.0)
    c400 = analyze_monolith_coupling(g, 400e3, 800.0)

    push!(rows, (
        geometry = s.name,
        porosity_pct = round(s.porosity_pct, digits=2),
        wall_or_strut_um = round(s.wall_or_strut_um, digits=1),
        mass_g = round(s.mass_g, digits=2),
        surface_area_cm2 = round(s.surface_area_cm2, digits=1),
        Sv_m2_m3 = round(s.Sv_m2_m3, digits=0),
        dh_mm = round(s.hydraulic_diam_mm, digits=2),
        F_macro_50kHz = round(c50.F_macro_cyl, digits=3),
        F_macro_200kHz = round(c200.F_macro_cyl, digits=3),
        F_macro_400kHz = round(c400.F_macro_cyl, digits=3),
        eta_macro_200kHz_pct = round(c200.eta_macro_pct, digits=1),
        regime_200kHz = String(c200.macro_regime),
    ))
end

df_geom = DataFrame(rows)
csv_path = joinpath(OUTPUT_DIR, "phase2_geometry_comparison.csv")
CSV.write(csv_path, df_geom)

println("\n=== Tabela Comparativa de Geometrias (D = 2 cm, L = 3 cm, Inconel 600 @ 800 °C) ===")
for r in eachrow(df_geom)
    println(
        "  ", rpad(r.geometry, 30),
        " | eps=", lpad(r.porosity_pct, 5), "%",
        " | m=", lpad(r.mass_g, 5), " g",
        " | Sv=", lpad(Int(r.Sv_m2_m3), 4), " m2/m3",
        " | dh=", lpad(r.dh_mm, 4), " mm",
        " | F(50/200/400kHz)=", r.F_macro_50kHz, "/", r.F_macro_200kHz, "/", r.F_macro_400kHz,
        " | ", r.regime_200kHz
    )
end

# Gráfico 3: Curvas de acoplamento macroscópico vs. frequência + trade-off Sv vs. massa
f_sweep = 10.0 .^ range(log10(10.0), log10(500.0), length=80)
p3_left = plot(
    xlabel = "Induction Frequency (kHz)",
    ylabel = "Macroscopic Coupling Factor, F_cyl (-)",
    title = "Inductive Coupling vs. Frequency (800°C)",
    xscale = :log10, ylims = (0.0, 1.02),
    legend = :bottomright,
)
for g in geometries
    F_vals = [analyze_monolith_coupling(g, f * 1000.0, 800.0).F_macro_cyl for f in f_sweep]
    plot!(p3_left, f_sweep, F_vals, label = g.name)
end
hline!(p3_left, [0.65], label = "Strong coupling knee (D/delta_eff ~ 3.5)",
    color = :gray, linestyle = :dot)

p3_right = scatter(
    df_geom.mass_g, df_geom.Sv_m2_m3,
    xlabel = "Inconel 600 Mass (g)",
    ylabel = "Specific Catalytic Area, Sv (m2/m3)",
    title = "Catalytic Area vs. Metal Mass (2x3 cm)",
    markersize = 8, color = :steelblue,
    legend = false,
)
for r in eachrow(df_geom)
    annotate!(p3_right, r.mass_g + 0.4, r.Sv_m2_m3, text(r.geometry, 7, :left))
end
xlims!(p3_right, 0.0, maximum(df_geom.mass_g) * 1.45)

fig3 = plot(p3_left, p3_right, layout = (1, 2), size = (1420, 500))
fig3_path = joinpath(OUTPUT_DIR, "fig_p2_3_geometry_comparison.png")
savefig(fig3, fig3_path)
println("\n      Salvo: ", csv_path, " e ", fig3_path)

elapsed = round(time() - T_START, digits=2)
println("Concluido em ", elapsed, " s!")