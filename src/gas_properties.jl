"""
    ThermoPolynomial

Single temperature-range polynomial loaded from a Cantera YAML file.
Supports both `:NASA7` and `:NASA9` models.
"""
struct ThermoPolynomial
    model::Symbol
    T_min_K::Float64
    T_max_K::Float64
    coeffs::Vector{Float64}
end

"""
    SpeciesData

Thermodynamic and atomic record parsed directly from official Cantera/NASA YAML files.
"""
struct SpeciesData
    symbol::Symbol
    yaml_name::String
    phase::Symbol  # :gas or :solid
    atoms::NamedTuple{(:C, :H, :O, :N), NTuple{4, Int}}
    ranges::Vector{ThermoPolynomial}
    source_file::String
end

const BASE_SPECIES = [:CH4, :CO2, :H2, :CO, :H2O, :N2, :O2, :C_s]
const EXTENDED_SPECIES = vcat(BASE_SPECIES, [:C2H6, :C2H4, :CH3OH, :CH3OCH3])

# Maps internal symbols to (phase, candidate names in Cantera YAML files)
const _SPECIES_MAPPING = Dict{Symbol, Tuple{Symbol, Vector{String}}}(
    :CH4     => (:gas,   ["CH4"]),
    :CO2     => (:gas,   ["CO2"]),
    :CO      => (:gas,   ["CO"]),
    :H2      => (:gas,   ["H2"]),
    :H2O     => (:gas,   ["H2O"]),
    :N2      => (:gas,   ["N2"]),
    :O2      => (:gas,   ["O2"]),
    :C_s     => (:solid, ["C(gr)"]),
    :C2H6    => (:gas,   ["C2H6"]),
    :C2H4    => (:gas,   ["C2H4"]),
    :CH3OH   => (:gas,   ["CH3OH"]),
    :CH3OCH3 => (:gas,   ["CH3OCH3", "CH3-O-CH3"]),
)

const _DEFAULT_DATA_DIR = normpath(joinpath(@__DIR__, "..", "data"))
const _CACHED_DB = Ref{Union{Nothing, Dict{Symbol, SpeciesData}}}(nothing)

function _parse_species_entry(
    sym::Symbol,
    entry::AbstractDict,
    phase::Symbol,
    source_file::String,
)::SpeciesData
    comp = entry["composition"]
    allowed = Set(["C", "H", "O", "N"])
    for elem in keys(comp)
        String(elem) in allowed || throw(ArgumentError("Unsupported element $elem in species $sym"))
    end
    atoms = (
        C = Int(get(comp, "C", 0)),
        H = Int(get(comp, "H", 0)),
        O = Int(get(comp, "O", 0)),
        N = Int(get(comp, "N", 0)),
    )

    thermo = entry["thermo"]
    model_str = uppercase(strip(String(thermo["model"])))
    model = model_str == "NASA9" ? :NASA9 :
            model_str == "NASA7" ? :NASA7 :
            throw(ArgumentError("Unsupported thermo model '$model_str' for $sym"))

    t_bounds = Float64.(thermo["temperature-ranges"])
    coeff_blocks = thermo["data"]
    length(t_bounds) == length(coeff_blocks) + 1 ||
        throw(ArgumentError("Mismatched temperature-ranges and data blocks for $sym"))

    expected_len = model === :NASA9 ? 9 : 7
    ranges = ThermoPolynomial[]
    for i in eachindex(coeff_blocks)
        c = Float64.(coeff_blocks[i])
        length(c) == expected_len ||
            throw(ArgumentError("Expected $expected_len coefficients for $model ($sym), got $(length(c))"))
        push!(ranges, ThermoPolynomial(model, t_bounds[i], t_bounds[i + 1], c))
    end

    return SpeciesData(sym, String(entry["name"]), phase, atoms, ranges, basename(source_file))
end

"""
    load_thermo_database(data_dir::AbstractString = _DEFAULT_DATA_DIR) -> Dict{Symbol, SpeciesData}

Load official NASA/Cantera YAML files (`nasa_gas.yaml` and `nasa_condensed.yaml`)
from `data_dir` and index all required species.
"""
function load_thermo_database(data_dir::AbstractString=_DEFAULT_DATA_DIR)::Dict{Symbol, SpeciesData}
    gas_file = joinpath(data_dir, "nasa_gas.yaml")
    solid_file = joinpath(data_dir, "nasa_condensed.yaml")

    isfile(gas_file) || error("Missing official thermodynamic file: $gas_file")
    isfile(solid_file) || error("Missing official thermodynamic file: $solid_file")

    gas_raw = YAML.load_file(gas_file)
    solid_raw = YAML.load_file(solid_file)

    gas_by_name = Dict(String(s["name"]) => s for s in gas_raw["species"])
    solid_by_name = Dict(String(s["name"]) => s for s in solid_raw["species"])

    db = Dict{Symbol, SpeciesData}()
    for (sym, (phase, candidates)) in _SPECIES_MAPPING
        lookup = phase === :gas ? gas_by_name : solid_by_name
        src_file = phase === :gas ? gas_file : solid_file

        found_entry = nothing
        for name in candidates
            if haskey(lookup, name)
                found_entry = lookup[name]
                break
            end
        end

        if found_entry === nothing
            if sym in BASE_SPECIES
                error("Required base species $sym (candidates: $candidates) not found in $src_file")
            end
            continue
        end

        db[sym] = _parse_species_entry(sym, found_entry, phase, src_file)
    end

    return db
end

"""
    default_thermo_db() -> Dict{Symbol, SpeciesData}

Return the cached thermodynamic database loaded from `data/`.
"""
function default_thermo_db()::Dict{Symbol, SpeciesData}
    if _CACHED_DB[] === nothing
        _CACHED_DB[] = load_thermo_database()
    end
    return _CACHED_DB[]::Dict{Symbol, SpeciesData}
end

"""
    get_species_data(species::Symbol; db = default_thermo_db()) -> SpeciesData
"""
function get_species_data(
    species::Symbol;
    db::Dict{Symbol, SpeciesData}=default_thermo_db(),
)::SpeciesData
    haskey(db, species) || throw(KeyError(species))
    return db[species]
end

"""
    species_list(mode::Symbol = :base; db = default_thermo_db()) -> Vector{Symbol}

Return the active chemical species for `:base` (8 species) or `:extended` (12 species).
"""
function species_list(
    mode::Symbol=:base;
    db::Dict{Symbol, SpeciesData}=default_thermo_db(),
)::Vector{Symbol}
    if mode === :base
        return copy(BASE_SPECIES)
    elseif mode === :extended
        missing_sp = [sp for sp in EXTENDED_SPECIES if !haskey(db, sp)]
        isempty(missing_sp) || error("Missing extended species in YAML database: $missing_sp")
        return copy(EXTENDED_SPECIES)
    else
        throw(ArgumentError("Unsupported species mode: $mode. Use :base or :extended."))
    end
end

function _select_poly(data::SpeciesData, T_K::Float64)::ThermoPolynomial
    for p in data.ranges
        if p.T_min_K <= T_K <= p.T_max_K
            return p
        end
    end
    t_lo = first(data.ranges).T_min_K
    t_hi = last(data.ranges).T_max_K
    throw(DomainError(T_K, "Temperature $T_K K outside valid range [$t_lo, $t_hi] K for $(data.symbol)."))
end

function _eval_cp_over_R(p::ThermoPolynomial, T::Float64)::Float64
    c = p.coeffs
    if p.model === :NASA9
        # NASA/TP-2002-211556 Eq. (1)
        return c[1] / T^2 + c[2] / T + c[3] + c[4] * T + c[5] * T^2 + c[6] * T^3 + c[7] * T^4
    else
        # NASA-7
        return c[1] + c[2] * T + c[3] * T^2 + c[4] * T^3 + c[5] * T^4
    end
end

function _eval_h_over_RT(p::ThermoPolynomial, T::Float64)::Float64
    c = p.coeffs
    if p.model === :NASA9
        # NASA/TP-2002-211556 Eq. (2)
        return -c[1] / T^2 + c[2] * log(T) / T + c[3] +
               (c[4] * T) / 2.0 + (c[5] * T^2) / 3.0 + (c[6] * T^3) / 4.0 + (c[7] * T^4) / 5.0 +
               c[8] / T
    else
        # NASA-7
        return c[1] + (c[2] * T) / 2.0 + (c[3] * T^2) / 3.0 +
               (c[4] * T^3) / 4.0 + (c[5] * T^4) / 5.0 + c[6] / T
    end
end

function _eval_s_over_R(p::ThermoPolynomial, T::Float64)::Float64
    c = p.coeffs
    if p.model === :NASA9
        # NASA/TP-2002-211556 Eq. (3)
        return -c[1] / (2.0 * T^2) - c[2] / T + c[3] * log(T) +
               c[4] * T + (c[5] * T^2) / 2.0 + (c[6] * T^3) / 3.0 + (c[7] * T^4) / 4.0 +
               c[9]
    else
        # NASA-7
        return c[1] * log(T) + c[2] * T + (c[3] * T^2) / 2.0 +
               (c[4] * T^3) / 3.0 + (c[5] * T^4) / 4.0 + c[7]
    end
end

"""
    heat_capacity_J_mol_K(species::Symbol, T_K::Real; db = default_thermo_db()) -> Float64

Standard molar heat capacity Cp°(T) in J/(mol·K).
"""
function heat_capacity_J_mol_K(
    species::Symbol,
    T_K::Real;
    db::Dict{Symbol, SpeciesData}=default_thermo_db(),
)::Float64
    T = Float64(T_K)
    p = _select_poly(get_species_data(species; db=db), T)
    return R_GAS * _eval_cp_over_R(p, T)
end

"""
    enthalpy_J_mol(species::Symbol, T_K::Real; db = default_thermo_db()) -> Float64

Total standard molar enthalpy H°(T) in J/mol (includes ΔfH° at 298.15 K).
"""
function enthalpy_J_mol(
    species::Symbol,
    T_K::Real;
    db::Dict{Symbol, SpeciesData}=default_thermo_db(),
)::Float64
    T = Float64(T_K)
    p = _select_poly(get_species_data(species; db=db), T)
    return R_GAS * T * _eval_h_over_RT(p, T)
end

"""
    sensible_enthalpy_J_mol(species::Symbol, T_K::Real; T_ref_K = 298.15, db = default_thermo_db()) -> Float64

Sensible enthalpy change H°(T) - H°(T_ref_K) in J/mol.
"""
function sensible_enthalpy_J_mol(
    species::Symbol,
    T_K::Real;
    T_ref_K::Real=298.15,
    db::Dict{Symbol, SpeciesData}=default_thermo_db(),
)::Float64
    return enthalpy_J_mol(species, T_K; db=db) - enthalpy_J_mol(species, T_ref_K; db=db)
end

"""
    entropy_J_mol_K(species::Symbol, T_K::Real; db = default_thermo_db()) -> Float64

Standard molar absolute entropy S°(T) in J/(mol·K) at standard pressure (1 bar for NASA-9).
"""
function entropy_J_mol_K(
    species::Symbol,
    T_K::Real;
    db::Dict{Symbol, SpeciesData}=default_thermo_db(),
)::Float64
    T = Float64(T_K)
    p = _select_poly(get_species_data(species; db=db), T)
    return R_GAS * _eval_s_over_R(p, T)
end

"""
    gibbs_energy_J_mol(species::Symbol, T_K::Real; db = default_thermo_db()) -> Float64

Standard molar Gibbs free energy G°(T) = H°(T) - T * S°(T) in J/mol.
"""
function gibbs_energy_J_mol(
    species::Symbol,
    T_K::Real;
    db::Dict{Symbol, SpeciesData}=default_thermo_db(),
)::Float64
    T = Float64(T_K)
    return enthalpy_J_mol(species, T; db=db) - T * entropy_J_mol_K(species, T; db=db)
end