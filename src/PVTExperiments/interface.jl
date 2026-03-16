"""
    generate_pvt_tables(eos, z, T_res; <keyword arguments>)

High-level interface that generates complete black oil PVT tables from a
compositional fluid description.

Goes from a fluid sample, reservoir temperature and surface conditions to
PVTO/PVDG (or PVTG/PVDG) tables + surface densities.

# Arguments
- `eos`: Equation of state (e.g., `GenericCubicEOS(mixture, PengRobinson())`)
- `z`: Overall mole fractions of reservoir fluid
- `T_res`: Reservoir temperature (K)

# Keyword Arguments
- `disgas::Union{Bool,Nothing} = nothing`: Allow dissolved gas in oil phase.
  If `nothing`, auto-detected from fluid type.
- `vapoil::Union{Bool,Nothing} = nothing`: Allow vaporized oil in gas phase.
  If `nothing`, auto-detected from fluid type.
- `p_range`: Pressure range for tables. Default: (50e6, 1e5)
- `p_sc`: Standard condition pressure (Pa). Default: 101325.0
- `T_sc`: Standard condition temperature (K). Default: 288.706
- `separator_stages`: Separator stage definitions. If provided, MSS is used
  for surface reference; otherwise DLE is used.
- `n_pvto`: Number of Rs levels for PVTO. Default: 15
- `n_pvtg`: Number of pressure levels for PVTG. Default: 15
- `n_pvdg`: Number of pressure points for PVDG. Default: 20
- `n_undersaturated`: Undersaturated points per level. Default: 5

# Returns
A NamedTuple with fields:
- `pvto`: PVTOTable or `nothing`
- `pvtg`: PVTGTable or `nothing`
- `pvdg`: PVDGTable
- `surface_densities`: SurfaceDensities
- `saturation_pressure`: Saturation pressure (Pa)
- `is_bubblepoint`: Whether the fluid has a bubble point (oil) or dew point (gas)
"""
function generate_pvt_tables(eos, z, T_res;
        disgas::Union{Bool,Nothing} = nothing,
        vapoil::Union{Bool,Nothing} = nothing,
        p_range = (50e6, 1e5),
        p_sc = 101325.0,
        T_sc = 288.706,
        separator_stages = nothing,
        n_pvto = 15,
        n_pvtg = 15,
        n_pvdg = 20,
        n_undersaturated = 5
    )
    z = collect(Float64, z)
    z ./= sum(z)

    # Find saturation pressure and fluid type
    p_sat, is_bp = find_saturation_pressure(eos, z, T_res;
        p_min = p_range[2], p_max = p_range[1])

    # Auto-detect table types if not specified
    if isnothing(disgas)
        disgas = is_bp  # Oil system → dissolved gas
    end
    if isnothing(vapoil)
        vapoil = !is_bp  # Gas system → vaporized oil
    end

    # Surface stages
    if isnothing(separator_stages)
        surf_stages = [SeparatorStage(p_sc, T_sc)]
    else
        surf_stages = separator_stages
    end

    # Generate tables
    pvto_result = nothing
    pvtg_result = nothing

    if disgas
        pvto_result = pvto_table(eos, z, T_res;
            p_range = p_range,
            n_rs = n_pvto,
            n_undersaturated = n_undersaturated,
            p_sc = p_sc,
            T_sc = T_sc,
            separator_stages = isnothing(separator_stages) ? nothing : separator_stages
        )
    end

    if vapoil
        pvtg_result = pvtg_table(eos, z, T_res;
            p_range = p_range,
            n_rv = n_pvtg,
            n_undersaturated = 3,
            p_sc = p_sc,
            T_sc = T_sc
        )
    end

    # PVDG is always generated
    # For oil systems, use the gas composition from the first liberation step
    # For gas systems, use the overall composition
    if is_bp
        # Get gas composition at bubble point
        props_bp = flash_and_properties(eos, p_sat * 0.95, T_res, z)
        z_gas = props_bp.y
    else
        z_gas = z
    end

    pvdg_result = pvdg_table(eos, z_gas, T_res;
        p_range = p_range,
        n_points = n_pvdg,
        p_sc = p_sc,
        T_sc = T_sc
    )

    # Surface densities
    sd = surface_densities(eos, z, T_res;
        p_sc = p_sc,
        T_sc = T_sc,
        separator_stages = separator_stages
    )

    return (
        pvto = pvto_result,
        pvtg = pvtg_result,
        pvdg = pvdg_result,
        surface_densities = sd,
        saturation_pressure = p_sat,
        is_bubblepoint = is_bp
    )
end
