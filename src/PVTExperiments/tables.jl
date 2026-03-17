"""
    pvto_table(eos, z, T; p_range, n_rs, n_undersaturated, p_sc, T_sc, separator_stages)

Generate a PVTO (live oil) table from DLE or MSS experiment data.

# Arguments
- `eos`: Equation of state
- `z`: Overall composition (mole fractions)
- `T`: Reservoir temperature (K)
- `p_range`: Pressure range (p_max, p_min). Default: (50e6, 1e5)
- `n_rs`: Number of Rs (dissolved gas) levels. Default: 15
- `n_undersaturated`: Number of undersaturated pressure steps per Rs level. Default: 5
- `p_sc`: Standard condition pressure (Pa). Default: 101325.0
- `T_sc`: Standard condition temperature (K). Default: 288.706
- `separator_stages`: Optional separator stages. If provided, MSS is used for surface conditions.
"""
function pvto_table(eos, z, T;
        p_range = (50e6, 1e5),
        n_rs = 15,
        n_undersaturated = 5,
        p_sc = 101325.0,
        T_sc = 288.706,
        separator_stages = nothing
    )
    z = collect(Float64, z)
    z ./= sum(z)

    # Find saturation pressure
    p_sat, _ = find_saturation_pressure(eos, z, T; p_min = p_range[2], p_max = p_range[1])

    # Set up surface flash conditions
    if isnothing(separator_stages)
        surface_stages = [SeparatorStage(p_sc, T_sc)]
    else
        surface_stages = separator_stages
    end

    # DLE-like liberation: step from p_sat to p_min
    pressures_lib = collect(range(p_sat, p_range[2], length = n_rs + 1))

    # Storage for the table
    Rs_values = Float64[]
    p_bub_values = Float64[]
    p_table = Vector{Vector{Float64}}()
    Bo_table = Vector{Vector{Float64}}()
    mu_o_table = Vector{Vector{Float64}}()

    z_oil = copy(z)
    n_oil = 1.0
    V_gas_sc_per_mol = phase_molar_volume_from_Z(p_sc, T_sc, 1.0)

    # Track liberation data
    lib_data = []

    for (i, p_lib) in enumerate(pressures_lib)
        if i == 1
            # At bubble point
            props = flash_and_properties(eos, p_lib, T, z_oil)
            push!(lib_data, (p = p_lib, z_oil = copy(z_oil), n_oil = n_oil, props = props))
        else
            # Flash and liberate gas
            props = flash_and_properties(eos, p_lib, T, z_oil)
            n_gas = n_oil * props.V
            n_oil = n_oil * (1.0 - props.V)
            z_oil = copy(props.x)
            push!(lib_data, (p = p_lib, z_oil = copy(z_oil), n_oil = n_oil, props = props))
        end
    end

    # Compute stock tank oil reference from final liquid
    sep_res = separator_flash(eos, z_oil, surface_stages)
    z_st_oil = sep_res.oil_compositions[end]
    props_st = flash_and_properties(eos, p_sc, T_sc, z_st_oil)
    V_mol_st_oil = props_st.V_mol_l
    n_oil_at_st = n_oil
    for V_i in sep_res.V_stages
        n_oil_at_st *= (1.0 - V_i)
    end

    # Build the table in reverse (from dead oil up to live oil)
    for (idx, ld) in enumerate(reverse(lib_data))
        p_bubble = ld.p
        z_oil_at_level = ld.z_oil

        # Flash this oil through separators to get surface products
        sep = separator_flash(eos, z_oil_at_level, surface_stages)
        z_st = sep.oil_compositions[end]
        props_st_level = flash_and_properties(eos, p_sc, T_sc, z_st)

        # Liquid fraction through separator
        liq_frac = 1.0
        gas_total_sc = 0.0
        for (j, V_j) in enumerate(sep.V_stages)
            gas_at_j = liq_frac * V_j
            gas_total_sc += gas_at_j * V_gas_sc_per_mol
            liq_frac *= (1.0 - V_j)
        end
        V_oil_st = liq_frac * props_st_level.V_mol_l

        if V_oil_st ≤ 0.0
            continue
        end

        Rs = max(0.0, gas_total_sc / V_oil_st)

        # Properties at bubble point
        props_bub = flash_and_properties(eos, p_bubble, T, z_oil_at_level)
        Bo_bub = props_bub.V_mol_l / (V_oil_st / liq_frac)

        # Saturated entry
        p_entries = Float64[p_bubble]
        Bo_entries = Float64[Bo_bub]
        mu_entries = Float64[props_bub.μ_l]

        # Undersaturated entries: higher pressures
        if p_bubble < p_range[1]
            p_max_us = min(p_range[1], 2.0 * p_bubble)
            dp = (p_max_us - p_bubble) / n_undersaturated
            for k in 1:n_undersaturated
                p_us = p_bubble + k * dp
                props_us = flash_and_properties(eos, p_us, T, z_oil_at_level)
                Bo_us = props_us.V_mol_l / (V_oil_st / liq_frac)
                push!(p_entries, p_us)
                push!(Bo_entries, Bo_us)
                push!(mu_entries, props_us.μ_l)
            end
        end

        push!(Rs_values, Rs)
        push!(p_bub_values, p_bubble)
        push!(p_table, p_entries)
        push!(Bo_table, Bo_entries)
        push!(mu_o_table, mu_entries)
    end

    # Sort by Rs
    order = sortperm(Rs_values)
    Rs_values = Rs_values[order]
    p_bub_values = p_bub_values[order]
    p_table = p_table[order]
    Bo_table = Bo_table[order]
    mu_o_table = mu_o_table[order]

    return PVTOTable(Rs_values, p_bub_values, p_table, Bo_table, mu_o_table)
end

"""
    pvdg_table(eos, z, T; p_range, n_points, p_sc, T_sc)

Generate a PVDG (dry gas) table.

# Arguments
- `eos`: Equation of state
- `z`: Gas composition (mole fractions)
- `T`: Temperature (K)
- `p_range`: Pressure range. Default: (50e6, 1e5)
- `n_points`: Number of pressure points. Default: 20
- `p_sc`: Standard condition pressure (Pa). Default: 101325.0
- `T_sc`: Standard condition temperature (K). Default: 288.706
"""
function pvdg_table(eos, z, T;
        p_range = (50e6, 1e5),
        n_points = 20,
        p_sc = 101325.0,
        T_sc = 288.706
    )
    z = collect(Float64, z)
    z ./= sum(z)

    pressures = collect(range(p_range[2], p_range[1], length = n_points))
    sort!(pressures)

    Bg_arr = zeros(n_points)
    mu_g_arr = zeros(n_points)

    V_mol_sc = phase_molar_volume_from_Z(p_sc, T_sc, 1.0)

    for (i, p) in enumerate(pressures)
        props = flash_and_properties(eos, p, T, z)
        V_mol_gas = props.V_mol_v
        Bg_arr[i] = V_mol_gas / V_mol_sc
        mu_g_arr[i] = props.μ_v
    end

    return PVDGTable(pressures, Bg_arr, mu_g_arr)
end

"""
    pvtg_table(eos, z, T; p_range, n_rv, n_undersaturated, p_sc, T_sc)

Generate a PVTG (wet gas / gas condensate) table.

# Arguments
- `eos`: Equation of state
- `z`: Overall composition (mole fractions)
- `T`: Temperature (K)
- `p_range`: Pressure range. Default: (50e6, 1e5)
- `n_rv`: Number of pressure (Rv) levels. Default: 15
- `n_undersaturated`: Number of undersaturated entries per level. Default: 3
- `p_sc`: Standard condition pressure (Pa). Default: 101325.0
- `T_sc`: Standard condition temperature (K). Default: 288.706
"""
function pvtg_table(eos, z, T;
        p_range = (50e6, 1e5),
        n_rv = 15,
        n_undersaturated = 3,
        p_sc = 101325.0,
        T_sc = 288.706
    )
    z = collect(Float64, z)
    z ./= sum(z)

    # Find dew point
    p_sat, _ = find_saturation_pressure(eos, z, T; p_min = p_range[2], p_max = p_range[1])

    # Pressure steps
    pressures = collect(range(p_sat, p_range[2], length = n_rv + 1))

    V_mol_sc = phase_molar_volume_from_Z(p_sc, T_sc, 1.0)

    p_values = Float64[]
    Rv_sat_values = Float64[]
    Rv_sub_table = Vector{Vector{Float64}}()
    Bg_table = Vector{Vector{Float64}}()
    mu_g_table = Vector{Vector{Float64}}()

    # CVD-like approach for gas condensate
    z_gas = copy(z)

    for (i, p) in enumerate(pressures)
        props = flash_and_properties(eos, p, T, z_gas)

        # Saturated gas properties
        V_mol_gas = props.V_mol_v
        Bg_sat = V_mol_gas / V_mol_sc

        # Rv: liquid content in gas at surface conditions
        # Flash the vapor phase at surface conditions
        if props.state == MultiComponentFlash.two_phase_lv
            y = props.y
        else
            y = copy(z_gas)
        end

        props_gas_sc = flash_and_properties(eos, p_sc, T_sc, y)
        # Rv = volume of oil from gas at SC / volume of gas at SC
        if props_gas_sc.state == MultiComponentFlash.two_phase_lv
            V_oil_from_gas = (1.0 - props_gas_sc.V) * props_gas_sc.V_mol_l
            V_gas_at_sc = props_gas_sc.V * props_gas_sc.V_mol_v
            if V_gas_at_sc > 0
                Rv = V_oil_from_gas / V_gas_at_sc
            else
                Rv = 0.0
            end
        else
            Rv = 0.0
        end

        # Saturated entry
        rv_entries = Float64[Rv]
        bg_entries = Float64[Bg_sat]
        mu_entries = Float64[props.μ_v]

        # Undersaturated entries: reduce Rv (leaner gas)
        if Rv > 0
            for k in 1:n_undersaturated
                frac = 1.0 - k / (n_undersaturated + 1)
                # Create a leaner gas by mixing with pure gas component
                z_lean = frac .* y .+ (1.0 - frac) .* z_gas
                z_lean ./= sum(z_lean)
                props_lean = flash_and_properties(eos, p, T, z_lean)
                Bg_lean = props_lean.V_mol_v / V_mol_sc

                # Get Rv for lean gas
                props_lean_sc = flash_and_properties(eos, p_sc, T_sc, z_lean)
                if props_lean_sc.state == MultiComponentFlash.two_phase_lv
                    V_oil_l = (1.0 - props_lean_sc.V) * props_lean_sc.V_mol_l
                    V_gas_l = props_lean_sc.V * props_lean_sc.V_mol_v
                    Rv_lean = V_gas_l > 0 ? V_oil_l / V_gas_l : 0.0
                else
                    Rv_lean = 0.0
                end

                push!(rv_entries, Rv_lean)
                push!(bg_entries, Bg_lean)
                push!(mu_entries, props_lean.μ_v)
            end
        end

        # Add zero Rv entry
        pushfirst!(rv_entries, 0.0)
        # Dry gas Bg and viscosity
        # Use the lightest component (index 1) as a proxy for dry gas.
        # This assumes the mixture is ordered with the lightest component first.
        z_first = zeros(length(z))
        z_first[1] = 1.0
        props_dry = flash_and_properties(eos, p, T, z_first)
        pushfirst!(bg_entries, props_dry.V_mol_v / V_mol_sc)
        pushfirst!(mu_entries, props_dry.μ_v)

        push!(p_values, p)
        push!(Rv_sat_values, Rv)
        push!(Rv_sub_table, rv_entries)
        push!(Bg_table, bg_entries)
        push!(mu_g_table, mu_entries)

        # Update gas composition (CVD-like)
        if props.state == MultiComponentFlash.two_phase_lv
            z_gas = copy(props.y)
        end
    end

    return PVTGTable(p_values, Rv_sat_values, Rv_sub_table, Bg_table, mu_g_table)
end

"""
    pvdo_table(eos, z, T; p_range, n_points, p_sc, T_sc)

Generate a PVDO (dead oil) table.

# Arguments
- `eos`: Equation of state
- `z`: Oil composition (mole fractions)
- `T`: Temperature (K)
- `p_range`: Pressure range. Default: (50e6, 1e5)
- `n_points`: Number of pressure points. Default: 20
- `p_sc`: Standard condition pressure (Pa). Default: 101325.0
- `T_sc`: Standard condition temperature (K). Default: 288.706
"""
function pvdo_table(eos, z, T;
        p_range = (50e6, 1e5),
        n_points = 20,
        p_sc = 101325.0,
        T_sc = 288.706
    )
    z = collect(Float64, z)
    z ./= sum(z)

    pressures = collect(range(p_range[2], p_range[1], length = n_points))
    sort!(pressures)

    Bo_arr = zeros(n_points)
    mu_o_arr = zeros(n_points)

    props_sc = flash_and_properties(eos, p_sc, T_sc, z)
    V_mol_sc = props_sc.V_mol_l

    for (i, p) in enumerate(pressures)
        props = flash_and_properties(eos, p, T, z)
        Bo_arr[i] = props.V_mol_l / V_mol_sc
        mu_o_arr[i] = props.μ_l
    end

    return PVDOTable(pressures, Bo_arr, mu_o_arr)
end

"""
    surface_densities(eos, z, T; p_sc, T_sc, separator_stages)

Compute surface densities of oil and gas at standard conditions.

# Arguments
- `eos`: Equation of state
- `z`: Overall composition (mole fractions)
- `T`: Reservoir temperature (K)
- `p_sc`: Standard condition pressure (Pa). Default: 101325.0
- `T_sc`: Standard condition temperature (K). Default: 288.706
- `separator_stages`: Optional separator stages for surface conditions.
"""
function surface_densities(eos, z, T;
        p_sc = 101325.0,
        T_sc = 288.706,
        separator_stages = nothing
    )
    z = collect(Float64, z)
    z ./= sum(z)

    if isnothing(separator_stages)
        surface_stages = [SeparatorStage(p_sc, T_sc)]
    else
        surface_stages = separator_stages
    end

    # Flash through separator stages
    sep = separator_flash(eos, z, surface_stages)

    # Oil density at stock tank
    z_oil_st = sep.oil_compositions[end]
    props_oil = flash_and_properties(eos, p_sc, T_sc, z_oil_st)
    ρ_oil = props_oil.ρ_l

    # Gas density at stock tank
    z_gas_st = sep.gas_compositions[end]
    mw_gas = mixture_molar_mass(z_gas_st, eos)
    V_mol_gas = phase_molar_volume_from_Z(p_sc, T_sc, 1.0)
    ρ_gas = mw_gas / V_mol_gas

    return SurfaceDensities(ρ_oil, ρ_gas)
end
