"""
    find_saturation_pressure(eos, z, T; p_min, p_max, n_points)

Find the saturation pressure (bubble or dew point) for a given composition
and temperature using bisection.

Returns `(p_sat, is_bubblepoint)`.
"""
function find_saturation_pressure(eos, z, T;
        p_min = 1e5,
        p_max = 100e6,
        n_points = 50,
        tol = 1e4
    )
    # Scan from high to low pressure to find where the phase transition occurs
    pressures = range(p_max, p_min, length = n_points)
    z_work = copy(z)
    method = SSIFlash()
    S = flash_storage(eos, (p = p_max, T = T, z = z_work), method = method)
    K = initial_guess_K(eos, (p = p_max, T = T, z = z_work))

    # Find the first pressure where two phases exist
    p_upper = p_max
    p_lower = p_min
    found_two_phase = false
    found_single = false

    for p in pressures
        cond = (p = p, T = T, z = z_work)
        K .= 1.0
        initial_guess_K!(K, eos, cond)
        V, K_out, rep = flash_2ph!(S, K, eos, cond; extra_out = true, check = false)
        if isnan(V)
            if found_two_phase
                p_lower = p
                break
            end
            found_single = true
            p_upper = p
        else
            if !found_two_phase && found_single
                p_upper = pressures[max(1, findfirst(==(p), pressures) - 1)]
            end
            found_two_phase = true
            p_lower = p
        end
    end

    if !found_two_phase
        # Mixture is single-phase everywhere in the range
        # Use label to determine phase type
        label = single_phase_label(eos, (p = p_max, T = T, z = z_work))
        return (p_min, label < 0.5)
    end

    # Bisection to refine
    for _ in 1:60  # max bisection iterations
        p_mid = 0.5 * (p_upper + p_lower)
        cond = (p = p_mid, T = T, z = z_work)
        K .= 1.0
        initial_guess_K!(K, eos, cond)
        V, _, _ = flash_2ph!(S, K, eos, cond; extra_out = true, check = false)
        if isnan(V)
            p_upper = p_mid
        else
            p_lower = p_mid
        end
        if abs(p_upper - p_lower) < tol
            break
        end
    end

    p_sat = 0.5 * (p_upper + p_lower)

    # Determine if bubble point or dew point
    # Check composition: flash just below p_sat
    cond = (p = p_sat * 0.99, T = T, z = z_work)
    K .= 1.0
    initial_guess_K!(K, eos, cond)
    V, K_out, _ = flash_2ph!(S, K, eos, cond; extra_out = true, check = false)
    if isnan(V)
        # Barely two-phase, use label
        label = single_phase_label(eos, (p = p_sat, T = T, z = z_work))
        is_bubblepoint = label < 0.5
    else
        is_bubblepoint = V < 0.5
    end

    return (p_sat, is_bubblepoint)
end

"""
    flash_and_properties(eos, p, T, z)

Perform a flash calculation and extract all relevant properties.
Returns a NamedTuple with phase properties or `nothing` if flash fails.
"""
function flash_and_properties(eos, p, T, z)
    z_work = copy(z)
    result = flashed_mixture_2ph(eos, (p = p, T = T, z = z_work); method = SSIFlash())
    state = result.state
    V = result.V

    if state == MultiComponentFlash.two_phase_lv
        x = copy(result.liquid.mole_fractions)
        y = copy(result.vapor.mole_fractions)
        Z_L = result.liquid.Z
        Z_V = result.vapor.Z
        ρ_l = mass_density(eos, p, T, result.liquid)
        ρ_v = mass_density(eos, p, T, result.vapor)
        μ_l = lbc_viscosity(eos, p, T, result.liquid)
        μ_v = lbc_viscosity(eos, p, T, result.vapor)
        V_mol_l = molar_volume(eos, p, T, result.liquid)
        V_mol_v = molar_volume(eos, p, T, result.vapor)
    elseif state == MultiComponentFlash.single_phase_l
        x = copy(z_work)
        y = copy(z_work)
        ph = result.liquid
        Z_L = ph.Z
        Z_V = ph.Z
        ρ_l = mass_density(eos, p, T, ph)
        ρ_v = ρ_l
        μ_l = lbc_viscosity(eos, p, T, ph)
        μ_v = μ_l
        V_mol_l = molar_volume(eos, p, T, ph)
        V_mol_v = V_mol_l
        V = 0.0
    else
        x = copy(z_work)
        y = copy(z_work)
        ph = result.vapor
        Z_L = ph.Z
        Z_V = ph.Z
        ρ_l = mass_density(eos, p, T, ph)
        ρ_v = ρ_l
        μ_l = lbc_viscosity(eos, p, T, ph)
        μ_v = μ_l
        V_mol_l = molar_volume(eos, p, T, ph)
        V_mol_v = V_mol_l
        V = 1.0
    end

    return (
        state = state,
        V = V,
        x = x,
        y = y,
        Z_L = Z_L,
        Z_V = Z_V,
        ρ_l = ρ_l,
        ρ_v = ρ_v,
        μ_l = μ_l,
        μ_v = μ_v,
        V_mol_l = V_mol_l,
        V_mol_v = V_mol_v
    )
end

"""
    mixture_molar_mass(z, eos)

Compute the molar mass of a mixture given composition z.
"""
function mixture_molar_mass(z, eos)
    mw = 0.0
    props = eos.mixture.properties
    for (i, prop) in enumerate(props)
        mw += z[i] * prop.mw
    end
    return mw
end

"""
    phase_molar_volume_from_Z(p, T, Z)

Compute molar volume from compressibility factor: V = ZRT/p
"""
function phase_molar_volume_from_Z(p, T, Z)
    return Z * R_GAS * T / p
end

"""
    separator_flash(eos, z_feed, stages::Vector{SeparatorStage})

Run a multistage separator train. Returns the compositions and properties at each stage.
"""
function separator_flash(eos, z_feed, stages::Vector{SeparatorStage})
    n_stages = length(stages)
    nc = number_of_components(eos)

    gas_compositions = Vector{Vector{Float64}}(undef, n_stages)
    oil_compositions = Vector{Vector{Float64}}(undef, n_stages)
    V_stages = zeros(n_stages)

    z_current = copy(z_feed)

    for (i, stage) in enumerate(stages)
        props = flash_and_properties(eos, stage.p, stage.T, z_current)
        if props.state == MultiComponentFlash.two_phase_lv
            gas_compositions[i] = copy(props.y)
            oil_compositions[i] = copy(props.x)
            V_stages[i] = props.V
            z_current = copy(props.x)  # liquid goes to next stage
        else
            # Single phase - no separation
            gas_compositions[i] = copy(z_current)
            oil_compositions[i] = copy(z_current)
            V_stages[i] = props.V
            # z_current stays the same
        end
    end

    return (
        gas_compositions = gas_compositions,
        oil_compositions = oil_compositions,
        V_stages = V_stages
    )
end

# Print functions
function Base.show(io::IO, r::CCEResult)
    println(io, "CCE Experiment Results")
    println(io, "=" ^ 60)
    @printf(io, "Temperature: %.2f K (%.2f °C)\n", r.T, r.T - 273.15)
    @printf(io, "Saturation pressure: %.4e Pa (%.2f bar)\n", r.p_sat, r.p_sat / 1e5)
    println(io, r.is_bubblepoint ? "Type: Bubble point" : "Type: Dew point")
    println(io, "-" ^ 60)
    @printf(io, "%14s %10s %10s %10s %12s\n", "P (Pa)", "V_rel", "V_frac", "ρ_oil", "ρ_gas")
    @printf(io, "%14s %10s %10s %10s %12s\n", "", "", "", "(kg/m³)", "(kg/m³)")
    println(io, "-" ^ 60)
    for i in eachindex(r.pressures)
        @printf(io, "%14.4e %10.4f %10.4f %10.2f %12.2f\n",
            r.pressures[i], r.relative_volume[i], r.V_factors[i],
            r.density_oil[i], r.density_gas[i])
    end
end

function Base.show(io::IO, r::DLEResult)
    println(io, "DLE Experiment Results")
    println(io, "=" ^ 60)
    @printf(io, "Temperature: %.2f K (%.2f °C)\n", r.T, r.T - 273.15)
    @printf(io, "Bubble point pressure: %.4e Pa (%.2f bar)\n", r.p_sat, r.p_sat / 1e5)
    println(io, "-" ^ 60)
    @printf(io, "%14s %10s %10s %10s %10s\n", "P (Pa)", "Bo", "Rs", "ρ_oil", "Bg")
    @printf(io, "%14s %10s %10s %10s %10s\n", "", "(m³/m³)", "(m³/m³)", "(kg/m³)", "(m³/m³)")
    println(io, "-" ^ 60)
    for i in eachindex(r.pressures)
        @printf(io, "%14.4e %10.4f %10.4f %10.2f %10.6f\n",
            r.pressures[i], r.Bo[i], r.Rs[i],
            r.density_oil[i], r.Bg[i])
    end
end

function Base.show(io::IO, r::CVDResult)
    println(io, "CVD Experiment Results")
    println(io, "=" ^ 60)
    @printf(io, "Temperature: %.2f K (%.2f °C)\n", r.T, r.T - 273.15)
    @printf(io, "Dew point pressure: %.4e Pa (%.2f bar)\n", r.p_sat, r.p_sat / 1e5)
    println(io, "-" ^ 60)
    @printf(io, "%14s %10s %12s %10s\n", "P (Pa)", "Z_gas", "Liq dropout", "Gas prod")
    println(io, "-" ^ 60)
    for i in eachindex(r.pressures)
        @printf(io, "%14.4e %10.4f %12.4f %10.4f\n",
            r.pressures[i], r.Z_gas[i],
            r.liquid_dropout[i], r.gas_produced_cumulative[i])
    end
end

function Base.show(io::IO, r::MSSResult)
    println(io, "Multistage Separator Test Results")
    println(io, "=" ^ 60)
    @printf(io, "Number of stages: %d\n", length(r.stages))
    @printf(io, "Bo = %.4f m³/m³\n", r.Bo)
    @printf(io, "Rs = %.4f m³/m³\n", r.Rs)
    @printf(io, "Stock tank oil density: %.2f kg/m³\n", r.density_oil_st)
    @printf(io, "API gravity: %.1f\n", r.API_gravity)
    println(io, "-" ^ 60)
    @printf(io, "%6s %14s %10s %12s\n", "Stage", "P (Pa)", "T (K)", "GOR (m³/m³)")
    println(io, "-" ^ 60)
    for i in eachindex(r.stages)
        @printf(io, "%6d %14.4e %10.2f %12.4f\n",
            i, r.stages[i].p, r.stages[i].T, r.GOR_stage[i])
    end
end

function Base.show(io::IO, t::PVTOTable)
    println(io, "PVTO Table")
    println(io, "=" ^ 60)
    @printf(io, "%10s %14s %12s %12s\n", "Rs", "P (Pa)", "Bo", "μ_o (Pa·s)")
    println(io, "-" ^ 60)
    for i in eachindex(t.Rs)
        for j in eachindex(t.p[i])
            rs_str = j == 1 ? @sprintf("%10.4f", t.Rs[i]) : " " ^ 10
            @printf(io, "%s %14.4e %12.6f %12.4e\n",
                rs_str, t.p[i][j], t.Bo[i][j], t.mu_o[i][j])
        end
        println(io, "/")
    end
end

function Base.show(io::IO, t::PVDGTable)
    println(io, "PVDG Table")
    println(io, "=" ^ 60)
    @printf(io, "%14s %12s %12s\n", "P (Pa)", "Bg", "μ_g (Pa·s)")
    println(io, "-" ^ 60)
    for i in eachindex(t.p)
        @printf(io, "%14.4e %12.6f %12.4e\n", t.p[i], t.Bg[i], t.mu_g[i])
    end
end

function Base.show(io::IO, t::PVTGTable)
    println(io, "PVTG Table")
    println(io, "=" ^ 60)
    @printf(io, "%14s %10s %12s %12s\n", "P (Pa)", "Rv", "Bg", "μ_g (Pa·s)")
    println(io, "-" ^ 60)
    for i in eachindex(t.p)
        for j in eachindex(t.Rv_sub[i])
            p_str = j == 1 ? @sprintf("%14.4e", t.p[i]) : " " ^ 14
            @printf(io, "%s %10.6f %12.6f %12.4e\n",
                p_str, t.Rv_sub[i][j], t.Bg[i][j], t.mu_g[i][j])
        end
        println(io, "/")
    end
end

function Base.show(io::IO, s::SurfaceDensities)
    println(io, "Surface Densities")
    @printf(io, "  Oil: %.2f kg/m³\n", s.oil)
    @printf(io, "  Gas: %.4f kg/m³\n", s.gas)
end
