"""
    dle(eos, z, T; p_range, n_points, p_sc, T_sc)

Perform a Differential Liberation Expansion (DLE) experiment.

Starting at the bubble point, the pressure is reduced in steps. At each step,
the equilibrium gas is completely removed and the liquid is re-equilibrated at
the next lower pressure. Properties (Bo, Rs, etc.) are referenced to the
final stock-tank oil volume at standard conditions.

# Arguments
- `eos`: Equation of state
- `z`: Overall mole fractions of initial fluid
- `T`: Reservoir temperature (K)
- `p_range`: (p_max, p_min) for the experiment. Default: (50e6, 1e5)
- `n_points`: Number of pressure steps below the saturation point. Default: 20
- `p_sc`: Standard condition pressure (Pa). Default: 101325.0
- `T_sc`: Standard condition temperature (K). Default: 288.706 (60°F)
"""
function dle(eos, z, T;
        p_range = (50e6, 1e5),
        n_points = 20,
        p_sc = 101325.0,
        T_sc = 288.706
    )
    z = collect(Float64, z)
    z ./= sum(z)
    nc = number_of_components(eos)

    # Find saturation pressure
    p_sat, is_bp = find_saturation_pressure(eos, z, T;
        p_min = p_range[2], p_max = p_range[1])

    if !is_bp
        @warn "DLE is designed for bubble point fluids. Found dew point fluid."
    end

    # Generate pressure steps from p_sat down to p_min
    pressures = collect(range(p_sat, p_range[2], length = n_points + 1))

    # Storage
    Bo_arr = zeros(n_points + 1)
    Rs_arr = zeros(n_points + 1)
    ρ_oil_arr = zeros(n_points + 1)
    ρ_gas_arr = zeros(n_points + 1)
    μ_oil_arr = zeros(n_points + 1)
    μ_gas_arr = zeros(n_points + 1)
    Bg_arr = zeros(n_points + 1)
    Z_gas_arr = zeros(n_points + 1)
    gas_comp = Vector{Vector{Float64}}(undef, n_points + 1)
    oil_comp = Vector{Vector{Float64}}(undef, n_points + 1)

    # Track remaining oil and gas
    z_oil = copy(z)  # Current oil composition
    n_oil = 1.0      # Moles of oil remaining (normalized)

    # Flash the residual oil to standard conditions to get reference volume
    # First, step through the liberation
    moles_gas_removed = Float64[]  # Gas removed at each step
    oil_molar_volumes = Float64[]  # Oil molar volume at reservoir conditions
    gas_molar_volumes_sc = Float64[]  # Gas molar volume at standard conditions

    for (i, p) in enumerate(pressures)
        if i == 1
            # At bubble point - single phase liquid
            props = flash_and_properties(eos, p, T, z_oil)
            gas_comp[i] = copy(z_oil)
            oil_comp[i] = copy(z_oil)
            ρ_oil_arr[i] = props.ρ_l
            ρ_gas_arr[i] = props.ρ_v
            μ_oil_arr[i] = props.μ_l
            μ_gas_arr[i] = props.μ_v
            Z_gas_arr[i] = props.Z_V
            Bg_arr[i] = 0.0
            push!(oil_molar_volumes, props.V_mol_l)
            push!(moles_gas_removed, 0.0)
        else
            # Flash at this pressure
            props = flash_and_properties(eos, p, T, z_oil)

            gas_comp[i] = copy(props.y)
            oil_comp[i] = copy(props.x)
            ρ_oil_arr[i] = props.ρ_l
            ρ_gas_arr[i] = props.ρ_v
            μ_oil_arr[i] = props.μ_l
            μ_gas_arr[i] = props.μ_v
            Z_gas_arr[i] = props.Z_V
            Bg_arr[i] = phase_molar_volume_from_Z(p, T, props.Z_V) /
                         phase_molar_volume_from_Z(p_sc, T_sc, 1.0)
            push!(oil_molar_volumes, props.V_mol_l)

            # Gas removed at this step
            n_gas = n_oil * props.V
            push!(moles_gas_removed, n_gas)

            # Update oil: only liquid remains
            n_oil = n_oil * (1.0 - props.V)
            z_oil = copy(props.x)
        end
    end

    # Get stock tank oil properties: flash remaining oil at standard conditions
    props_sc = flash_and_properties(eos, p_sc, T_sc, z_oil)
    V_mol_oil_sc = props_sc.V_mol_l
    mw_oil_sc = mixture_molar_mass(props_sc.x, eos)

    # Reference: 1 mol of original fluid at bubble point
    # Volume of stock tank oil per mole of original fluid
    V_oil_sc_per_mol_init = n_oil * (1.0 - props_sc.V) * V_mol_oil_sc

    # Compute Bo and Rs
    # Bo = volume of oil at reservoir conditions / volume of stock tank oil
    # Rs = volume of gas at standard conditions / volume of stock tank oil
    V_gas_sc_per_mol = phase_molar_volume_from_Z(p_sc, T_sc, 1.0)  # ideal gas at SC

    # Cumulative gas remaining in solution from step i downward
    cum_gas_removed_below = zeros(n_points + 1)
    for i in (n_points + 1):-1:2
        cum_gas_removed_below[i] = 0.0
    end
    # Gas still in solution at step i = total gas from step i+1 to end
    gas_in_solution = zeros(n_points + 1)
    running_gas = 0.0
    for i in (n_points + 1):-1:1
        gas_in_solution[i] = running_gas
        if i > 1
            running_gas += moles_gas_removed[i]
        end
    end

    for (i, p) in enumerate(pressures)
        if V_oil_sc_per_mol_init > 0.0
            # Bo: volume of oil at (p,T) per unit stock tank oil volume
            # At step i, we have n_oil_at_i moles of oil
            # Compute n_oil_at_i
            n_oil_i = 1.0
            for j in 2:i
                n_oil_i *= (1.0 - flash_and_properties(eos, pressures[j], T, oil_comp[max(1, j - 1)]).V)
            end
            if i == 1
                Bo_arr[i] = oil_molar_volumes[i] / V_oil_sc_per_mol_init
            else
                Bo_arr[i] = n_oil_i * oil_molar_volumes[i] / V_oil_sc_per_mol_init
            end

            # Rs: gas in solution at standard conditions / stock tank oil volume
            Rs_arr[i] = gas_in_solution[i] * V_gas_sc_per_mol / V_oil_sc_per_mol_init
        end
    end

    return DLEResult(T, collect(Float64, z ./ sum(z)), pressures,
        Bo_arr, Rs_arr, ρ_oil_arr, ρ_gas_arr, μ_oil_arr, μ_gas_arr,
        Bg_arr, Z_gas_arr, gas_comp, oil_comp, p_sat)
end
