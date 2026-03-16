"""
    cvd(eos, z, T; p_range, n_points, p_sc, T_sc)

Perform a Constant Volume Depletion (CVD) experiment.

Starting at the dew point, the pressure is reduced in steps. At each step,
gas is removed to restore the cell to its original volume. The liquid dropout
and produced gas properties are recorded.

This experiment is typical for gas condensate systems.

# Arguments
- `eos`: Equation of state
- `z`: Overall mole fractions
- `T`: Temperature (K)
- `p_range`: (p_max, p_min) for the experiment. Default: (50e6, 1e5)
- `n_points`: Number of pressure steps. Default: 20
- `p_sc`: Standard condition pressure (Pa). Default: 101325.0
- `T_sc`: Standard condition temperature (K). Default: 288.706 (60°F)
"""
function cvd(eos, z, T;
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

    if is_bp
        @warn "CVD is designed for dew point (gas condensate) fluids. Found bubble point fluid."
    end

    # Pressure steps from p_sat down
    pressures = collect(range(p_sat, p_range[2], length = n_points + 1))

    # Storage
    Z_factors = zeros(n_points + 1)
    liq_dropout = zeros(n_points + 1)
    gas_produced_cum = zeros(n_points + 1)
    ρ_gas_arr = zeros(n_points + 1)
    μ_gas_arr = zeros(n_points + 1)
    Z_gas_arr = zeros(n_points + 1)
    gas_comp = Vector{Vector{Float64}}(undef, n_points + 1)

    # Initial state at dew point: single phase gas
    props_init = flash_and_properties(eos, p_sat, T, z)
    V_cell = props_init.V_mol_v  # Reference cell volume per mole

    n_total = 1.0  # Total moles in cell
    z_cell = copy(z)  # Current overall composition in cell
    cum_gas_produced = 0.0  # Cumulative gas produced (moles)

    for (i, p) in enumerate(pressures)
        if i == 1
            # At dew point - single phase gas
            Z_factors[i] = props_init.Z_V
            liq_dropout[i] = 0.0
            gas_produced_cum[i] = 0.0
            ρ_gas_arr[i] = props_init.ρ_v
            μ_gas_arr[i] = props_init.μ_v
            Z_gas_arr[i] = props_init.Z_V
            gas_comp[i] = copy(z)
        else
            # Flash at lower pressure
            props = flash_and_properties(eos, p, T, z_cell)

            Z_gas_arr[i] = props.Z_V
            gas_comp[i] = copy(props.y)
            ρ_gas_arr[i] = props.ρ_v
            μ_gas_arr[i] = props.μ_v

            # Two-phase Z factor
            Z_factors[i] = (1.0 - props.V) * props.Z_L + props.V * props.Z_V

            # Volumes
            V_liq = n_total * (1.0 - props.V) * props.V_mol_l
            V_gas = n_total * props.V * props.V_mol_v
            V_total = V_liq + V_gas

            # Liquid dropout = liquid volume / cell volume
            liq_dropout[i] = V_liq / V_cell

            # Gas to remove to restore cell volume
            V_excess = V_total - V_cell
            if V_excess > 0 && props.V_mol_v > 0
                n_gas_remove = V_excess / props.V_mol_v
            else
                n_gas_remove = 0.0
            end

            cum_gas_produced += n_gas_remove
            gas_produced_cum[i] = cum_gas_produced

            # Update cell: remove gas, keep liquid + remaining gas
            n_gas_remaining = n_total * props.V - n_gas_remove
            n_liq = n_total * (1.0 - props.V)

            if n_gas_remaining < 0
                n_gas_remaining = 0.0
            end

            n_total_new = n_liq + n_gas_remaining

            # New composition
            if n_total_new > 0
                z_cell = (n_liq .* props.x .+ n_gas_remaining .* props.y) ./ n_total_new
                z_cell ./= sum(z_cell)
            end
            n_total = n_total_new
        end
    end

    # Normalize cumulative gas produced
    gas_produced_cum ./= max(gas_produced_cum[end], 1e-30)

    return CVDResult(T, collect(Float64, z ./ sum(z)), pressures,
        Z_factors, liq_dropout, gas_produced_cum,
        ρ_gas_arr, μ_gas_arr, Z_gas_arr, gas_comp, p_sat)
end
