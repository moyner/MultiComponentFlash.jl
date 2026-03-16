"""
    cce(eos, z, T; p_range, n_points)

Perform a Constant Composition Expansion (CCE) experiment.

The mixture with composition `z` is held at temperature `T` (K) and the
pressure is decreased from a high value through the saturation pressure.
The volume is measured relative to the volume at the saturation point.

# Arguments
- `eos`: Equation of state
- `z`: Overall mole fractions
- `T`: Temperature (K)
- `p_range`: Tuple (p_max, p_min) in Pa. Default: (50e6, 1e5)
- `n_points`: Number of pressure steps. Default: 40
"""
function cce(eos, z, T;
        p_range = (50e6, 1e5),
        n_points = 40
    )
    z = collect(Float64, z)
    z ./= sum(z)
    nc = number_of_components(eos)

    # Find saturation pressure
    p_sat, is_bp = find_saturation_pressure(eos, z, T;
        p_min = p_range[2], p_max = p_range[1])

    # Generate pressure steps
    pressures = collect(range(p_range[1], p_range[2], length = n_points))
    sort!(pressures, rev = true)

    # Compute properties at saturation point for reference volume
    props_sat = flash_and_properties(eos, p_sat, T, z)
    V_mol_sat = (1.0 - props_sat.V) * props_sat.V_mol_l + props_sat.V * props_sat.V_mol_v

    # Storage
    Z_factors = zeros(n_points)
    V_factors = zeros(n_points)
    rel_vol = zeros(n_points)
    ρ_oil = zeros(n_points)
    ρ_gas = zeros(n_points)
    μ_oil = zeros(n_points)
    μ_gas = zeros(n_points)

    for (i, p) in enumerate(pressures)
        props = flash_and_properties(eos, p, T, z)
        V_factors[i] = props.V
        Z_factors[i] = (1.0 - props.V) * props.Z_L + props.V * props.Z_V

        # Total molar volume at this pressure
        V_mol = (1.0 - props.V) * props.V_mol_l + props.V * props.V_mol_v
        rel_vol[i] = V_mol / V_mol_sat

        ρ_oil[i] = props.ρ_l
        ρ_gas[i] = props.ρ_v
        μ_oil[i] = props.μ_l
        μ_gas[i] = props.μ_v
    end

    return CCEResult(T, z, pressures, Z_factors, V_factors, rel_vol,
        ρ_oil, ρ_gas, μ_oil, μ_gas, p_sat, is_bp)
end
