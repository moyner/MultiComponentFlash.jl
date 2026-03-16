"""
    mss(eos, z_feed, T_res, stages::Vector{SeparatorStage}; p_res)

Perform a Multistage Separator Test (MSS).

The feed fluid at reservoir conditions is flashed through a sequence of
separator stages at specified (p, T) conditions. Gas is removed at each stage,
and the liquid proceeds to the next stage. The final stage is the stock tank.

Returns Bo, Rs (referenced to stock tank oil), and properties at each stage.

# Arguments
- `eos`: Equation of state
- `z_feed`: Feed composition (overall mole fractions at reservoir conditions)
- `T_res`: Reservoir temperature (K)
- `stages`: Vector of `SeparatorStage(p, T)` (last stage is stock tank)
- `p_res`: Reservoir pressure (Pa). Default: will be determined from saturation pressure.
"""
function mss(eos, z_feed, T_res, stages::Vector{SeparatorStage};
        p_res = nothing
    )
    z_feed = collect(Float64, z_feed)
    z_feed ./= sum(z_feed)
    nc = number_of_components(eos)
    n_stages = length(stages)

    # Find reservoir pressure if not given
    if isnothing(p_res)
        p_res_val, _ = find_saturation_pressure(eos, z_feed, T_res)
    else
        p_res_val = p_res
    end

    # Flash feed through separator stages
    sep_result = separator_flash(eos, z_feed, stages)

    # Compute Bo: oil volume at reservoir conditions / stock tank oil volume
    # Oil at reservoir conditions (1 mole basis)
    props_res = flash_and_properties(eos, p_res_val, T_res, z_feed)
    V_mol_oil_res = props_res.V_mol_l

    # Track liquid mole fraction through stages
    n_oil = 1.0  # Start with 1 mole of feed
    liquid_fraction_remaining = 1.0
    GOR_stage = zeros(n_stages)

    p_sc = stages[end].p
    T_sc = stages[end].T

    for i in 1:n_stages
        V_i = sep_result.V_stages[i]
        n_gas_i = liquid_fraction_remaining * V_i
        n_liq_i = liquid_fraction_remaining * (1.0 - V_i)

        # GOR at this stage: gas volume at SC / oil volume at SC
        if n_liq_i > 0
            V_gas_sc = n_gas_i * phase_molar_volume_from_Z(p_sc, T_sc, 1.0)
            GOR_stage[i] = V_gas_sc
        end

        liquid_fraction_remaining = n_liq_i
    end

    # Stock tank oil properties
    z_st = sep_result.oil_compositions[end]
    props_st = flash_and_properties(eos, p_sc, T_sc, z_st)
    V_mol_oil_sc = props_st.V_mol_l
    ρ_oil_st = props_st.ρ_l

    # Stock tank oil volume per initial mole
    V_oil_st = liquid_fraction_remaining * V_mol_oil_sc

    # Bo = oil volume at reservoir / stock tank oil volume
    Bo = V_mol_oil_res / V_oil_st

    # Rs = total gas volume at SC / stock tank oil volume
    V_gas_total_sc = sum(GOR_stage)
    Rs = V_gas_total_sc / V_oil_st

    # Normalize GOR per stage by stock tank oil volume
    GOR_stage ./= V_oil_st

    # Gas density at stock tank
    z_gas_st = sep_result.gas_compositions[end]
    mw_gas_st = mixture_molar_mass(z_gas_st, eos)
    V_mol_gas_sc = phase_molar_volume_from_Z(p_sc, T_sc, 1.0)
    ρ_gas_st = mw_gas_st / V_mol_gas_sc

    # API gravity
    sg_oil = ρ_oil_st / WATER_DENSITY_SC  # specific gravity relative to water
    API = 141.5 / sg_oil - 131.5

    return MSSResult(stages, z_feed, Bo, Rs,
        sep_result.gas_compositions,
        sep_result.oil_compositions,
        GOR_stage, ρ_oil_st, ρ_gas_st, API)
end
