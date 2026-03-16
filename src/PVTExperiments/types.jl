"""
    SeparatorStage(p, T)

A single separator stage defined by pressure `p` (Pa) and temperature `T` (K).
"""
struct SeparatorStage
    p::Float64
    T::Float64
end

"""
    CCEResult

Results from a Constant Composition Expansion experiment.

# Fields
- `T`: Temperature (K)
- `z`: Overall composition (mole fractions)
- `pressures`: Pressure steps (Pa)
- `Z_factors`: Two-phase compressibility factors
- `V_factors`: Vapor mole fraction at each step
- `relative_volume`: Relative volume V/V_sat
- `density_oil`: Oil (liquid) density at each step (kg/m³)
- `density_gas`: Gas (vapor) density at each step (kg/m³)
- `viscosity_oil`: Oil viscosity (Pa·s)
- `viscosity_gas`: Gas viscosity (Pa·s)
- `p_sat`: Bubble/dew point pressure (Pa)
- `is_bubblepoint`: true if bubble point, false if dew point
"""
struct CCEResult
    T::Float64
    z::Vector{Float64}
    pressures::Vector{Float64}
    Z_factors::Vector{Float64}
    V_factors::Vector{Float64}
    relative_volume::Vector{Float64}
    density_oil::Vector{Float64}
    density_gas::Vector{Float64}
    viscosity_oil::Vector{Float64}
    viscosity_gas::Vector{Float64}
    p_sat::Float64
    is_bubblepoint::Bool
end

"""
    DLEResult

Results from a Differential Liberation Expansion experiment.

# Fields
- `T`: Temperature (K)
- `z_init`: Initial overall composition (mole fractions)
- `pressures`: Pressure steps (Pa) (from p_sat down to final pressure)
- `Bo`: Oil formation volume factor (m³/m³ at standard conditions)
- `Rs`: Solution gas-oil ratio (m³/m³ at standard conditions)
- `density_oil`: Oil density at each step (kg/m³)
- `density_gas`: Gas density at each step (kg/m³)
- `viscosity_oil`: Oil viscosity (Pa·s)
- `viscosity_gas`: Gas viscosity (Pa·s)
- `Bg`: Gas formation volume factor (m³/m³ at standard conditions)
- `Z_gas`: Gas compressibility factor
- `gas_composition`: Gas composition at each liberation step
- `oil_composition`: Oil composition at each liberation step
- `p_sat`: Bubble point pressure (Pa)
"""
struct DLEResult
    T::Float64
    z_init::Vector{Float64}
    pressures::Vector{Float64}
    Bo::Vector{Float64}
    Rs::Vector{Float64}
    density_oil::Vector{Float64}
    density_gas::Vector{Float64}
    viscosity_oil::Vector{Float64}
    viscosity_gas::Vector{Float64}
    Bg::Vector{Float64}
    Z_gas::Vector{Float64}
    gas_composition::Vector{Vector{Float64}}
    oil_composition::Vector{Vector{Float64}}
    p_sat::Float64
end

"""
    CVDResult

Results from a Constant Volume Depletion experiment.

# Fields
- `T`: Temperature (K)
- `z_init`: Initial overall composition (mole fractions)
- `pressures`: Pressure steps (Pa)
- `Z_factors`: Two-phase Z factors
- `liquid_dropout`: Liquid volume fraction (fraction of cell volume)
- `gas_produced_cumulative`: Cumulative gas produced (mole fraction of initial)
- `density_gas`: Gas density at each step (kg/m³)
- `viscosity_gas`: Gas viscosity (Pa·s)
- `Z_gas`: Gas phase compressibility factor
- `gas_composition`: Gas composition at each step
- `p_sat`: Dew point pressure (Pa)
"""
struct CVDResult
    T::Float64
    z_init::Vector{Float64}
    pressures::Vector{Float64}
    Z_factors::Vector{Float64}
    liquid_dropout::Vector{Float64}
    gas_produced_cumulative::Vector{Float64}
    density_gas::Vector{Float64}
    viscosity_gas::Vector{Float64}
    Z_gas::Vector{Float64}
    gas_composition::Vector{Vector{Float64}}
    p_sat::Float64
end

"""
    MSSResult

Results from a Multistage Separator Test.

# Fields
- `stages`: Separator stages (p, T pairs)
- `z_feed`: Feed composition (mole fractions)
- `Bo`: Oil formation volume factor from reservoir to stock tank (m³/m³)
- `Rs`: Solution GOR from reservoir to stock tank (m³/m³)
- `gas_composition`: Gas composition at each stage
- `oil_composition`: Oil composition exiting each stage
- `GOR_stage`: Gas-oil ratio at each separator stage (m³/m³)
- `density_oil_st`: Stock tank oil density (kg/m³)
- `density_gas_st`: Gas density at stock tank (kg/m³)
- `API_gravity`: API gravity of stock tank oil
"""
struct MSSResult
    stages::Vector{SeparatorStage}
    z_feed::Vector{Float64}
    Bo::Float64
    Rs::Float64
    gas_composition::Vector{Vector{Float64}}
    oil_composition::Vector{Vector{Float64}}
    GOR_stage::Vector{Float64}
    density_oil_st::Float64
    density_gas_st::Float64
    API_gravity::Float64
end

"""
    PVTOTable

Black oil PVTO (live oil) table for dissolved gas oil systems.

# Fields
- `Rs`: Solution gas-oil ratio values (m³/m³ at surface conditions)
- `p_bub`: Bubble point pressures for each Rs (Pa)
- `p`: Pressure values for undersaturated oil at each Rs level (Pa)
- `Bo`: Oil formation volume factor (m³/m³)
- `mu_o`: Oil viscosity (Pa·s)
"""
struct PVTOTable
    Rs::Vector{Float64}
    p_bub::Vector{Float64}
    p::Vector{Vector{Float64}}
    Bo::Vector{Vector{Float64}}
    mu_o::Vector{Vector{Float64}}
end

"""
    PVDGTable

Black oil PVDG (dry gas) table.

# Fields
- `p`: Pressure values (Pa)
- `Bg`: Gas formation volume factor (m³/m³)
- `mu_g`: Gas viscosity (Pa·s)
"""
struct PVDGTable
    p::Vector{Float64}
    Bg::Vector{Float64}
    mu_g::Vector{Float64}
end

"""
    PVTGTable

Black oil PVTG (wet gas / gas condensate) table.

# Fields
- `p`: Pressure values (Pa)
- `Rv`: Vaporized oil-gas ratio at each pressure (m³/m³ at surface conditions)
- `Rv_sub`: Undersaturated Rv values at each pressure level
- `Bg`: Gas formation volume factor for each (p, Rv) entry (m³/m³)
- `mu_g`: Gas viscosity for each (p, Rv) entry (Pa·s)
"""
struct PVTGTable
    p::Vector{Float64}
    Rv::Vector{Float64}
    Rv_sub::Vector{Vector{Float64}}
    Bg::Vector{Vector{Float64}}
    mu_g::Vector{Vector{Float64}}
end

"""
    SurfaceDensities

Surface densities of oil and gas at standard conditions.

# Fields
- `oil`: Oil density at surface conditions (kg/m³)
- `gas`: Gas density at surface conditions (kg/m³)
"""
struct SurfaceDensities
    oil::Float64
    gas::Float64
end
