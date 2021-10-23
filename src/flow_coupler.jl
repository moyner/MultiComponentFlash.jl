abstract type AbstractPhaseState end
abstract type AbstractTwoPhaseState <: AbstractPhaseState end
abstract type AbstractSinglePhaseState <: AbstractPhaseState end

"Single-phase liquid state for dispatch"
struct SinglePhaseLiquid <: AbstractSinglePhaseState end
"Single-phase vapor state for dispatch"
struct SinglePhaseVapor <: AbstractSinglePhaseState end
"Two-phase liquid-vapor state for dispatch"
struct TwoPhaseLiquidVapor <: AbstractTwoPhaseState end
"Unknown phase state (not initialized)"
struct UnknownPhaseState <: AbstractPhaseState end

"""
    phase_is_present(label, phase_state)

Check if a phase (symbol :liquid/:vapor) is present with the provided phase state.
"""
phase_is_present(label, phase_state) = false
phase_is_present(label, ::TwoPhaseLiquidVapor) = true
phase_is_present(label, ::SinglePhaseLiquid) = label == :liquid
phase_is_present(label, ::SinglePhaseVapor) = label == :vapor

"Type that holds values for a flashed phase (mole fractions + compressibility factor)"
struct FlashedPhase{T}
    mole_fractions::AbstractArray{T}
    Z::T
    function FlashedPhase(mole_fractions::AbstractArray, Z::Real)
        new{typeof(Z)}(mole_fractions, Z)
    end
end

function FlashedPhase(n::Integer, T::DataType = Float64)
    x = zeros(T, n)
    Z = zero(T)
    return FlashedPhase(x, Z)
end

"Type that holds liquid and vapor phase states together with their state"
struct FlashedMixture2Phase
    state::AbstractPhaseState
    K # Equilibrium constants
    V # Vapor mole fractions
    liquid::FlashedPhase
    vapor::FlashedPhase
    function FlashedMixture2Phase(state::AbstractPhaseState, K, V, liquid, vapor)
        new(state, K, V, liquid, vapor)
    end
end

function FlashedMixture2Phase(state, K, V, x, y, Z_L, Z_V)
    liquid = FlashedPhase(x, Z_L)
    vapor = FlashedPhase(y, Z_V)
    return FlashedMixture2Phase(state, K, V, liquid, vapor)
end

function FlashedMixture2Phase(eos::AbstractEOS, T = Float64)
    n = number_of_components(eos)
    V = zero(T)
    # K values are always doubles
    K = zeros(Float64, n)
    liquid = FlashedPhase(n, T)
    vapor = FlashedPhase(n, T)

    return FlashedMixture2Phase(UnknownPhaseState(), K, V, liquid, vapor)
end

"""
    phase_saturations(flashed_mixture)

Compute phase saturations for a flashed two-phase mixture.

Always returns a named tuple of (S_l, S_v), even if the mixture is single-phase.

The value in the absent phase will be zero.
"""
phase_saturations(r::FlashedMixture2Phase) = phase_saturations(r, r.state)
phase_saturations(f, ::SinglePhaseVapor) = (S_l = 0.0, S_v = 1.0)
phase_saturations(f, ::SinglePhaseLiquid) = (S_l = 1.0, S_v = 0.0)
phase_saturations(f, ::UnknownPhaseState) = error("Phase state is not known. No flash results available.")

function phase_saturations(f, ::TwoPhaseLiquidVapor)
    Z_l = f.liquid.Z
    Z_v = f.vapor.Z
    V = f.V
    S_v = Z_v*V/(Z_l*(1-V) + Z_v*V)
    return (S_l = 1 - S_v, S_v = S_v)
end

"Compute molar volume of a flashed phase"
molar_volume(eos, p, T, ph::FlashedPhase) = molar_volume(IDEAL_GAS_CONSTANT, p, T, ph.Z)

"Compute mass density of a flashed phase"
function mass_density(eos, p, T, ph::FlashedPhase)
    props = eos.mixture.properties
    z = ph.mole_fractions
    t = 0.0
    for (i, p) in enumerate(props)
        t += z[i]*p.mw
    end
    ρ = t/molar_volume(eos, p, T, ph)
    return ρ
end

"""
    mass_densities(eos, p, T, flashed_mixture)

Compute mass densities for a flashed two-phase mixture.

Always returns a named tuple of (ρ_l, ρ_v), even if the mixture is single-phase.

The value in the absent phase will be zero.
"""
mass_densities(eos, p, T, r::FlashedMixture2Phase) = mass_densities(eos, p, T, r, r.state)

function mass_densities(eos, p, T, f, ::SinglePhaseVapor)
    v = mass_density(eos, p, T, f.vapor)
    (ρ_l = 0.0, ρ_v = v)
end

function mass_densities(eos, p, T, f, ::SinglePhaseLiquid)
    l = mass_density(eos, p, T, f.liquid)
    (ρ_l = l, ρ_v = 0.0)
end

function mass_densities(eos, p, T, f, ::TwoPhaseLiquidVapor)
    l = mass_density(eos, p, T, f.liquid)
    v = mass_density(eos, p, T, f.vapor)
    (ρ_l = l, ρ_v = v)
end

mass_densities(eos, p, T, f, ::UnknownPhaseState) = error("Phase state is not known. No flash results available.")

"""
    lbc_viscosities(eos, p, T, flashed_mixture)

Compute phase viscosities for a flashed two-phase mixture using the LBC correlation.

Always returns a named tuple of (μ_l, μ_v), even if the mixture is single-phase.

The value in the absent phase will be zero.
"""
lbc_viscosities(eos, p, T, r::FlashedMixture2Phase) = lbc_viscosities(eos, p, T, r, r.state)

function lbc_viscosities(eos, p, T, f, ::SinglePhaseLiquid)
    l = lbc_viscosity(eos, p, T, f.liquid)
    (μ_l = l, μ_v = 0.0)
end

function lbc_viscosities(eos, p, T, f, ::SinglePhaseVapor)
    v = lbc_viscosity(eos, p, T, f.vapor)
    (μ_l = 0.0, μ_v = v)
end

function lbc_viscosities(eos, p, T, f, ::TwoPhaseLiquidVapor)
    l = lbc_viscosity(eos, p, T, f.liquid)
    v = lbc_viscosity(eos, p, T, f.vapor)
    (μ_l = l, μ_v = v)
end

lbc_viscosities(eos, p, T, f, ::UnknownPhaseState) =  error("Phase state is not known. No flash results available.")

