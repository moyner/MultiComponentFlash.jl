export PhaseState2Phase, liquid_phase_present, vapor_phase_present, two_phases_present

@enum PhaseState2Phase two_phase_lv single_phase_l single_phase_v unknown_phase_state_lv

@inline liquid_phase_present(state::PhaseState2Phase) = state == two_phase_lv || state == single_phase_l
@inline vapor_phase_present(state::PhaseState2Phase) = state == two_phase_lv || state == single_phase_v
@inline two_phases_present(state::PhaseState2Phase) = state == two_phase_lv

"""
    phase_is_present(label, phase_state)

Check if a phase (symbol :liquid/:vapor) is present with the provided phase state.
"""
function phase_is_present(label, phase_state::PhaseState2Phase)
    if label == :liquid
        present = liquid_phase_present(phase_state)
    elseif label == :vapor
        present = vapor_phase_present(phase_state)
    else
        present = false
    end
    return present
end

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
struct FlashedMixture2Phase{T, A}
    state::PhaseState2Phase
    K::A # Equilibrium constants
    V::T # Vapor mole fraction
    liquid::FlashedPhase{T}
    vapor::FlashedPhase{T}
    function FlashedMixture2Phase(state::PhaseState2Phase, K, V, liquid, vapor)
        new{typeof(V), typeof(K)}(state, K, V, liquid, vapor)
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

    return FlashedMixture2Phase(unknown_phase_state_lv, K, V, liquid, vapor)
end

"""
    phase_saturations(flashed_mixture)

Compute phase saturations for a flashed two-phase mixture.

Always returns a named tuple of (S_l, S_v), even if the mixture is single-phase.

The value in the absent phase will be zero.
"""
@inline function phase_saturations(f::FlashedMixture2Phase{T}) where T
    state = f.state
    @assert state != unknown_phase_state_lv "Phase state is not known. Cannot compute saturations. Has flash been called?."
    if state == two_phase_lv
        Z_l = f.liquid.Z
        Z_v = f.vapor.Z
        V = f.V
        S_v = Z_v*V/(Z_l*(1-V) + Z_v*V)
    elseif state == single_phase_v
        S_v = one(T)
    else
        S_v = zero(T)
    end
    S_l = 1 - S_v
    return (S_l = S_l::T, S_v = S_v::T)
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
@inline function mass_densities(eos, p, temperature, f::FlashedMixture2Phase{T}) where T
    state = f.state
    @assert state != unknown_phase_state_lv "Phase state is not known. Cannot compute densities. Has flash been called?."
    if liquid_phase_present(state)
        l = mass_density(eos, p, temperature, f.liquid)
    else
        l = zero(T)
    end
    if vapor_phase_present(state)
        v = mass_density(eos, p, temperature, f.vapor)
    else
        v = zero(T)
    end
    return (ρ_l = l::T, ρ_v = v::T)
end

"""
    lbc_viscosities(eos, p, T, flashed_mixture)

Compute phase viscosities for a flashed two-phase mixture using the LBC correlation.

Always returns a named tuple of (μ_l, μ_v), even if the mixture is single-phase.

The value in the absent phase will be a tiny value (eps of the numeric type) to make
division for mobilities etc. safe.
"""
@inline function lbc_viscosities(eos, p, temperature, f::FlashedMixture2Phase{T}) where T
    state = f.state
    @assert state != unknown_phase_state_lv "Phase state is not known. Cannot compute viscosities. Has flash been called?."
    if liquid_phase_present(state)
        l = lbc_viscosity(eos, p, temperature, f.liquid)
    else
        l = convert(T, eps(T))
    end
    if vapor_phase_present(state)
        v = lbc_viscosity(eos, p, temperature, f.vapor)
    else
        v = convert(T, eps(T))
    end
    return (μ_l = l::T, μ_v = v::T)
end
