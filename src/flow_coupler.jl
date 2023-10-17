export PhaseState2Phase, liquid_phase_present, vapor_phase_present, two_phases_present, phase_data

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
struct FlashedPhase{T, A<:AbstractVector{T}}
    mole_fractions::A
    Z::T
    function FlashedPhase(mole_fractions::AbstractVector, Z::Real)
        new{typeof(Z), typeof(mole_fractions)}(mole_fractions, Z)
    end
end

function FlashedPhase(n::Integer, T::DataType = Float64)
    x = zeros(T, n)
    Z = zero(T)
    return FlashedPhase(x, Z)
end

"Type that holds liquid and vapor phase states together with their state"
struct FlashedMixture2Phase{T, A<:AbstractVector{T}, E}
    state::PhaseState2Phase
    K::E # Equilibrium constants
    V::T # Vapor mole fraction
    liquid::FlashedPhase{T, A}
    vapor::FlashedPhase{T, A}
    function FlashedMixture2Phase(state::PhaseState2Phase, K::K_t, V::V_t, liquid, vapor; vec_type = Vector{V_t}) where {V_t, K_t}        
        new{V_t, vec_type, K_t}(state, K, V, liquid, vapor)
    end
end

function FlashedMixture2Phase(state, K, V, x, y, Z_L, Z_V)
    liquid = FlashedPhase(x, Z_L)
    vapor = FlashedPhase(y, Z_V)
    return FlashedMixture2Phase(state, K, V, liquid, vapor)
end

function FlashedMixture2Phase(eos::AbstractEOS, T = Float64, T_num = Float64)
    n = number_of_components(eos)
    V = zero(T)
    # K values are always doubles
    K = zeros(T_num, n)
    liquid = FlashedPhase(n, T)
    vapor = FlashedPhase(n, T)

    return FlashedMixture2Phase(unknown_phase_state_lv, K, V, liquid, vapor)
end

function flashed_mixture_2ph(eos, cond, K = initial_guess_K(eos, cond); method = SSIFlash(), kwarg...)
    # Convenience function for getting flashed phases
    S = flash_storage(eos, cond, method = method)
    return flashed_mixture_2ph!(S, eos, cond, K; method = method, kwarg...)
end

function flashed_mixture_2ph!(storage, eos, conditions, K; kwarg...)
    V, K, rep = flash_2ph!(storage, K, eos, conditions; kwarg..., extra_out = true)
    V = single_phase_label(eos.mixture, conditions)
    (; p, T, z) = conditions
    forces = storage.forces
    x = @. liquid_mole_fraction(z, K, V)
    y = @. vapor_mole_fraction(x, K)
    if V == 0 || V == 1
        if V == 0
            state = single_phase_l
        else
            state = single_phase_v
        end
        liquid = vapor = cond
        Z_L = Z_V = mixture_compressibility_factor(eos, conditions, forces)
    else
        state = two_phase_lv
        liquid = (p = p, T = T, z = x)
        vapor = (p = p, T = T, z = y)
        Z_L = mixture_compressibility_factor(eos, liquid, forces)
        Z_V = mixture_compressibility_factor(eos, vapor, forces)
    end
    return FlashedMixture2Phase(state, K, V, x, y, Z_L, Z_V)
end

function phase_data(mix::FlashedMixture2Phase, phase)
    if phase == :liquid
        return mix.liquid
    else
        return mix.vapor
    end
end

"""
    phase_saturations(eos, p, T, flashed_mixture)

Compute phase saturations for a flashed two-phase mixture.

Always returns a named tuple of (S_l, S_v), even if the mixture is single-phase.

The value in the absent phase will be zero.
"""
@inline function phase_saturations(eos, p, Temp, f::FlashedMixture2Phase{T}) where T
    state = f.state
    # @assert state != unknown_phase_state_lv "Phase state is not known. Cannot compute saturations. Has flash been called?."
    if state == two_phase_lv
        S_v = two_phase_vapor_saturation(eos, p, Temp, f)
    elseif state == single_phase_v
        S_v = one(T)
    else
        S_v = zero(T)
    end
    S_l = one(T) - S_v
    return (S_l = S_l::T, S_v = S_v::T)
end

@inline function two_phase_vapor_saturation(eos, p, Temp, f::FlashedMixture2Phase{T}) where T
    state = f.state
    # A faster definition that doesn't go via molar volume, but assumes no shifts:
    # Z_l = f.liquid.Z
    # Z_v = f.vapor.Z
    # S_v = Z_v*V/(Z_l*(1-V) + Z_v*V)
    V = f.V
    L = one(V) - V
    vol_v = V*molar_volume(eos, p, Temp, f.vapor)
    vol_l = L*molar_volume(eos, p, Temp, f.liquid)
    S_v = vol_v/(vol_v + vol_l)
    return S_v
end

"Compute molar volume of a flashed phase"
@inline function molar_volume(eos, p, T, ph::FlashedPhase) 
    V = molar_volume(IDEAL_GAS_CONSTANT, p, T, ph.Z)
    V = correct_volume(V, eos, p, T, ph, eos.volume_shift)
    return V
end

function correct_volume(V, eos, p, T, ph, volume_shift)
    R = IDEAL_GAS_CONSTANT
    xy = ph.mole_fractions
    cond = (p = p, T = T, z = xy)
    c = zero(V)
    @inbounds for i in eachindex(xy)
        prp = molecular_property(eos.mixture, i)
        T_ci = critical_temperature(prp)
        P_ci = critical_pressure(prp)
        ω_bi = weight_bi(eos, cond, i);
        c += volume_shift[i]*xy[i]*ω_bi*R*T_ci/P_ci 
    end
    return V - c
end

@inline correct_volume(V, eos, p, T, ph, volume_shift::Nothing) = V

"Compute mass density of a flashed phase"
function mass_density(eos, p, T, ph::FlashedPhase{V}) where V
    props = eos.mixture.properties
    z = ph.mole_fractions
    t = zero(V)
    @inbounds for (i, p) in enumerate(props)
        t += z[i]*p.mw
    end
    ρ = t/molar_volume(eos, p, T, ph)
    return convert(V, ρ)
end

"""
    mass_densities(eos, p, T, flashed_mixture)

Compute mass densities for a flashed two-phase mixture.

Always returns a named tuple of (ρ_l, ρ_v), even if the mixture is single-phase.

The value in the absent phase will be zero.
"""
@inline function mass_densities(eos, p, temperature, f::FlashedMixture2Phase{T}) where T
    state = f.state
    # @assert state != unknown_phase_state_lv "Phase state is not known. Cannot compute densities. Has flash been called?."
    if liquid_phase_present(state)
        l = mass_density(eos, p, temperature, f.liquid)::T
    else
        l = zero(T)
    end
    if vapor_phase_present(state)
        v = mass_density(eos, p, temperature, f.vapor)::T
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
    # @assert state != unknown_phase_state_lv "Phase state is not known. Cannot compute viscosities. Has flash been called?."
    if liquid_phase_present(state)
        l = lbc_viscosity(eos, p, temperature, f.liquid)::T
    else
        l = convert(T, eps(T))
    end
    if vapor_phase_present(state)
        v = lbc_viscosity(eos, p, temperature, f.vapor)::T
    else
        v = convert(T, eps(T))
    end
    return (μ_l = l::T, μ_v = v::T)
end
