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

function flashed_mixture_2ph(eos, cond, K = initial_guess_K(eos, cond); method = SSIFlash(), kwarg...)
    # Convenience function for getting flashed phases
    S = flash_storage(eos, cond, method = method)
    return flashed_mixture_2ph!(S, eos, cond, K; method = method, kwarg...)
end

function flashed_mixture_2ph!(storage, eos, conditions, K; kwarg...)
    V, K, rep = flash_2ph!(storage, K, eos, conditions; kwarg..., extra_out = true)
    p = conditions.p
    T = conditions.T
    z = conditions.z
    forces = storage.forces
    if isnan(V)
        V = single_phase_label(eos.mixture, conditions)
        if V == 0
            state = single_phase_l
        else
            state = single_phase_v
        end
        liquid = vapor = cond
        Z_L = Z_V = mixture_compressibility_factor(eos, conditions, forces)
        x = copy(z)
        y = copy(z)
    else
        x = @. liquid_mole_fraction(z, K, V)
        y = @. vapor_mole_fraction(x, K)
        state = two_phase_lv
        liquid = (p = p, T = T, z = x)
        vapor = (p = p, T = T, z = y)
        Z_L = mixture_compressibility_factor(eos, liquid, forces)
        Z_V = mixture_compressibility_factor(eos, vapor, forces)
    end
    return FlashedMixture2Phase(state, K, V, x, y, Z_L, Z_V)
end

function phase_data(mix::FlashedMixture2Phase{T, A, E}, phase) where {T, A, E}
    if phase == :liquid
        out = mix.liquid
    else
        out = mix.vapor
    end
    return out::FlashedPhase{T, A}
end

"""
    phase_saturations(eos, p, T, flashed_mixture)

Compute phase saturations for a flashed two-phase mixture.

Always returns a named tuple of (S_l, S_v), even if the mixture is single-phase.

The value in the absent phase will be zero.
"""
@inline function phase_saturations(eos, p, Temp, f::FlashedMixture2Phase{T}) where T
    state = f.state::PhaseState2Phase
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
        ω_bi = weight_bi(eos, cond, i)
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

The value in the absent phase will be that of the present phase for single-phase conditions.
"""
@inline function mass_densities(eos, p, temperature, f::FlashedMixture2Phase{T}) where T
    state = f.state
    # @assert state != unknown_phase_state_lv "Phase state is not known. Cannot compute densities. Has flash been called?."
    has_liquid = liquid_phase_present(state)
    has_vapor = vapor_phase_present(state)
    if has_liquid && has_vapor
        l = mass_density(eos, p, temperature, f.liquid)::T
        v = mass_density(eos, p, temperature, f.vapor)::T
    elseif has_liquid
        l = v = mass_density(eos, p, temperature, f.liquid)::T
    else
        l = v = mass_density(eos, p, temperature, f.vapor)::T
    end
    return (ρ_l = l::T, ρ_v = v::T)
end

"""
    lbc_viscosities(eos, p, T, flashed_mixture)

Compute phase viscosities for a flashed two-phase mixture using the LBC correlation.

Always returns a named tuple of (μ_l, μ_v), even if the mixture is single-phase.

The value in the absent phase will be that of the present phase for single-phase conditions.
"""
@inline function lbc_viscosities(eos, p, temperature, f::FlashedMixture2Phase{T}) where T
    state = f.state
    # @assert state != unknown_phase_state_lv "Phase state is not known. Cannot compute viscosities. Has flash been called?."
    has_liquid = liquid_phase_present(state)
    has_vapor = vapor_phase_present(state)
    if has_liquid && has_vapor
        l = lbc_viscosity(eos, p, temperature, f.liquid)::T
        v = lbc_viscosity(eos, p, temperature, f.vapor)::T
    elseif has_liquid
        l = v = lbc_viscosity(eos, p, temperature, f.liquid)::T
    else
        l = v = lbc_viscosity(eos, p, temperature, f.vapor)::T
    end
    return (μ_l = l::T, μ_v = v::T)
end
