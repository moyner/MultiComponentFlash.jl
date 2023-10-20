export PhaseState2Phase, liquid_phase_present, vapor_phase_present, two_phases_present, phase_data

@enum PhaseState2Phase two_phase_lv single_phase_l single_phase_v unknown_phase_state_lv

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
