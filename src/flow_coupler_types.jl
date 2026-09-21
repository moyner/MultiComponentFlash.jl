export PhaseState2Phase, liquid_phase_present, vapor_phase_present, two_phases_present, phase_data

@enum PhaseState2Phase two_phase_lv single_phase_l single_phase_v unknown_phase_state_lv

"Type that holds values for a flashed phase (mole fractions + compressibility factor)"
struct FlashedPhase{T, A<:AbstractVector{T}}
    mole_fractions::A
    Z::T
    function FlashedPhase(mole_fractions::AbstractVector, Z::Tz) where Tz
        T = Base.promote_type(Tz, eltype(mole_fractions))
        mole_fractions = map(x -> convert(T, x), mole_fractions)
        Z = convert(T, Z)
        new{T, typeof(mole_fractions)}(mole_fractions, Z)
    end
end

function FlashedPhase(n::Integer, T::DataType = Float64)
    x = zeros(T, n)
    Z = zero(T)
    return FlashedPhase(x, Z)
end

function Base.convert(::Type{FlashedPhase{T, Vector{T}}}, ph::FlashedPhase{K, Vector{K}}) where {T<:AbstractFloat, K<:ForwardDiff.Dual}
    F = x -> convert(T, value(x))
    mf = map(F, ph.mole_fractions)
    Z = F(ph.Z)
    return FlashedPhase(mf, Z)
end

function Base.convert(::Type{FlashedPhase{T, Vector{T}}}, ph::FlashedPhase{K, Vector{K}}) where {T, K}
    F = x -> convert(T, x)
    mf = map(F, ph.mole_fractions)
    Z = F(ph.Z)
    return FlashedPhase(mf, Z)
end

function Base.convert(::Type{FlashedPhase{T, A}}, ph::FlashedPhase) where {T, A<:AbstractVector{T}}
    fractions = convert(A, map(x -> convert(T, x), ph.mole_fractions))
    return FlashedPhase(fractions, convert(T, ph.Z))
end

"Type that holds liquid and vapor phase states together with their state"
struct FlashedMixture2Phase{T, A<:AbstractVector{T}, E}
    state::PhaseState2Phase
    K::E # Equilibrium constants
    V::T # Vapor mole fraction
    liquid::FlashedPhase{T, A}
    vapor::FlashedPhase{T, A}
    critical_distance::eltype(E)
    flash_cond::NamedTuple{(:p, :T, :z), Tuple{eltype(E), eltype(E), E}}
    flash_stability::StabilityReport
    function FlashedMixture2Phase(state::PhaseState2Phase, K::K_t, V::V_t,
            liquid::FlashedPhase{V_t, A}, vapor::FlashedPhase{V_t, A};
            vec_type = A,
            critical_distance = NaN,
            cond = missing,
            stability_report = StabilityReport()
        ) where {V_t, K_t, A}
        # The phase vector type is part of the result type. Infer it from the
        # phases instead of using the keyword value as a type parameter.
        vec_type === A || throw(ArgumentError("vec_type must match the phase vectors"))
        R = eltype(K_t)
        if ismissing(cond)
            z0 = convert(K_t, fill(NaN, length(K)))
            cond = (p = R(NaN), T = R(NaN), z = z0)
        end
        cond = (p = R(cond.p), T = R(cond.T), z = convert(K_t, cond.z))
        new{V_t, A, K_t}(state, K, V, liquid, vapor,
            R(critical_distance), cond, stability_report)
    end
end

function Base.convert(::Type{FlashedMixture2Phase{T, Vector{T}, F}}, mixture::FlashedMixture2Phase{K, Vector{K}, F}) where {T<:ForwardDiff.Dual, K, F}
    T_phase = FlashedPhase{T, Vector{T}}
    liquid = convert(T_phase, mixture.liquid)
    vapor = convert(T_phase, mixture.vapor)

    to_T = x -> convert(T, x)
    V = to_T(mixture.V)
    Kv = mixture.K
    converted_mixture = FlashedMixture2Phase(
        mixture.state,
        Kv,
        V,
        liquid,
        vapor;
        vec_type = Vector{T},
        critical_distance = mixture.critical_distance,
        cond = mixture.flash_cond,
        stability_report = mixture.flash_stability
    )
    return converted_mixture
end

function Base.convert(::Type{FlashedMixture2Phase{T, A, E}},
        mixture::FlashedMixture2Phase) where {T, A<:AbstractVector{T}, E}
    liquid = convert(FlashedPhase{T, A}, mixture.liquid)
    vapor = convert(FlashedPhase{T, A}, mixture.vapor)
    K = convert(E, mixture.K)
    cond = mixture.flash_cond
    flash_cond = (p = cond.p, T = cond.T, z = convert(E, cond.z))
    return FlashedMixture2Phase(mixture.state, K, convert(T, mixture.V),
        liquid, vapor;
        vec_type = A,
        critical_distance = mixture.critical_distance,
        cond = flash_cond,
        stability_report = mixture.flash_stability)
end

function FlashedMixture2Phase(state, K, V, x, y, Z_L, Z_V, b = NaN, cond = missing, stability = StabilityReport())
    T = typeof(V)
    liquid = FlashedPhase(map(a -> convert(T, a), x), convert(T, Z_L))
    vapor = FlashedPhase(map(a -> convert(T, a), y), convert(T, Z_V))
    return FlashedMixture2Phase(state, K, V, liquid, vapor,
        vec_type = typeof(liquid.mole_fractions), critical_distance = b, cond = cond,
        stability_report = stability)
end

function FlashedMixture2Phase(eos::AbstractEOS, T = Float64, T_num = Float64, b = NaN, cond = missing, stability = StabilityReport())
    n = number_of_components(eos)
    V = zero(T)
    K = zeros(T_num, n)
    liquid = FlashedPhase(n, T)
    vapor = FlashedPhase(n, T)
    return FlashedMixture2Phase(unknown_phase_state_lv, K, V, liquid, vapor, critical_distance = b, cond = cond, stability_report = stability)
end
