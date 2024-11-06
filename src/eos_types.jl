abstract type AbstractEOS end
abstract type AbstractCubicEOS <: AbstractEOS end

"""
GenericCubicEOS is an implementation of generalized cubic equations of state.

Many popular cubic equations can be written in a single form with a few changes in
definitions for the terms (they are, after all, all cubic in form). References:

 1. [Cubic Equations of State-Which? by J.J. Martin](https://doi.org/10.1021/i160070a001)
 2. [Simulation of Gas Condensate Reservoir Performance  by K.H. Coats](https://doi.org/10.2118/10512-PA)

"""
struct GenericCubicEOS{T, R, N, V} <: AbstractCubicEOS
    type::T
    mixture::MultiComponentMixture{R, N}
    m_1::R
    m_2::R
    ω_a::R
    ω_b::R
    volume_shift::V
end

"""
Abstract base type for generalized cubics.

Any subtype may override the any of following internal functions to alter the EOS:

    * static_coefficients
    * weight_ai
    * weight_bi

In addition, the GenericCubicEOS is parametric on the specific cubic type for further
dispatch modifications.
"""
abstract type AbstractGeneralizedCubic end


abstract type AbstractPengRobinson <: AbstractGeneralizedCubic end

"""
    pr = PengRobinson()

Specializes the GenericCubicEOS to the Peng-Robinson cubic equation of state.
"""
struct PengRobinson <: AbstractPengRobinson end


"""
    prc = PengRobinsonCorrected()

Specializes the GenericCubicEOS to Peng-Robinson modified for large acentric factors
"""

struct PengRobinsonCorrected <: AbstractPengRobinson end
"""
    zj = ZudkevitchJoffe(; F_a = (T, c_i) -> 1.0, F_b = (T, c_i) -> 1.0)

Specializes the GenericCubicEOS to the Zudkevitch-Joffe cubic equation of state.

The Zudkevitch-Joffe equations of state allows for per-component functions of 
temperature that modify the `weight_ai` and `weight_bi` functions. These additional
fitting parameters allows for more flexibility when matching complex mixtures.
"""
struct ZudkevitchJoffe{A, B} <: AbstractGeneralizedCubic
    F_a::A
    F_b::B
    function ZudkevitchJoffe(; F_a = (T, c_i) -> 1.0, F_b = (T, c_i) -> 1.0)
        new{typeof(F_a), typeof(F_b)}(F_a, F_b)
    end
end

"""
    srk = SoaveRedlichKwong()

Specializes the GenericCubicEOS to the Soave-Redlich-Kwong cubic equation of state.
"""
struct SoaveRedlichKwong <: AbstractGeneralizedCubic end
"""
    rk = RedlichKwong()

Specializes the GenericCubicEOS to the Redlich-Kwong cubic equation of state.
"""
struct RedlichKwong <: AbstractGeneralizedCubic end

"""
    GenericCubicEOS(mixture)
    GenericCubicEOS(mixture, PengRobinson())


Instantiate a generic cubic equation-of-state for a `MultiComponentMixture` and
a specified EOS.

Currently supported choices for second argument:

    1. `PengRobinson()` (default)
    2. `PengRobinsonCorrected()`
    3. `ZudkevitchJoffe()`
    4. `RedlichKwong()`
    5. `SoaveRedlichKwong()`
"""
function GenericCubicEOS(mixture, type = PengRobinson(); kwarg...)
    ω_a, ω_b, m_1, m_2 = static_coefficients(type)
    setup = (ω_a = ω_a, ω_b = ω_b, m_1 = m_1, m_2 = m_2, type = type)
    return GenericCubicEOS(setup, mixture; kwarg...)
end

function GenericCubicEOS(setup::NamedTuple, mixture; volume_shift = nothing)
    if !isnothing(volume_shift)
        length(volume_shift) == number_of_components(mixture) || throw(ArgumentError("Volume shift must have one value per component."))
    end
    return GenericCubicEOS(setup.type, mixture, setup.m_1, setup.m_2, setup.ω_a, setup.ω_b, volume_shift)
end

struct KValuesEOS{T, R, N, V} <: AbstractEOS
    "Callable on the form `cond -> V` or a set of constants (Tuple/AbstractVector)"
    K_values_evaluator::T
    mixture::MultiComponentMixture{R, N}
    volume_shift::V
end

function KValuesEOS(K, mixture; volume_shift = nothing)
    if !isnothing(volume_shift)
        length(volume_shift) == number_of_components(mixture) || throw(ArgumentError("Volume shift must have one value per component."))
    end
    return KValuesEOS(K, mixture, volume_shift)
end

function eostype(::KValuesEOS)
    return KValuesEOS
end
