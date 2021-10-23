"""
    MolecularProperty(molar_mass, p_crit, T_crit, V_crit, acentric_factor = 0.0)

Type that defines the static properties of a molecular species.
"""
struct MolecularProperty{R}
    mw::R
    p_c::R
    T_c::R
    V_c::R
    ω::R
    function MolecularProperty(mw, p_c, T_c, V_c, ω = 0.0)
        @assert p_c >= 0
        @assert mw > 0
        @assert T_c > 0
        @assert V_c > 0
        @assert -1.0 <= ω "Expected -1 ≤ ω, was $ω"
        new{typeof(mw)}(mw, p_c, T_c, V_c, ω)
    end
end

"""
    MolecularProperty("Name")

Convenience constructor that looks up molecular properties from a table in `tabulated_properties`.

The properties are taken from the wonderful MIT-licensed [CoolProp](http://coolprop.org). Please note that
the equations of state included in this module may not be approprioate for all the available fluids, especially for mixtures!

See list of species [at CoolProp website](http://coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids).

"""
function MolecularProperty(s::String)
    t = tabulated_properties();
    @assert haskey(t, s)
    return t[s]
end
"""Get the critical pressure for a species."""
@inline critical_pressure(M::MolecularProperty{R}) where R = M.p_c::R
"""Get the critical temperature for a species."""
@inline critical_temperature(M::MolecularProperty{R}) where R = M.T_c::R
"""Get the critical volume for a species."""
@inline critical_volume(M::MolecularProperty{R}) where R = M.V_c::R
"""Get the molar weight for a species."""
@inline molar_weight(M::MolecularProperty{R}) where R = M.mw::R
"""Get the acentric factorfor a species."""
@inline acentric_factor(M::MolecularProperty{R}) where R = M.ω::R

"""
    MultiComponentMixture(properties::NTuple{N, MolecularProperty}; A_ij = nothing, names = ["C1", "C2", ...], name = "UnnamedMixture")

Create a multicomponent mixture with an optional binary interaction coefficient matrix `A_ij`.
"""
struct MultiComponentMixture{R, N}
    name::String
    component_names::Vector{String}
    properties::NTuple{N, MolecularProperty{R}}
    binary_interaction::Union{Matrix{R}, Nothing}
    function MultiComponentMixture(properties; A_ij = nothing, names = ["C$d" for d in 1:length(properties)], name = "UnnamedMixture")
        properties = Tuple(properties)
        realtype = typeof(properties[1].mw)
        n = length(properties)
        if !isnothing(A_ij)
            @assert size(A_ij) == (n, n)
        end
        new{realtype, n}(name, names, properties, A_ij)
    end
end

"""
    number_of_components(mixture)

Return number of components in the `MultiComponentMixture`.
"""
number_of_components(m::MultiComponentMixture{R, N}) where {R, N} = N
function molecular_properties(mixture::MultiComponentMixture{R, N}) where {R, N}
    return mixture.properties::NTuple{N, MolecularProperty{R}}
end
@inline molecular_property(mixture::MultiComponentMixture{R}, i) where R = @inbounds molecular_properties(mixture)[i]::MolecularProperty{R}
@inline reduced_pressure(mixture::MultiComponentMixture, cond, i) = cond.p/critical_pressure(molecular_property(mixture, i))
@inline reduced_temperature(mixture::MultiComponentMixture, cond, i) = cond.T/critical_temperature(molecular_property(mixture, i))
@inline reduced_pT(arg...) = (reduced_pressure(arg...), reduced_temperature(arg...))


abstract type AbstractEOS end
abstract type AbstractCubicEOS <: AbstractEOS end

@inline reduced_pT(mix, c, i) = (reduced_pressure(mix, c, i), reduced_temperature(mix, c, i))
@inline reduced_pT(eos::AbstractEOS, c, i) = (reduced_pressure(eos.mixture, c, i), reduced_temperature(eos.mixture, c, i))

"""
GenericCubicEOS is an implementation of generalized cubic equations of state.

Many popular cubic equations can be written in a single form with a few changes in
definitions for the terms (they are, after all, all cubic in form). References:

 1. [Cubic Equations of State-Which? by J.J. Martin](https://doi.org/10.1021/i160070a001)
 2. [Simulation of Gas Condensate Reservoir Performance  by K.H. Coats](https://doi.org/10.2118/10512-PA)

"""
struct GenericCubicEOS{T, R, N} <: AbstractCubicEOS
    type::T
    mixture::MultiComponentMixture{R, N}
    m_1::R
    m_2::R
    ω_a::R
    ω_b::R
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

"""
Specializes the GenericCubicEOS to the Peng-Robinson cubic equation of state.
"""
struct PengRobinson <: AbstractGeneralizedCubic end
"""
Specializes the GenericCubicEOS to the Zudkevitch-Joffe cubic equation of state.

The Zudkevitch-Joffe equations of state allows for per-component functions of 
temperature that modify the `weight_ai` and `weight_bi` functions. These additional
fitting parameters allows for more flexibility when matching complex mixtures.
"""
struct ZudkevitchJoffe <: AbstractGeneralizedCubic
    F_a
    F_b
    function ZudkevitchJoffe(; F_a = (T, c_i) -> 1.0, F_b = (T, c_i) -> 1.0)
        new(F_a, F_b)
    end
end

"""
Specializes the GenericCubicEOS to the Soave-Redlich-Kwong cubic equation of state.
"""
struct SoaveRedlichKwong <: AbstractGeneralizedCubic end
"""
Specializes the GenericCubicEOS to the Redlich-Kwong cubic equation of state.
"""
struct RedlichKwong <: AbstractGeneralizedCubic end

"""
    GenericCubicEOS(mixture, [type = PengRobinson()])

Instantiate a generic cubic equation-of-state for a `MultiComponentMixture` and 
a specified EOS.

Currently supported choices for type:

    1. `PengRobinson` (default)
    2. `ZudkevitchJoffe`
    3. `RedlichKwong`
    4. `SoaveRedlichKwong`
"""
function GenericCubicEOS(mixture, type = PengRobinson())
    ω_a, ω_b, m_1, m_2 = static_coefficients(type)
    setup = (ω_a = ω_a, ω_b = ω_b, m_1 = m_1, m_2 = m_2, type = type)
    return GenericCubicEOS(setup, mixture)
end

function GenericCubicEOS(setup::NamedTuple, mixture)
    return GenericCubicEOS(setup.type, mixture, setup.m_1, setup.m_2, setup.ω_a, setup.ω_b)
end
