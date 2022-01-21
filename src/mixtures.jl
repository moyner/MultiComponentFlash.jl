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
            @assert size(A_ij) == (n, n) "Binary interaction coefficients must be a $n x $n matrix."
            @assert issymmetric(A_ij) "Binary interaction coefficients must be symmetric if provided."
            A_ij = Symmetric(A_ij)
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


