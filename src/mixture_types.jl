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
    MolecularProperty(; mw, p_c, T_c, V_c, acentric_factor = 0.0)

Keyword constructor version of `MolecularProperty`. Except for the acentric
factor, all properties must be specified. Explanation of inputs:

mw: Molar mass (kg / mol)
p_c: Critical pressure (Pa)
T_c: Critical temperature (°K)
V_c: Critical volume (m^3 / mol)
acentric_factor (dimensionless)
"""
function MolecularProperty(; mw, p_c, T_c, V_c, acentric_factor = 0.0)
    return MolecularProperty(mw, p_c, T_c, V_c, acentric_factor)
end

"""
    MolecularProperty("Name")

Convenience constructor that looks up molecular properties from a table in
`tabulated_properties`.

The properties are taken from the wonderful MIT-licensed
[CoolProp](http://coolprop.org). Please note that the equations of state
included in this module may not be appropriate for all the available species
included in CoolProp, especially for mixtures!

See list of species [at CoolProp
website](http://coolprop.org/fluid_properties/PurePseudoPure.html#list-of-fluids).

"""
function MolecularProperty(s::String)
    t = tabulated_properties();
    haskey(t, s) || throw(ArgumentError("$s not found in `tabulated_properties()`."))
    return t[s]
end

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
    MultiComponentMixture(names::Vector{String}; A_ij = nothing, name = "UnnamedMixture")

Create a multicomponent mixture using name lookup for species.
"""
function MultiComponentMixture(names::Vector{String}; kwarg...)
    props = map(MolecularProperty, names)
    return MultiComponentMixture(props; names = names, kwarg...)
end
