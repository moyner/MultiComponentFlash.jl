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


