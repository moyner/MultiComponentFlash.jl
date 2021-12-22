
"""
    inverse_flash_update!(storage, eos, c, V)

Update internal matrix of partial derivatives for a converged flash result.
"""
function inverse_flash_update!(storage, eos, c, V)
    x, y = storage.x, storage.y
    J_p, = update_primary_jacobian!(storage, eos, c, storage.forces, x, y, V)
    J_s, = update_secondary_jacobian!(storage, eos, c, storage.forces_secondary, x, y, V)
    invert_sp!(J_s, J_p)
    return J_s
end

function invert_sp!(J_s, J_p)
    F = lu!(J_p)
    ldiv!(F, J_s)
end

function invert_sp!(J_s::MMatrix, J_p::MMatrix)
    tmp = SMatrix(J_p)\SMatrix(J_s)
    J_s .= tmp
end


"""
    set_partials(v, storage, eos, c, index)

Modify a value `v` of ForwardDiff.Dual type to get the correct derivatives.

If the mixture is two-phase, the partial derivatives of phase molar fractions and vapor fraction
with respect to the flash conditions (pressure, temperature and overall mole fractions) are not
trivial. This function is the main gateway for setting these values.

# Notes
Experimental interface, subject to change. You most likely want to use either
[`set_partials_phase_mole_fractions!`](@ref) or [`set_partials_vapor_fraction`](@ref)

In order for this routine to work, the storage must be initialized using `flash_storage`
with the following options enabled:
- `inc_jac = true`
- `diff_externals = true`
and `inverse_flash_update!` must be called after a successful flash. The partial derivatives
with respect to p, T, z is then contained in `storage.buf_inv` with a negative sign. This
function then performs the requisite chain rule operations for the input.
"""
function set_partials(v, storage, eos, c, index)
    M = storage.J_inv
    p, z, T = c.p, c.z, c.T
    buf = storage.buf_inv
    @. buf = 0
    is_dual(x) = typeof(x)<:ForwardDiff.Dual
    if is_dual(p)
        @inbounds @. buf -= M[index, 1].*p.partials
    end
    if is_dual(T)
        @inbounds @. buf -= M[index, 2].*T.partials
    end
    if eltype(z)<:ForwardDiff.Dual
        for (ind, zi) = enumerate(z)
            @inbounds ∂z = -M[index, ind+2]
            for i in eachindex(zi.partials)
                @inbounds buf[i] += ∂z*zi.partials[i]
            end
        end
    end
    ∂T = typeof(v)
    P = typeof(v.partials)

    return ∂T(v.value, P(Tuple(buf)))
end

"""
    set_partials_vapor_fraction(V, storage, eos, ∂c)

Set partial derivatives to a vapor mole fraction `ForwardDiff.Dual` instance with correct value, but
missing partial derivatives.
"""
set_partials_vapor_fraction(V, storage, eos, ∂c) = set_partials(V, storage, eos, ∂c, 2*number_of_components(eos) + 1)

"""
set_partials_phase_mole_fractions!(xy, storage, eos, ∂c, [phase_symbol])

Set partial derivatives to phase mole fraction vector with type `ForwardDiff.Dual` with correct values, but
missing partial derivatives.
"""
function set_partials_phase_mole_fractions!(xy, storage, eos, ∂c, phase = :liquid)
    n = length(xy)
    if phase == :liquid
        offset = 0
    else
        offset = n
    end
    for i = 1:n
        xy[i] = set_partials(xy[i], storage, eos, ∂c, offset + i)
    end
    return xy
end
