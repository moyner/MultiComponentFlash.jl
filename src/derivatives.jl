import ForwardDiff: Dual, Partials, value
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
    return set_partials(v, M, buf, p, T, z, index)
end

function set_partials(v::Dual{T,V,N}, M, buf, p, temp, z, index) where {T,V,N}
    # Zero out buffer just in case
    @. buf = 0
    # ∂v/∂p
    if isa(p, Dual)
        set_partials_scalar!(buf, p, M, index, 1)
    end
    # ∂v/∂T
    if isa(temp, Dual)
        set_partials_scalar!(buf, temp, M, index, 2)
    end
    # ∂v/∂z_i for all i
    if eltype(z)<:Dual
        set_partials_vector!(buf, z, M, index, 2)
    end
    ∂ = Partials{N, V}(Tuple(buf))
    val = value(v)
    return Dual{T, V, N}(val, ∂)
end

set_partials_scalar!(buf, X::AbstractFloat, M, index, pos) = nothing

function set_partials_scalar!(buf, X::ForwardDiff.Dual{T,V,N}, M, index, pos) where {T, V, N}
    @inbounds @. buf -= M[index, pos]*X.partials
end

set_partials_vector!(buf, z, M, index, offset) = nothing

function set_partials_vector!(buf, z::AbstractVector{ForwardDiff.Dual{T,V,N}}, M, index, offset) where {T, V, N}
    for (ind, zi) = enumerate(z)
        @inbounds ∂z = -M[index, ind+offset]
        for i in eachindex(zi.partials)
            @inbounds buf[i] += ∂z*zi.partials[i]
        end
    end
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
    @inbounds for i = 1:n
        xy[i] = set_partials(xy[i], storage, eos, ∂c, offset + i)
    end
    return xy
end
