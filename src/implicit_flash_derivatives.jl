module ImplicitFlashDerivatives

using LinearAlgebra
using ForwardDiff
using StaticArrays: SVector, SMatrix
using ..MultiComponentFlash: GenericCubicEOS, KValuesEOS,
    liquid_mole_fraction, vapor_mole_fraction, objectiveRR, objectiveRR_dV,
    static_force_coefficients, static_fugacity_coefficients, prep

export implicit_flash_derivatives

"""
    implicit_flash_derivatives(eos::GenericCubicEOS, numeric_cond, ad_cond,
                               V_numeric, K_numeric)

Differentiate an already converged, interior two-phase cubic flash. Return
`(V_ad, K_ad)`, with derivatives of the vapor fraction and equilibrium K-values
with respect to the AD inputs in `ad_cond`.

`numeric_cond` and `ad_cond` each contain `p`, `T`, and `z`. Their primal values
must agree: `numeric_cond` is the numeric condition used to solve the flash,
while any of `ad_cond.p`, `ad_cond.T`, and `ad_cond.z` may carry AD derivatives.
Both `z` fields and `K_numeric` must be `SVector`s of the EOS component count.
`V_numeric` and `K_numeric` are the converged numeric vapor fraction and
equilibrium K-values at `numeric_cond`; this function does not solve a flash.

The function differentiates the fugacity equilibrium and Rachford-Rice
residuals, then solves their static `(N+1) × (N+1)` Jacobian system. It requires
`0 < V_numeric < 1`, positive compositions and K-values, and a nonsingular
Jacobian. It does not differentiate stability decisions or phase boundaries.
For accelerator use, pass an immutable EOS from [`make_eos_immutable`](@ref).
"""
@inline function implicit_flash_derivatives(
        eos::GenericCubicEOS{E, R, N}, numeric_cond, ad_cond,
        V_numeric, K_numeric::SVector{N}) where {E, R, N}
    u = static_flash_unknowns(K_numeric, V_numeric)

    # Seed only the equilibrium unknowns for the primary Jacobian. The
    # condition here is numeric, so no nested AD types are needed.
    u_seeded = seed_static_flash_unknowns(u)
    residual_seeded = static_equilibrium_residual(eos, numeric_cond, u_seeded)
    jacobian = static_flash_jacobian(residual_seeded)
    residual_numeric = static_flash_residual_values(residual_seeded)

    # The constant residual is subtracted so the returned primals remain the
    # converged numeric solution, even when its tolerance is finite.
    T = promote_type(typeof(ad_cond.p), typeof(ad_cond.T), eltype(ad_cond.z))
    ad_cond = (p = convert(T, ad_cond.p), T = convert(T, ad_cond.T),
        z = SVector{N, T}(ad_cond.z))
    u_ad = SVector{N+1, T}(u)
    residual_ad = static_equilibrium_residual(eos, ad_cond, u_ad)
    correction = jacobian \ (residual_ad - residual_numeric)
    updated = SVector{N+1, T}(u_ad - correction)
    return updated[N+1], static_flash_K(updated)
end

"""
    implicit_flash_derivatives(eos::KValuesEOS, numeric_cond, ad_cond,
                               V_numeric, K_numeric, K_ad)

Differentiate an already converged, interior two-phase K-value flash. Return
`(V_ad, K_ad)` using the implicit derivative of the Rachford-Rice equation.

`numeric_cond` and `ad_cond` each contain `p`, `T`, and `z` with matching primal
values. Both `z` fields, `K_numeric`, and `K_ad` must be `SVector`s of the EOS
component count. `V_numeric` is the converged numeric vapor fraction for
`numeric_cond` and `K_numeric`. Evaluate `K_ad` at `ad_cond` before calling this
function; its primal values must match `K_numeric`. Pressure and temperature
derivatives enter through `K_ad`, while composition derivatives enter through
`ad_cond.z`.

This function does not solve a flash or evaluate K-values. It requires
`0 < V_numeric < 1`, positive K-values, and a nonzero Rachford-Rice slope. It
does not differentiate phase boundaries.
"""
@inline function implicit_flash_derivatives(
        eos::KValuesEOS{E, R, N}, numeric_cond, ad_cond,
        V_numeric, K_numeric::SVector{N}, K_ad::SVector{N}) where {E, R, N}
    residual_numeric = objectiveRR(V_numeric, K_numeric, numeric_cond.z)
    T = promote_type(typeof(ad_cond.p), typeof(ad_cond.T),
        eltype(ad_cond.z), eltype(K_ad))
    K = SVector{N, T}(K_ad)
    z = SVector{N, T}(ad_cond.z)
    residual_ad = objectiveRR(V_numeric, K, z)
    slope = objectiveRR_dV(V_numeric, K_numeric, numeric_cond.z)
    V = convert(T, V_numeric) - convert(T,
        (residual_ad - residual_numeric)/slope)
    return V, K
end

@generated function static_flash_unknowns(K::SVector{N, F}, V) where {N, F}
    entries = Any[:(K[$i]) for i in 1:N]
    push!(entries, :(V))
    return :(SVector{$(N+1)}(($(entries...),)))
end

@generated function seed_static_flash_unknowns(
        u::SVector{M, F}) where {M, F}
    entries = [:(D(u[$i], ForwardDiff.single_seed(
        ForwardDiff.Partials{$M, F}, Val($i)))) for i in 1:M]
    return quote
        D = ForwardDiff.Dual{Nothing, F, $M}
        SVector{$M, D}(($(entries...),))
    end
end

@generated function static_flash_jacobian(
        residual::SVector{M, <:ForwardDiff.Dual}) where M
    entries = [:(residual[$row].partials[$col])
        for col in 1:M for row in 1:M]
    return :(SMatrix{$M, $M}(($(entries...),)))
end

@generated function static_flash_residual_values(
        residual::SVector{M, <:ForwardDiff.Dual}) where M
    entries = [:(ForwardDiff.value(residual[$i])) for i in 1:M]
    return :(SVector{$M}(($(entries...),)))
end

@generated function static_flash_K(u::SVector{M}) where M
    entries = [:(u[$i]) for i in 1:(M-1)]
    return :(SVector{$(M-1)}(($(entries...),)))
end

@generated function static_equilibrium_residual(
        eos::GenericCubicEOS{E, R, N}, cond,
        unknowns::SVector{M, F}) where {E, R, N, M, F}
    M == N + 1 || error("Expected N K-values and one vapor fraction")
    k_entries = [:(unknowns[$i]) for i in 1:N]
    x_entries = [:(liquid_mole_fraction(z[$i], K[$i], V)) for i in 1:N]
    y_entries = [:(vapor_mole_fraction(x[$i], K[$i])) for i in 1:N]
    residual_entries = [:(log(K[$i]) + lnphi_v[$i] - lnphi_l[$i])
        for i in 1:N]
    rr = :(zero(F))
    for i in 1:N
        rr = :($rr + z[$i]*(K[$i] - one(F))/
            (one(F) + V*(K[$i] - one(F))))
    end
    push!(residual_entries, rr)
    return quote
        V = unknowns[$M]
        K = SVector{$N, F}(($(k_entries...),))
        z = cond.z
        x = SVector{$N, F}(($(x_entries...),))
        y = SVector{$N, F}(($(y_entries...),))
        liquid = (p = cond.p, T = cond.T, z = x, phase = Val(:liquid))
        vapor = (p = cond.p, T = cond.T, z = y, phase = Val(:vapor))
        forces = static_force_coefficients(eos, cond, F)
        Z_l, scalars_l = prep(eos, liquid, forces)
        Z_v, scalars_v = prep(eos, vapor, forces)
        lnphi_l = static_fugacity_coefficients(
            eos, liquid, Z_l, forces, scalars_l, F)
        lnphi_v = static_fugacity_coefficients(
            eos, vapor, Z_v, forces, scalars_v, F)
        SVector{$M, F}(($(residual_entries...),))
    end
end

end # module ImplicitFlashDerivatives
