"""
    stability_2ph(eos, c, [K])

Determine if mixture is single-phase stable under conditions `c`.

This is done using a version of Michelsen's stability test.

Reference: [The isothermal flash problem. Part I. Stability](https://doi.org/10.1016/0378-3812(82)85001-2)
"""
function stability_2ph(eos, c, K = initial_guess_K(eos, c); kwarg...)
    return stability_2ph(eos, c, K, FlashConfig(); kwarg...)
end

function stability_2ph(eos, c, K, config::FlashConfig; kwarg...)
    storage = flash_storage(eos, c, SSIFlash(), config)
    stability_2ph!(storage, K, eos, c, config; kwarg...)
end

function stability_2ph!(storage, K, eos, c; kwarg...)
    return stability_2ph!(storage, K, eos, c, FlashConfig(); kwarg...)
end

function stability_2ph!(storage, K, eos, c, config::FlashConfig;
        verbose::Bool = false,
        extra_out::Bool = false,
        check_vapor::Bool = true,
        check_liquid::Bool = true,
        kwarg...
    )
    forces = storage.forces
    f_z = storage.buffer1
    f_xy = storage.buffer2
    x, y = storage.x, storage.y
    z, p, T = c.z, c.p, c.T
    liquid_phase = phase_value(config, Val(:liquid))
    vapor_phase = phase_value(config, Val(:vapor))
    liquid = (p = p, T = T, z = x, phase = liquid_phase)
    vapor = (p = p, T = T, z = y, phase = vapor_phase)
    current_as_liquid = (p = p, T = T, z = z, phase = liquid_phase)
    current_as_vapor = (p = p, T = T, z = z, phase = vapor_phase)
    mixture_fugacities!(f_z, eos, current_as_vapor, forces)
    if check_vapor
        wilson_estimate!(K, eos, p, T)
        v = michelsen_test!(vapor, f_z, f_xy, vapor.z, z, K, eos, c, forces, Val(true), config; kwarg...)
    else
        v = (true, true, 0)
    end
    stable_vapor, trivial_vapor, i_v = v
    if check_liquid
        if forces_per_phase(eos)
            # Need to recalculate fugacities for the liquid phase if the flash
            # uses e.g. different bic coefficients for each phase. Otherwise,
            # these are already ok.
            mixture_fugacities!(f_z, eos, current_as_liquid, forces)
        end
        wilson_estimate!(K, eos, p, T)
        l = michelsen_test!(liquid, f_z, f_xy, liquid.z, z, K, eos, c, forces, Val(false), config; kwarg...)
    else
        l = (true, true, 0)
    end
    stable_liquid, trivial_liquid, i_l = l
    report = StabilityReport(
        stable_liquid = stable_liquid,
        trivial_liquid = trivial_liquid,
        stable_vapor = stable_vapor,
        trivial_vapor = trivial_vapor
    )
    stable = report.stable
    if !stable
        @. K = y/x
    end
    if print_output(config) && verbose
        @info "Stability done. Iterations:\nV: $i_v\nL: $i_l" stable_vapor stable_liquid stable
    end
    if extra_out
        out = (stable, report)
    else
        out = stable
    end
    return out
end

f_ratio(f_z, f_xy, ::Val{true}) = f_z/f_xy
f_ratio(f_z, f_xy, ::Val{false}) = f_xy/f_z
xy_value(z, K, ::Val{true}) = z*K
xy_value(z, K, ::Val{false}) = z/K

stability_phase(::Val{true}) = Val(:vapor)
stability_phase(::Val{false}) = Val(:liquid)

@generated function stability_xy(z::SVector{N, F}, K::SVector{N, F}, phase) where {N, F}
    values = [:(xy_value(z[$i], K[$i], phase)) for i in 1:N]
    return :(SVector{N, F}(($(values...),)))
end

@generated function static_scale(v::SVector{N, F}, scale::F) where {N, F}
    values = [:(v[$i]/scale) for i in 1:N]
    return :(SVector{N, F}(($(values...),)))
end

@generated function stability_ratios(f_z::SVector{N, F}, f_xy::SVector{N, F},
        scale::F, phase) where {N, F}
    values = [:(f_ratio(f_z[$i], scale*f_xy[$i], phase)) for i in 1:N]
    return :(SVector{N, F}(($(values...),)))
end

@generated function static_multiply(a::SVector{N, F}, b::SVector{N, F}) where {N, F}
    values = [:(a[$i]*b[$i]) for i in 1:N]
    return :(SVector{N, F}(($(values...),)))
end

@generated function static_divide(a::SVector{N, F}, b::SVector{N, F}) where {N, F}
    values = [:(a[$i]/b[$i]) for i in 1:N]
    return :(SVector{N, F}(($(values...),)))
end

@inline function michelsen_test_stack(f_z, z::SVector{N, F}, K::SVector{N, F},
        eos, cond, forces, inside_is_vapor;
        tol_equil = 1e-10,
        tol_trivial = tol_equil,
        tol_sat = tol_trivial,
        maxiter = 1000
    ) where {N, F}
    trivial = false
    S = one(F)
    iter = 0
    xy = zero(SVector{N, F})
    while true
        iter += 1
        unnormalized = stability_xy(z, K, inside_is_vapor)
        S = zero(F)
        @inbounds for i in 1:N
            S += unnormalized[i]
        end
        xy = static_scale(unnormalized, S)
        inside = (p = cond.p, T = cond.T, z = xy,
            phase = stability_phase(inside_is_vapor))
        f_xy = static_fugacities(eos, inside, forces, F)
        ratios = stability_ratios(f_z, f_xy, S, inside_is_vapor)
        K = static_multiply(K, ratios)
        R_norm = zero(F)
        K_norm = zero(F)
        @inbounds for i in 1:N
            R_norm += (ratios[i] - one(F))^2
            K_norm += log(K[i])^2
        end
        trivial = K_norm < tol_trivial
        converged = R_norm < tol_equil
        if trivial || converged
            break
        elseif iter == maxiter
            # Match the checked host path: failed stability iterations are
            # conservatively treated as a trivial, stable solution.
            trivial = true
            break
        end
    end
    stable = trivial || S <= one(F) + tol_sat
    return stable, trivial, iter, K, xy
end

@inline function stability_2ph_stack(K::SVector{N, F}, eos, cond, forces;
        check_vapor::Bool = true,
        check_liquid::Bool = true,
        kwarg...
    ) where {N, F}
    vapor_phase = (p = cond.p, T = cond.T, z = cond.z, phase = Val(:vapor))
    f_z_vapor = static_fugacities(eos, vapor_phase, forces, F)
    K_wilson = initial_guess_K_static(eos, cond, FlashConfig{false, false}())
    if check_vapor
        stable_vapor, trivial_vapor, i_v, K_vapor, y = michelsen_test_stack(
            f_z_vapor, cond.z, K_wilson, eos, cond, forces, Val(true); kwarg...)
    else
        stable_vapor, trivial_vapor, i_v, K_vapor, y = true, true, 0, K_wilson, cond.z
    end
    if check_liquid
        liquid_phase = (p = cond.p, T = cond.T, z = cond.z, phase = Val(:liquid))
        f_z_liquid = forces_per_phase(eos) ?
            static_fugacities(eos, liquid_phase, forces, F) : f_z_vapor
        stable_liquid, trivial_liquid, i_l, K_liquid, x = michelsen_test_stack(
            f_z_liquid, cond.z, K_wilson, eos, cond, forces, Val(false); kwarg...)
    else
        stable_liquid, trivial_liquid, i_l, K_liquid, x = true, true, 0, K_wilson, cond.z
    end
    report = StabilityReport(stable_liquid, trivial_liquid,
        stable_vapor, trivial_vapor)
    K_out = report.stable ? K_liquid : static_divide(y, x)
    return report.stable, report, K_out
end

"""
    stability_2ph!(storage, eos, c, [K])

In-place version of [`stability_2ph`](@ref). `storage` should be allocated by `flash_storage`.
"""
function michelsen_test!(c_inside, f_z, f_xy, xy, z, K, eos, cond, forces, inside_is_vapor; kwarg...)
    return michelsen_test!(c_inside, f_z, f_xy, xy, z, K, eos, cond, forces,
        inside_is_vapor, FlashConfig(); kwarg...)
end

function michelsen_test!(c_inside, f_z, f_xy, xy, z, K, eos, cond, forces,
        inside_is_vapor, config::FlashConfig;
        tol_equil = 1e-10,
        tol_trivial = tol_equil,
        tol_sat = tol_trivial,
        maxiter = 1000
    )
    trivial = false
    S = 1.0
    iter = 0
    done = false
    while !done
        iter += 1
        S = 0.0
        @inbounds for c in eachindex(xy)
            xy_i = xy_value(z[c], K[c], inside_is_vapor)
            xy[c] = xy_i
            S += xy_i
        end
        @. xy /= S
        mixture_fugacities!(f_xy, eos, c_inside, forces)

        R_norm = 0.0
        K_norm = 0.0
        @inbounds for c in eachindex(K)
            R = f_ratio(f_z[c], S*f_xy[c], inside_is_vapor)
            K[c] *= R

            R_norm += (R-1)^2
            K_norm += log(K[c])^2
        end
        # Two convergence criteria:
        # - Approaching trivial solution (K-values are all 1)
        # - Equilibrium for a small amount of the "other" phase,
        #   the single-phase conditions are not stable.
        trivial = K_norm < tol_trivial
        converged = R_norm < tol_equil

        # Termination of loop
        ok = trivial || converged
        done = ok || iter == maxiter
        if done && !ok
            trivial = true
            if print_output(config)
                @warn "Stability test failed to converge in $maxiter iterations. Assuming stability." cond xy K_norm R_norm K
            end
        end
    end
    stable = trivial || S <= 1.0 + tol_sat
    return (stable, trivial, iter)
end
