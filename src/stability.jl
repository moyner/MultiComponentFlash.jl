"""
    stability_2ph(eos, c, [K])

Determine if mixture is single-phase stable under conditions `c`.

This is done using a version of Michelsen's stability test.

Reference: [The isothermal flash problem. Part I. Stability](https://doi.org/10.1016/0378-3812(82)85001-2)
"""
function stability_2ph(eos, c, K = initial_guess_K(eos, c); kwarg...)
    storage = flash_storage(eos, c)
    stability_2ph!(storage, K, eos, c)
end

function stability_2ph!(storage, K, eos, c; verbose = false, kwarg...)
    forces = storage.forces
    f_z = storage.buffer1
    f_xy = storage.buffer2
    x, y = storage.x, storage.y
    z, p, T = c.z, c.p, c.T
    liquid = (p = p, T = T, z = x)
    vapor = (p = p, T = T, z = y)
    # Update fugacities for current conditions used in both tests
    mixture_fugacities!(f_z, eos, c, forces)
    wilson_estimate!(K, eos, p, T)
    stable_vapor, i_v = michelsen_test!(vapor, f_z, f_xy, vapor.z, z, K, eos, c, forces, Val(true); kwarg...)
    wilson_estimate!(K, eos, p, T)
    stable_liquid, i_l = michelsen_test!(liquid, f_z, f_xy, liquid.z, z, K, eos, c, forces, Val(false); kwarg...)
    stable = stable_vapor && stable_liquid
    if !stable
        @. K = y/x
    end
    if verbose
        @info "Stability done. Iterations:\nV: $i_v\nL: $i_l" stable_vapor stable_liquid stable
    end
    return stable
end

f_ratio(f_z, f_xy, ::Val{true}) = f_z/f_xy
f_ratio(f_z, f_xy, ::Val{false}) = f_xy/f_z
xy_value(z, K, ::Val{true}) = z*K
xy_value(z, K, ::Val{false}) = z/K

"""
    stability_2ph!(storage, eos, c, [K])

In-place version of [`stability_2ph`](@ref). `storage` should be allocated by `flash_storage`.
"""
function michelsen_test!(c_inside, f_z, f_xy, xy, z, K, eos, cond, forces, inside_is_vapor;
    tol_equil = 1e-10, tol_trivial = 1e-5, maxiter = 1000)
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
            @warn "Stability test failed to converge in $maxiter iterations. Assuming stability." cond xy K_norm R_norm K
        end
    end
    stable = trivial || S <= 1 + tol_trivial
    return (stable, iter)
end
