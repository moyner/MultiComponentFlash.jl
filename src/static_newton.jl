# Phase-amount Gibbs minimization for the short reuse attempt described in
# Belkadi et al. (2011), SPE 142132, Section 2.1. The caller handles fallback.
@inline function static_newton_flash(eos, cond, K, V, tolerance, maxiter, dMax)
    if isinf(V) || !(
            all(k -> isfinite(k) && k > 0, K) &&
                all(z -> isfinite(z) && z > 0, cond.z)
        )
        return V, K, false, 0
    end
    # Always restore the material-balance root before the first EOS evaluation.
    # A caller's finite V may belong to a previous composition or set of K values.
    initial_V = solve_rachford_rice(
        K, cond.z, V;
        tol = max(1.0e-12, eps(eltype(K)))
    )
    if !(isfinite(initial_V) && zero(initial_V) < initial_V < one(initial_V))
        return initial_V, K, false, 0
    end
    liquid = map(cond.z, K) do zi, ki
        (one(initial_V) - initial_V) * liquid_mole_fraction(zi, ki, initial_V)
    end
    vapor = map(cond.z, K) do zi, ki
        initial_V * vapor_mole_fraction(liquid_mole_fraction(zi, ki, initial_V), ki)
    end
    liquid, vapor, converged, iterations = static_gibbs_newton(
        eos, cond, liquid, vapor, tolerance; maxiter = maxiter, dMax = dMax
    )
    total_liquid, total_vapor = sum(liquid), sum(vapor)
    next_K = (vapor / total_vapor) ./ (liquid / total_liquid)
    return total_vapor / (total_liquid + total_vapor), next_K, converged, iterations
end

@inline function static_newton_phase_potentials(eos, cond, amounts, phase, forces)
    composition = amounts / sum(amounts)
    phase_cond = (p = cond.p, T = cond.T, z = composition, phase = phase)
    Z, scalars = prep(eos, phase_cond, forces)
    lnphi = static_fugacity_coefficients(
        eos, phase_cond, Z, forces, scalars, eltype(amounts)
    )
    return log.(composition) + lnphi
end

@inline function static_newton_system(
        eos, cond, liquid::SVector{N, F},
        vapor::SVector{N, F}
    ) where {N, F}
    L, V = sum(liquid), sum(vapor)
    K = (vapor / V) ./ (liquid / L)
    residual, jacobian = ImplicitFlashDerivatives.static_flash_system(
        eos, cond, K, V / (L + V)
    )
    gradient = SVector{N, F}(ntuple(i -> residual[i], Val(N)))
    scale = sqrt.(liquid .* (vapor ./ cond.z))
    # Reuse the (K, beta) Jacobian of the two-phase implicit derivatives.
    # For v_i = z_i - l_i, d(log K_i)/dv_j = delta_ij*(1/v_i + 1/l_i)
    # - 1/V - 1/L, and d(beta)/dv_j = 1/(L + V). Apply this chain rule
    # directly to the scaled Hessian to avoid large trace-component reciprocals.
    row_correction = SVector{N, F}(
        ntuple(Val(N)) do i
            total = sum(j -> jacobian[i, j] * K[j], 1:N)
            jacobian[i, N + 1] / (L + V) - (inv(L) + inv(V)) * total
        end
    )
    hessian = SMatrix{N, N, F}(
        ntuple(Val(N * N)) do index
            row = mod1(index, N)
            col = (index - 1) ÷ N + 1
            (scale[row] / scale[col]) * jacobian[row, col] * K[col] +
                scale[row] * scale[col] * row_correction[row]
        end
    )
    return gradient, (hessian + transpose(hessian)) / 2, scale
end

@inline static_positive_phase_amounts(amounts) = all(x -> isfinite(x) && x > 0, amounts)

@inline function static_gibbs_newton(
        eos, cond, liquid::SVector{N, F},
        vapor::SVector{N, F}, tolerance; maxiter = 3, dMax = 1.0
    ) where {N, F}
    forces = static_force_coefficients(eos, cond, F)
    # StaticArrays' small Cholesky is allocation-free and accelerator compatible.
    # Larger systems retain the ordinary flash instead of a dynamic factorization.
    N <= 24 || return liquid, vapor, false, 0
    for iteration in 0:maxiter
        if !(static_positive_phase_amounts(liquid) && static_positive_phase_amounts(vapor))
            return liquid, vapor, false, iteration
        end
        gradient, H, scale = static_newton_system(eos, cond, liquid, vapor)
        all(isfinite, gradient) && all(isfinite, H) ||
            return liquid, vapor, false, iteration
        factor = cholesky(Symmetric(H); check = false)
        if !issuccess(factor) || !all(
                i -> isfinite(factor.U[i, i]) &&
                    factor.U[i, i] > 0, 1:N
            )
            return liquid, vapor, false, iteration
        end
        if maximum(abs, gradient) <= tolerance
            return liquid, vapor, true, iteration
        elseif iteration == maxiter
            return liquid, vapor, false, iteration
        end
        increment = scale .* (factor \ (-scale .* gradient))
        max_increment = maximum(abs, increment)
        if max_increment > F(dMax)
            increment *= F(dMax) / max_increment
        end
        next_liquid = liquid - increment
        next_vapor = vapor + increment
        if !(static_positive_phase_amounts(next_liquid) && static_positive_phase_amounts(next_vapor))
            return liquid, vapor, false, iteration + 1
        end
        next_mu_l = static_newton_phase_potentials(eos, cond, next_liquid, Val(:liquid), forces)
        next_mu_v = static_newton_phase_potentials(eos, cond, next_vapor, Val(:vapor), forces)
        mu_l = static_newton_phase_potentials(eos, cond, liquid, Val(:liquid), forces)
        mu_v = static_newton_phase_potentials(eos, cond, vapor, Val(:vapor), forces)
        gibbs = sum(liquid .* mu_l) + sum(vapor .* mu_v)
        next_gibbs = sum(next_liquid .* next_mu_l) + sum(next_vapor .* next_mu_v)
        if !(isfinite(next_gibbs) && next_gibbs <= gibbs + 8eps(F) * max(1, abs(gibbs)))
            return liquid, vapor, false, iteration + 1
        end
        liquid, vapor = next_liquid, next_vapor
    end
    return liquid, vapor, false, maxiter
end
