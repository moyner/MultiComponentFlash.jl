"""
    solve_rachford_rice(K, z, [V]; <keyword arguments>)

Compute the physical vapor mole fraction for equilibrium constants `K` and
overall mole fractions `z`. Return `0` for liquid-only conditions and `1` for
vapor-only conditions. Use [`solve_rachford_rice_unconstrained`](@ref) when a
negative flash requires a root outside `[0, 1]`.

`V` is an optional initial guess for an interior two-phase root. Solver keyword
arguments are forwarded to `solve_rachford_rice_unconstrained` in that case.
"""
@inline function solve_rachford_rice(K, z, V = NaN;
        tol = 1e-12, maxiter = 1000, ad = false, analytical = true,
        verbose = false)
    # RR decreases on [0, 1] for positive K. Its endpoint signs distinguish
    # a physical split from a single-phase condition.
    r_liquid = r_vapor = zero(K[1]*z[1])
    @inbounds for i in eachindex(z)
        K_i, z_i = K[i], z[i]
        if !isfinite(K_i) || K_i <= zero(K_i) ||
                !isfinite(z_i) || z_i < zero(z_i)
            return oftype(r_liquid, NaN)
        end
        delta_K = K_i - one(K_i)
        r_liquid += z_i*delta_K
        r_vapor += z_i*delta_K/K_i
    end
    if r_liquid <= zero(r_liquid)
        return zero(r_liquid)
    elseif r_vapor >= zero(r_vapor)
        return one(r_vapor)
    end
    return solve_rachford_rice_unconstrained(K, z, V;
        tol = tol, maxiter = maxiter, ad = ad, analytical = analytical,
        verbose = verbose)
end

# Retain the previous internal spelling for callers that used it directly.
@inline physical_vapor_fraction(K, z, V = NaN) = solve_rachford_rice(K, z, V)

"""
    solve_rachford_rice_unconstrained(K, z, [V]; <keyword arguments>)

Compute a negative-flash vapor fraction for given equilibrium constants `K`
and mole fractions `z`. The root may lie outside `[0, 1]`, but both phase
compositions must remain nonnegative. Return `NaN` when there is no such root.

# Arguments
`K` - Equal length to `z`, containing the equilibrium constants for each component.
`z` - Mole fractions. Should sum up to unity.
`V` - Optional initial guess. `NaN` or `Inf` selects the midpoint of the admissible interval for iterative solves.

# Keyword arguments

- `tol = 1e-12`: Tolerance for solve.
- `maxiter=1000`: Maximum number of iterations
- `ad=false`: Use automatic differentiation (ForwardDiff) instead of analytical gradient.
- `analytical=true`: Use analytical solutions for 2 and 3 components.

# Examples
```julia-repl
julia> solve_rachford_rice_unconstrained([0.5, 1.5], [0.3, 0.7])
0.8000000000000002
```
"""
function solve_rachford_rice_unconstrained(K, z, V = NaN; tol = 1e-12, maxiter = 1000,
        ad = false, analytical = true, verbose = false)
    V_lo, V_hi = positive_rachford_rice_bounds(K, z)
    V_lo < V_hi || return oftype(K[1], NaN)
    if analytical
        n = length(z)
        if n == 2
            root = rachford_rice_analytic_2(K, z, V_lo, V_hi)
            if isfinite(root) &&
                    rachford_rice_balance_error(root, objectiveRR(root, K, z)) <= tol
                return root
            end
        elseif n == 3
            root = rachford_rice_analytic_3(K, z, V_lo, V_hi)
            if isfinite(root) &&
                    rachford_rice_balance_error(root, objectiveRR(root, K, z)) <= tol
                return root
            end
        end
    end
    return solve_rachford_rice_bounded(K, z, V, V_lo, V_hi;
        tol = tol, maxiter = maxiter, ad = ad, verbose = verbose)
end

@inline function rachford_rice_balance_error(V, residual)
    # RR = sum(y) - sum(x). For normalized z, the two normalization errors
    # are -V*RR and (1 - V)*RR. A small unscaled RR residual is insufficient
    # when a negative flash has a very large |V|.
    return abs(residual)*max(one(V), abs(V), abs(one(V) - V))
end

@inline function positive_rachford_rice_bounds(K, z)
    # Whitson and Michelsen, Fluid Phase Equilibria 53 (1989), 51-71:
    # the negative-flash window is the intersection of
    # 1 + V*(K_i - 1) > 0 for every present component. Other intervals
    # between poles can contain roots, but give negative phase compositions.
    K_min = K_max = one(K[1])
    has_below = has_above = false
    invalid = oftype(K[1], NaN)
    @inbounds for i in eachindex(z)
        z_i, K_i = z[i], K[i]
        if !isfinite(z_i) || z_i < zero(z_i) ||
                (z_i > zero(z_i) && (!isfinite(K_i) || K_i <= zero(K_i)))
            return invalid, invalid
        end
        if z_i > zero(z_i)
            if K_i < one(K_i)
                K_min = min(K_min, K_i)
                has_below = true
            elseif K_i > one(K_i)
                K_max = max(K_max, K_i)
                has_above = true
            end
        end
    end
    has_below && has_above || return invalid, invalid
    return inv(one(K_max) - K_max), inv(one(K_min) - K_min)
end

@inline function rachford_rice_analytic_2(K, z, V_lo, V_hi)
    k1, k2 = K
    if k1 == one(k1) || k2 == one(k2)
        return oftype(k1, NaN)
    end
    z1, z2 = z
    b1, b2 = inv(one(k1) - k1), inv(one(k2) - k2)
    root = (z1*b2 + z2*b1)/(z1 + z2)
    return V_lo < root < V_hi ? root : oftype(root, NaN)
end

@inline function rachford_rice_analytic_3(K, z, V_lo, V_hi)
    k1, k2, k3 = K
    if k1 == one(k1) || k2 == one(k2) || k3 == one(k3)
        return oftype(k1, NaN)
    end
    z1, z2, z3 = z
    b1, b2, b3 = inv(one(k1) - k1), inv(one(k2) - k2), inv(one(k3) - k3)
    a2 = z1 + z2 + z3
    a1 = -b1*(z2 + z3) - b2*(z1 + z3) - b3*(z1 + z2)
    a0 = b1*b2*z3 + b1*b3*z2 + b2*b3*z1
    discriminant = a1*a1 - 4*a0*a2
    if discriminant >= zero(discriminant)
        inv_2a2 = inv(2*a2)
        offset = sqrt(discriminant)*inv_2a2
        center = -a1*inv_2a2
        root1, root2 = center - offset, center + offset
        if V_lo < root1 < V_hi
            return root1
        elseif V_lo < root2 < V_hi
            return root2
        end
    end
    return oftype(k1, NaN)
end

@inline function solve_rachford_rice_bounded(K, z, V, V_lo, V_hi;
        tol = 1e-12, maxiter = 1000, ad = false, verbose = false)
    # A supplied guess may be on a pole or outside the window after an SSI
    # K-update. Never evaluate the Rachford-Rice function at that guess.
    if !(V_lo < V < V_hi) || !isfinite(V)
        V = V_lo/2 + V_hi/2
    end
    V += zero(V_lo)
    best_V = V
    best_error = oftype(V, Inf)
    verbose && println("Solving Rachford-Rice in ($V_lo, $V_hi) from $V")
    for iteration in 1:maxiter
        residual = zero(V)
        denominator = zero(V)
        @inbounds for i in eachindex(z)
            z_i = z[i]
            iszero(z_i) && continue
            delta_K = K[i] - one(K[i])
            term = muladd(V, delta_K, one(V))
            if !(term > zero(term)) || !isfinite(term)
                return oftype(V, NaN)
            end
            residual += z_i*delta_K/term
            denominator += z_i*delta_K^2/term^2
        end
        if !isfinite(residual) || !isfinite(denominator)
            return oftype(V, NaN)
        end
        balance_error = rachford_rice_balance_error(V, residual)
        if balance_error < best_error
            best_V, best_error = V, balance_error
        end
        balance_error <= tol && return V
        denominator > zero(denominator) || return oftype(V, NaN)
        if ad
            denominator = -ForwardDiff.derivative(v -> objectiveRR(v, K, z), V)
        end
        if residual > zero(residual)
            V_lo = V
        else
            V_hi = V
        end
        V_next = V + residual/denominator
        if !(V_lo < V_next < V_hi) || !isfinite(V_next) || V_next == V
            V_next = V_lo/2 + V_hi/2
        end
        if !(V_lo < V_next < V_hi) || V_next == V
            # A root arbitrarily close to a pole can exhaust Float64 spacing
            # before meeting the requested residual tolerance. Use the best
            # representable point only if both phase sums remain accurate.
            return best_error <= max(tol, 1e-8) ? best_V : oftype(V, NaN)
        end
        verbose && println("#$iteration V = $V_next, residual = $residual")
        V = V_next
    end
    return oftype(V, NaN)
end

function objectiveRR(V, K, z)
    eq = 0.0
    for (i, k) in enumerate(K)
        iszero(z[i]) && continue
        @inbounds eq = eq + ((k - 1.0)*z[i])/muladd(V, k - 1.0, 1.0)
    end
    return eq
end

function objectiveRR_dV(V, K, z)
    RR_dv = 0.0
    for i in eachindex(K)
        z_i = z[i]
        K_i = K[i]
        RR_dv -= (z_i*(K_i - 1)^2)/(1+V*(K_i-1))^2
    end
    return RR_dv
end

function objectiveRR_dK(V, K, z, i)
    return z[i]/(1.0 + V*(K[i]-1))^2
end

function objectiveRR_dz(V, K, z, i)
    return (K[i] - 1.0)/(1.0 + V*(K[i]-1))
end
