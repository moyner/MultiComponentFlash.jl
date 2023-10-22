"""
    solve_rachford_rice(K, z, [V]; <keyword arguments>)

Compute vapor mole fraction `V` for given equilibrium constants `K` and mole fractions `z`.

# Arguments
`K` - Equal length to `z`, containing the equilibrium constants for each component.
`z` - Mole fractions. Should sum up to unity.

# Keyword arguments
`tol = 1e-12` - Tolerance for solve.
`maxiter` - Maximum number of iterations
`ad` - Use automatic differentiation (ForwardDiff) instead of analytical gradient.

# Examples
```julia-repl
julia> solve_rachford_rice([0.5, 1.5], [0.3, 0.7])
0.8000000000000002
```
"""
function solve_rachford_rice(K, z, V = NaN; tol=1e-12, maxiter=1000, ad = false)
    # Lowest bound for solution
    V_min = 1/(1 - maximum(K))
    # Largest bound for solution
    V_max = 1/(1 - minimum(K))
    if isnan(V)
        V = (V_min + V_max)/2
    end
    i = 1
    n = length(z)
    if n == 2 #exact solution, linear
        z1,z2 = z
        k1,k2 = K
        b1,b2 = 1/(1-k1),1/(1-k2)
        a1 = (z1*b2 + z2*b1)
        a0 = z1+z2
        return a1/a0
    elseif n == 3 #exact solution, quadratic
        z1,z2,z3 = z
        k1,k2,k3 = K
        b1,b2,b3 = 1/(1-k1),1/(1-k2),1/(1-k3)
        a2 =(z1 + z2 + z3)
        a1 = -b1*(z2 + z3) - b2*(z1 + z3) - b3*(z1 + z2)
        a0 = b1*b2*z3 + b1*b3*z2 + b2*b3*z1
        Δ = a1*a1 - 4*a0*a2
        inva2 = 1/(2*a2)
        Δ2 = sqrt(Δ)*inva2
        x = -a1*inva2
        β1 = x-Δ2
        β2 = x+Δ2
        βmax = max(β1,β2)
        βmin = min(β1,β2)
        if 0 < β1 < 1
            return β1
        elseif 0 < β2 < 1
            return β2
        elseif isfinite(β1+β2)
            kmin = min(k1,k2,k3)
            kmax = max(k1,k2,k3)
            kmin > 1 && return βmax
            kmax < 1 && return βmin
        else
            return zero(Δ2)/zero(Δ2)
        end
    end
    while i < maxiter
        if ad
            f(V) = objectiveRR(V, K, z)
            r = f(V)
            Jo = ForwardDiff.derivative(f, V)
            V = V - Jo\r
        else
            denum = 0.0
            r = 0.0
            @inbounds for i = 1:n
                ΔK = K[i]-1
                a = z[i]*ΔK
                b = (1 + V*ΔK)
                r += a/b
                denum += z[i]*ΔK^2/b^2
            end
            V = V + r/denum
        end
        if abs(r) < tol
            break
        end
        if V < V_min || V > V_max
            obj = (v) -> objectiveRR(v, K, z)
            V = find_zero(obj, (V_min, V_max), Bisection())
            break
        end
        i += 1
    end
    return V
end

function objectiveRR(V, K, z)
    eq = 0.0
    for (i, k) in enumerate(K)
        @inbounds eq = eq + ((k - 1.0)*z[i])/(1.0 + V*(k - 1.0))
    end
    return eq
end
