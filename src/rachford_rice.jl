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
