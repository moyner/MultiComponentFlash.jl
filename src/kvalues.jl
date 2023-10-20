function wilson_estimate(c::MolecularProperty, p, T)
    return wilson_estimate(p, T, c.ω, c.p_c, c.T_c)
end

"""
    wilson_estimate(p, T, ω, p_c, T_c)

Estimate K-values for a given acentric factor ω and pressure and temperature at current and critical conditions.

Reference [Vapor-Liquid Equilibrium. XI. A New Expression for the Excess Free Energy of Mixing by GM Wilson](https://doi.org/10.1021/ja01056a002)
"""
function wilson_estimate(p::R, T::R, ω::R, p_c::R, T_c::R) where R<:Real
    K = exp(5.37*(1.0 + ω)*(1.0 - T_c/T))*(p_c/p)
    return K::R
end

"""
    wilson_estimate!(K, properties, p, T)

Update a vector K in-place with K-values from `wilson_estimate`.
"""
function wilson_estimate!(K::AbstractVector{R}, props, p::R, T::R) where R<:Real
    @inbounds for i in eachindex(K)
        K[i] = wilson_estimate(props[i], p, T)::R
    end
end

function wilson_estimate!(K::AbstractVector{R}, props::MultiComponentMixture, p::R, T::R) where R<:Real
    return wilson_estimate!(K, props.properties, p, T)
end

function wilson_estimate!(K::AbstractVector{R}, props::AbstractCubicEOS, p::R, T::R) where R<:Real
    return wilson_estimate!(K, props.mixture, p, T)
end

#used as a hook for extensions
function wilson_estimate!(K::AbstractVector{R}, props, p::R, T::R, storage) where R<:Real
    return wilson_estimate!(K, props, p, T)
end
"""
    wilson_estimate(properties, p, T)

Create vector of K-values that holds the `wilson_estimate` for each species.
"""
function wilson_estimate(mixture::MultiComponentMixture, p, T)
    K = zeros(number_of_components(mixture))
    wilson_estimate!(K, mixture.properties, p, T)
    return K
end

"""
    initial_guess_K(eos, cond)

Produce a plausible initial guess for K values for `eos` under current `cond`.
"""
function initial_guess_K(eos, cond)
    K = zeros(number_of_components(eos))
    initial_guess_K!(K, eos, cond)
    return K
end

"""
    initial_guess_K(eos, cond)

In-place version of `initial_guess_K`.
"""
initial_guess_K!(K, eos, cond) = wilson_estimate!(K, eos.mixture.properties, cond.p, cond.T)
