@inline reduced_pT(mix, c, i) = (reduced_pressure(mix, c, i), reduced_temperature(mix, c, i))
@inline reduced_pT(eos::AbstractEOS, c, i) = (reduced_pressure(eos.mixture, c, i), reduced_temperature(eos.mixture, c, i))

"""
    number_of_components(eos)

Return number of components for the underlying mixture of the EOS.
"""
number_of_components(e::AbstractEOS) = number_of_components(e.mixture)

forces_per_phase(eos::GenericCubicEOS) = false

function get_phase(cond)
    return phase_symbol(get(cond, :phase, :unknown))
end

@inline phase_symbol(phase::Symbol) = phase
@inline phase_symbol(::Val{phase}) where phase = phase

function set_phase(cond, phase::Symbol, throw::Bool = false)
    if throw && haskey(cond, :phase) && cond.phase != :unknown
        throw(ArgumentError("Phase state already set to $(cond.phase), cannot change to $phase."))
    end
    return (p = cond.p, T = cond.T, z = cond.z, phase = phase)
end

function static_coefficients(::AbstractGeneralizedCubic)
    return (0.4274802327, 0.08664035, 0.0, 1.0)
end

function weight_bi(eos::GenericCubicEOS{T, R}, cond, i) where {T<:AbstractGeneralizedCubic, R}
    return eos.ω_b
end

function weight_ai(eos::GenericCubicEOS{T, R}, cond, i) where {T<:AbstractGeneralizedCubic, R}
    mix = eos.mixture
    T_r = reduced_temperature(mix, cond, i)
    return eos.ω_a*T_r^(-1/2)
end

# Peng Robinson and variants
include("eos/peng_robinson.jl")
include("eos/peng_robinson_corrected.jl")
include("eos/soreide_whitson.jl")
# Soave Redlich-Kwong
include("eos/soave_redlich_kwong.jl")
# Zudkevitch-Joffe
include("eos/zudkevitch_joffe.jl")

# Generic part
Base.@propagate_inbounds function binary_interaction(eos::AbstractEOS, i::Int, j::Int, cond)
    return binary_interaction(eos.mixture, i, j)
end

Base.@propagate_inbounds function binary_interaction(mixture::MultiComponentMixture{R}, i, j) where {R}
    return binary_interaction(mixture.binary_interaction, i, j)::R
end

function binary_interaction(::Nothing, i, j)
    return 0.0
end

Base.@propagate_inbounds function binary_interaction(B::AbstractMatrix, i, j)
    return B[i, j]
end

"""
    mixture_compressibility_factor(eos, cond, [forces, scalars, phase])

Compute compressibility factor for current conditions. Forces and scalars can be 
provided if they are already known.

The compressibility factor adjusts the ideal gas law to account for non-linear behavior: ``pV = nRTZ``
"""

function mixture_compressibility_factor(
        eos::AbstractCubicEOS,
        cond,
        forces = force_coefficients(eos, cond),
        scalars = force_scalars(eos, cond, forces),
        phase = missing
    )
    if !ismissing(phase)
        cond = set_phase(cond, phase)
    end
    poly = eos_polynomial(eos, forces, scalars)
    roots = solve_roots(eos, poly)
    r = pick_root(eos, roots, cond, forces, scalars)
    return r
end

minimum_allowable_root(eos::AbstractCubicEOS, forces, scalars) = scalars.B
minimum_allowable_root(eos, forces, scalars) = 1e-16

@inline function pick_root(eos, roots::Real, cond, forces, scalars)
    return roots
end

@inline function root_bounds(roots, minimum_root)
    max_root = -Inf
    min_root = Inf
    for root in roots
        max_root = max(max_root, root)
        if root > minimum_root
            min_root = min(min_root, root)
        end
    end
    return min_root, max_root
end

@inline function pick_root(eos, roots, cond, forces, scalars)
    phase = get(cond, :phase, :unknown)
    return pick_root(eos, roots, cond, forces, scalars, phase)
end

@inline function pick_root(eos, roots, cond, forces, scalars, ::Val{:liquid})
    min_r, _ = root_bounds(roots, minimum_allowable_root(eos, forces, scalars))
    return min_r
end

@inline function pick_root(eos, roots, cond, forces, scalars, ::Val{:vapor})
    _, max_r = root_bounds(roots, minimum_allowable_root(eos, forces, scalars))
    return max_r
end

function pick_root(eos, roots, cond, forces, scalars, phase::Symbol)
    min_r, max_r = root_bounds(roots, minimum_allowable_root(eos, forces, scalars))
    if min_r == max_r || phase == :liquid
        return min_r
    elseif phase == :vapor
        return max_r
    end
    function Gibbs(Z)
        E = 0.0
        z = cond.z
        @inbounds for i in eachindex(z)
            ϕ = component_fugacity_coefficient(eos, cond, i, Z, forces, scalars)
            E += z[i]*ϕ
        end
        return E
    end
    return Gibbs(min_r) < Gibbs(max_r) ? min_r : max_r
end

"""
    force_coefficients(eos, cond)

Get coefficients for forces for a specific EOS (component interactions). For most cubics,
these are a set of attractive (linear and quadratic) forces and a set of linear repulsive forces.

Note that the current implementation of flash assumes that these are independent of the 
compositions themselves.

See also [`force_coefficients!`](@ref)
"""
function force_coefficients(eos::AbstractCubicEOS, cond; static_size = false)
    n = number_of_components(eos)
    eT = Base.promote_eltype(cond.p, cond.T, cond.z[1])
    if static_size
        A_ij = zero(MMatrix{n, n, eT})
        A_i = zero(MVector{n, eT})
        B_i = zero(MVector{n, eT})
    else
        A_ij = zeros(eT, n, n)
        A_i = zeros(eT, n)
        B_i = zeros(eT, n)
    end
    coeff = (A_ij = A_ij, A_i = A_i, B_i = B_i)
    update_force_coefficients!(coeff, eos, cond)
    return coeff
end

"""Construct inline force coefficients with the component count from the EOS type."""
function force_coefficients_static(eos::GenericCubicEOS{E, R, N}, cond) where {E, R, N}
    T = Base.promote_eltype(cond.p, cond.T, cond.z[1])
    return force_coefficients_static(eos, cond, T)
end

function force_coefficients_static(eos::GenericCubicEOS{E, R, N}, cond, ::Type{T}) where {E, R, N, T}
    coeff = (
        A_ij = zero(MMatrix{N, N, T}),
        A_i = zero(MVector{N, T}),
        B_i = zero(MVector{N, T})
    )
    return update_force_coefficients!(coeff, eos, cond)
end

"""Immutable force coefficients for heap-free accelerator kernels."""
@inline function force_coefficients_stack(eos::GenericCubicEOS{E, R, N}, cond, ::Type{T}) where {E, R, N, T}
    A_i_static = SVector{N, T}(ntuple(i -> A_i(eos, cond, i), Val(N)))
    B_i_static = SVector{N, T}(ntuple(i -> B_i(eos, cond, i), Val(N)))
    A_ij_static = SMatrix{N, N, T}(ntuple(Val(N*N)) do index
        i = mod1(index, N)
        j = (index - 1) ÷ N + 1
        sqrt(A_i_static[i]*A_i_static[j]) *
            (one(T) - binary_interaction(eos, i, j, cond))
    end)
    return (A_ij = A_ij_static, A_i = A_i_static, B_i = B_i_static)
end

function get_force_coefficients(forces, eos, cond)
    if forces_per_phase(eos)
        phase = get_phase(cond)
        if phase == :liquid
            return forces.liquid
        elseif phase == :vapor
            return forces.vapor
        else
            error("Forces per phase are only supported for liquid and vapor phases, not $phase.")
        end
    else
        return forces
    end
end

@inline get_force_coefficients(forces, eos::GenericCubicEOS, cond) = forces

"""
    force_coefficients!(coeff, eos, cond)

In-place update of force coefficients.

See also [`force_coefficients`](@ref)
"""
function force_coefficients!(coeff, eos::AbstractCubicEOS, cond)
    if forces_per_phase(eos)
        cond = set_phase(cond, :liquid)
        update_force_coefficients!(coeff.liquid, eos, cond)
        cond = set_phase(cond, :vapor)
        update_force_coefficients!(coeff.vapor, eos, cond)
    else
        update_force_coefficients!(coeff, eos, cond)
    end
    return coeff
end

function update_force_coefficients!(coeff, eos::AbstractCubicEOS, cond)
    update_attractive_linear!(coeff.A_i, eos, cond)
    update_attractive_quadratic!(coeff.A_ij, coeff.A_i, eos, cond)
    update_repulsive!(coeff.B_i, eos, cond)
    return coeff
end

"""
Update linear part of attractive term
"""
function update_attractive_linear!(A, eos::GenericCubicEOS, cond)
    @inbounds for i in eachindex(A)
        A[i] = A_i(eos, cond, i)
    end
end

function A_i(eos::GenericCubicEOS, cond, i)
    p_r, T_r = reduced_pT(eos, cond, i)
    ω_ai = weight_ai(eos, cond, i)
    Ai = ω_ai*p_r/(T_r^2)
    return Ai
end

"""
Update quadratic part of attractive term
"""
function update_attractive_quadratic!(A_ij, A_i, eos::AbstractCubicEOS, cond)
    N = number_of_components(eos)
    T = eltype(A_i)
    for i = 1:N
        @inbounds for j = i:N
            a = sqrt(A_i[i]*A_i[j])*(one(T) - binary_interaction(eos, i, j, cond))
            a::T
            A_ij[i, j] = a
            A_ij[j, i] = a
        end
    end
end

function update_repulsive!(B, eos::AbstractCubicEOS, cond)
    @inbounds for i in eachindex(B)
        B[i] = B_i(eos, cond, i)
    end
end

function B_i(eos::GenericCubicEOS, cond, i)
    p_r, T_r = reduced_pT(eos.mixture, cond, i)
    ω_bi = weight_bi(eos, cond, i)
    return ω_bi*p_r/T_r
end

"""
    force_scalars(eos, cond, forces)

Compute EOS specific scalars for the current conditions based on the forces.
"""
function force_scalars(eos::AbstractCubicEOS, cond, forces)
    f = get_force_coefficients(forces, eos, cond)
    return cubic_scalars(f.A_ij, f.B_i, cond.z)
end

function cubic_scalars(A_ij, Bv, z)
    T = promote_type(eltype(A_ij), eltype(Bv), eltype(z))
    A = zero(T)
    B = zero(T)
    @inbounds for i in eachindex(z)
        z_i = z[i]
        B += z_i*Bv[i]
        @inbounds for j in eachindex(z)
            A += z_i*z[j]*A_ij[i,j]
        end
    end
    return (A = A, B = B)
end


"""
Generic version (only depend on scalars, ignores forces)
"""
eos_polynomial(eos::GenericCubicEOS, forces, scalars) = cubic_polynomial(eos, scalars.A, scalars.B)

"""
(Expensive version)
"""
function eos_polynomial(eos::GenericCubicEOS, cond)
    forces = force_coefficients(eos, cond)
    A, B = force_scalars(eos, cond, forces)
    cubic_polynomial(eos, A, B)
end

"""
(Precompute version)
"""
function cubic_polynomial(eos::GenericCubicEOS, A::Real, B::Real)
    m1 = eos.m_1
    m2 = eos.m_2

    a = (m1 + m2 - 1)*B - 1
    b = A - (m1 + m2 - m1*m2)*B^2 - (m1 + m2)*B
    c = -(A*B + m1*m2*(B^2)*(B+1));
    return (a, b, c)
end

"""
Fugacity
"""
Base.@propagate_inbounds function component_fugacity_coefficient(eos::AbstractCubicEOS, cond, i, Z, forces, scalars)
    # NOTE: This returns ln(ψ), not ψ!
    m1 = eos.m_1
    m2 = eos.m_2
    f = get_force_coefficients(forces, eos, cond)
    return component_fugacity_coefficient_cubic(m1, m2, cond.z, Z, scalars.A, scalars.B, f.A_ij, f.B_i, i)
end

Base.@propagate_inbounds function component_fugacity_coefficient_cubic(m1, m2, x, Z, A, B, A_mat, B_i, i)
    Δm = m1 - m2
    T = promote_type(eltype(x), eltype(A_mat))
    A_s = zero(T)
    @inbounds for j in eachindex(x)
        A_s += A_mat[i, j]*x[j]
    end
    α = -log(Z - B) + (B_i[i]/B)*(Z - 1)
    β = log((Z + m2*B)/(Z + m1*B))*A/(Δm*B)
    γ = (2/A)*A_s - B_i[i]/B
    return α + β*γ
end

"""
    component_fugacity(eos, cond, i, Z, forces, scalars)

Get fugacity of component `i` in a phase with compressibility `Z` and EOS constants `scalars`.
"""
Base.@propagate_inbounds function component_fugacity(eos::GenericCubicEOS, cond, i, Z, forces, scalars)
    lnψ_i = component_fugacity_coefficient(eos, cond, i, Z, forces, scalars)
    tmp = cond.z[i]*cond.p*exp(lnψ_i)
    return tmp
end

"""
Allocating version of `mixture_fugacities!`
"""
function mixture_fugacities(eos, arg...)
    n = number_of_components(eos)
    f = zeros(n)
    return mixture_fugacities!(f, eos, arg...)
end

"""
Compute fugacities for all components.
"""
function mixture_fugacities!(f, eos, cond, forces = force_coefficients(eos, cond), scalars = force_scalars(eos, cond, forces))
    Z = mixture_compressibility_factor(eos, cond, forces, scalars)
    @inbounds for i in eachindex(f)
        f[i] = component_fugacity(eos, cond, i, Z, forces, scalars)
    end
    return f
end

"""
Specialization of solve_roots for cubics
"""
solve_roots(::AbstractCubicEOS, P) = solve_cubic_positive_roots(P)

solve_cubic_positive_roots(P) = solve_cubic_positive_roots(P...)
function solve_cubic_positive_roots(a, b, c)
    Q = (a^2 - 3*b)/9
    R = (2*a^3 - 9*a*b + 27*c)/54
    M = R^2 - Q^3
    single_root = M > 0
    if single_root
        # Single real roots
        S = -sign(R)*(abs(R) + sqrt(M))^(1/3)
        if S == 0
            T = 0
        else
            T = Q/S
        end
        return S + T - a/3
    else
        # Three real roots
        theta = acos(R/sqrt(Q^3))
        r1 = -(2*sqrt(Q)*cos(theta/3)) - a/3
        r2 = -(2*sqrt(Q)*cos((theta + 2*pi)/3)) - a/3
        r3 = -(2*sqrt(Q)*cos((theta - 2*pi)/3)) - a/3
        return (r1, r2, r3)
    end
end

"""
    molar_masses(eos)

Returns a list of molecular weights (in kg/mol) for each component present in the input equation of state)
"""
molar_masses(eos::AbstractEOS) = molar_masses(eos.mixture)
molar_masses(mixture::MultiComponentMixture) = map((x) -> x.mw, mixture.properties)

"""
    component_names(eos)

Returns the list of component names for the input equation of state)
"""
component_names(eos::AbstractEOS) = component_names(eos.mixture)
component_names(mixture::MultiComponentMixture) = mixture.component_names

"""
    eostype(eos::AbstractEOS)
Returns the underlying EOS type.
"""
eostype(eos::AbstractEOS) = eos.type

function Base.summary(eos::AbstractEOS)
    n = number_of_components(eos)
    cnames = join(component_names(eos), ", ")
    return "$(eostype(eos)) EOS with $n components: $cnames"
end

