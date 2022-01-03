abstract type AbstractEOS end
abstract type AbstractCubicEOS <: AbstractEOS end

@inline reduced_pT(mix, c, i) = (reduced_pressure(mix, c, i), reduced_temperature(mix, c, i))
@inline reduced_pT(eos::AbstractEOS, c, i) = (reduced_pressure(eos.mixture, c, i), reduced_temperature(eos.mixture, c, i))

"""
GenericCubicEOS is an implementation of generalized cubic equations of state.

Many popular cubic equations can be written in a single form with a few changes in
definitions for the terms (they are, after all, all cubic in form). References:

 1. [Cubic Equations of State-Which? by J.J. Martin](https://doi.org/10.1021/i160070a001)
 2. [Simulation of Gas Condensate Reservoir Performance  by K.H. Coats](https://doi.org/10.2118/10512-PA)

"""
struct GenericCubicEOS{T, R, N} <: AbstractCubicEOS
    type::T
    mixture::MultiComponentMixture{R, N}
    m_1::R
    m_2::R
    ω_a::R
    ω_b::R
end

"""
Abstract base type for generalized cubics.

Any subtype may override the any of following internal functions to alter the EOS:

    * static_coefficients
    * weight_ai
    * weight_bi

In addition, the GenericCubicEOS is parametric on the specific cubic type for further
dispatch modifications.
"""
abstract type AbstractGeneralizedCubic end

"""
Specializes the GenericCubicEOS to the Peng-Robinson cubic equation of state.
"""
struct PengRobinson <: AbstractGeneralizedCubic end
"""
Specializes the GenericCubicEOS to the Zudkevitch-Joffe cubic equation of state.

The Zudkevitch-Joffe equations of state allows for per-component functions of 
temperature that modify the `weight_ai` and `weight_bi` functions. These additional
fitting parameters allows for more flexibility when matching complex mixtures.
"""
struct ZudkevitchJoffe <: AbstractGeneralizedCubic
    F_a
    F_b
    function ZudkevitchJoffe(; F_a = (T, c_i) -> 1.0, F_b = (T, c_i) -> 1.0)
        new(F_a, F_b)
    end
end

"""
Specializes the GenericCubicEOS to the Soave-Redlich-Kwong cubic equation of state.
"""
struct SoaveRedlichKwong <: AbstractGeneralizedCubic end
"""
Specializes the GenericCubicEOS to the Redlich-Kwong cubic equation of state.
"""
struct RedlichKwong <: AbstractGeneralizedCubic end

"""
    GenericCubicEOS(mixture, [type = PengRobinson()])

Instantiate a generic cubic equation-of-state for a `MultiComponentMixture` and 
a specified EOS.

Currently supported choices for type:

    1. `PengRobinson` (default)
    2. `ZudkevitchJoffe`
    3. `RedlichKwong`
    4. `SoaveRedlichKwong`
"""
function GenericCubicEOS(mixture, type = PengRobinson())
    ω_a, ω_b, m_1, m_2 = static_coefficients(type)
    setup = (ω_a = ω_a, ω_b = ω_b, m_1 = m_1, m_2 = m_2, type = type)
    return GenericCubicEOS(setup, mixture)
end

function GenericCubicEOS(setup::NamedTuple, mixture)
    return GenericCubicEOS(setup.type, mixture, setup.m_1, setup.m_2, setup.ω_a, setup.ω_b)
end

"""
    number_of_components(eos)

Return number of components for the underlying mixture of the EOS.
"""
number_of_components(e::AbstractEOS) = number_of_components(e.mixture)


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

# PengRobinson specialization
function static_coefficients(::PengRobinson)
    return (0.4572355, 0.0779691, 1 + sqrt(2), 1 - sqrt(2))
end

function weight_ai(eos::GenericCubicEOS{PengRobinson}, cond, i)
    mix = eos.mixture
    m = molecular_property(mix, i)
    a = acentric_factor(m)
    T_r = reduced_temperature(mix, cond, i)
    return eos.ω_a*(1 + (0.37464 + 1.54226*a - 0.26992*a^2)*(1-T_r^(1/2)))^2;
end

# ZudkevitchJoffe
function weight_ai(eos::GenericCubicEOS{ZudkevitchJoffe}, cond, i)
    zj = eos.type
    mix = eos.mixture
    T = cond.T
    T_r = reduced_temperature(mix, cond, i)
    return eos.ω_a*zj.F_a(T, i)*T_r^(-1/2)
end

function weight_bi(eos::GenericCubicEOS{ZudkevitchJoffe}, cond, i)
    zj = eos.type
    T = cond.T
    return eos.ω_b*zj.F_b(T, i)
end

# SoaveRedlichKwong
function weight_ai(eos::GenericCubicEOS{SoaveRedlichKwong}, cond, i)
    mix = eos.mixture
    m = molecular_property(mix, i)
    a = acentric_factor(m)
    T_r = reduced_temperature(mix, cond, i)
    return eos.ω_a*(1 + (0.48 + 1.574*a - 0.176*a^2)*(1-T_r^(1/2)))^2;
end

# Generic part
binary_interaction(eos::AbstractEOS, i, j) = binary_interaction(eos.mixture, i, j)
binary_interaction(mixture::MultiComponentMixture{R}, i, j) where {R} = binary_interaction(mixture.binary_interaction, i, j)::R
binary_interaction(::Nothing, i, j) = 0.0
Base.@propagate_inbounds binary_interaction(B::AbstractMatrix, i, j) = B[i, j]

"""
    mixture_compressibility_factor(eos, cond, [forces, scalars])

Compute compressibility factor for current conditions. Forces and scalars can be 
provided if they are already known.

The compressibility factor adjusts the ideal gas law to account for non-linear behavior: ``pV = nRTZ``
"""

function mixture_compressibility_factor(eos::AbstractCubicEOS, cond,
                                            forces = force_coefficients(eos, cond),
                                            scalars = force_scalars(eos, cond, forces))
    poly = eos_polynomial(eos, forces, scalars)
    roots = solve_roots(eos, poly)
    r = pick_root(eos, roots, cond, forces, scalars)
    return r
end

minimum_allowable_root(eos::AbstractCubicEOS, forces, scalars) = scalars.B
minimum_allowable_root(eos, forces, scalars) = 1e-16

@inline pick_root(eos, roots::Real, cond, forces, scalars) = roots
function pick_root(eos, roots, cond, forces, scalars)
    r_ϵ = minimum_allowable_root(eos, forces, scalars)
    max_r = maximum(roots)
    min_r = minimum((x) -> x > r_ϵ ? x : Inf, roots)
    if min_r == max_r
        r = min_r
    else
        function Gibbs(Z)
            E = 0.0
            z = cond.z
            @inbounds for i in eachindex(z)
                ϕ = component_fugacity_coefficient(eos, cond, i, Z, forces, scalars)
                E += z[i]*ϕ
            end
            return E
        end
        if Gibbs(min_r) < Gibbs(max_r)
            r = min_r
        else
            r = max_r
        end
    end
    return r
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
        A_ij = @MMatrix zeros(eT, n, n)
        A_i = @MVector zeros(eT, n)
        B_i = @MVector zeros(eT, n)
    else
        A_ij = zeros(eT, n, n)
        A_i = zeros(eT, n)
        B_i = zeros(eT, n)
    end
    coeff = (A_ij = A_ij, A_i = A_i, B_i = B_i)
    force_coefficients!(coeff, eos, cond)
end

"""
    force_coefficients!(coeff, eos, cond)

In-place update of force coefficients.

See also [`force_coefficients`](@ref)
"""
function force_coefficients!(coeff, eos::AbstractCubicEOS, arg...)
    update_attractive_linear!(coeff.A_i, eos, arg...)
    update_attractive_quadratic!(coeff.A_ij, coeff.A_i, eos, arg...)
    update_repulsive!(coeff.B_i, eos, arg...)
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
            a = sqrt(A_i[i]*A_i[j])*(one(T) - binary_interaction(eos, i, j))
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
    return cubic_scalars(forces.A_ij, forces.B_i, cond.z)
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
function component_fugacity_coefficient(eos::AbstractCubicEOS, cond, i, Z, forces, scalars)
    # NOTE: This returns ln(ψ), not ψ!
    m1 = eos.m_1
    m2 = eos.m_2
    return component_fugacity_coefficient_cubic(m1, m2, cond.z, Z, scalars.A, scalars.B, forces.A_ij, forces.B_i, i)
end

function component_fugacity_coefficient_cubic(m1, m2, x, Z, A, B, A_mat, B_i, i)
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
