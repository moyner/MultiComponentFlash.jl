const CENTI_POISE_TO_PASCAL_SECOND = 1000.0
const MOL_TO_KMOL = 1000.0
const KELVIN_TO_RANKINE = 5/9
const PASCAL_TO_PSI = 1.450377377302092e-04

"""
    single_phase_label(mixture, cond)

Li's method for single-phase labeling of a mixture. Estimate of pure vapor/liquid.

Returns a vapor fraction that is either 1.0 (=pure vapor) or 0.0 (=pure liquid).
"""
function single_phase_label(mixture::MultiComponentMixture, cond)
    # Li's method for phase labeling.
    z = cond.z
    T_c = 0.0
    V_c = 0.0
    @inbounds for (i, m) in enumerate(mixture.properties)
        V = m.V_c*z[i]
        T_c += m.T_c*V
        V_c += V
    end
    T_c /= V_c
    return Float64(cond.T > T_c)
end

single_phase_label(eos::AbstractEOS,cond) = single_phase_label(eos.mixture,cond)

"""
    lbc_viscosity(eos, p, T, ph; <keyword arguments>)

Compute the viscosity of a mixture using the [Lohrenz-Bray-Clark](https://doi.org/10.2118/915-PA) correlation.
"""
function lbc_viscosity(eos, p, temperature, ph::FlashedPhase{T}; coeff = (0.1023, 0.023364, 0.058533, -0.040758, 0.0093324), shift = -1e-4) where T
    z = ph.mole_fractions
    properties = eos.mixture.properties
    mw_mix, P_pc, T_pc, V_pc = pseudo_critical_properties(properties, z)
    mu_atm = atmospheric_mu_estimate(properties, z, temperature)
    e_mix = mixture_viscosity_parameter(mw_mix, T_pc, P_pc)
    # From Jossi et al
    # coeffs = [0.1023, 0.023364, 0.058533, -0.040758, 0.0093724, 1e-4]
    # From LBC paper
    # coeffs = [0.1023, 0.023364, 0.058533, -0.040758, 0.0093324, -1e-4]
    # Compute reduced density
    V = molar_volume(eos, p, temperature, ph)
    rho_r = V_pc/V
    corr = zero(T)
    for (i, c) in enumerate(coeff)
        corr += c*rho_r^(i-1)
    end
    mu_correction = (corr^4 + shift)/e_mix
    # Final expression is compound and given in centi poise. We convert to Pa s
    # instead for strict SI outputs.
    mu = 1e-3*(mu_atm + mu_correction)
    return mu::T
end

function atmospheric_mu_estimate(props, z::V, temperature) where V<:AbstractVector{T} where T
    a = zero(T)
    b = zero(T)
    @inbounds for (prop, zi) in zip(props, z)
        mw = prop.mw
        T_c = prop.T_c
        p_c = prop.p_c
        # Add contributions to atmospheric mu
        T_r = temperature/T_c
        e_i = mixture_viscosity_parameter(mw, T_c, p_c)
        if T_r > 1.5
            mu_i = 17.78e-5*(4.58*T_r - 1.67)^0.625
        else
            mu_i = 34e-5*T_r^(0.94)
        end
        tmp = sqrt(1000.0*mw)*zi
        a += tmp*mu_i/e_i
        b += tmp
    end
    return a/b
end

@inline function mixture_viscosity_parameter(molar_mass, T_c, p_c)
    # Note: Original paper uses g/mol and pressure in atmospheres while we use
    # kg/mol and Pa. Convert accordingly
    A = T_c^(1/6)
    B = sqrt(1000.0*molar_mass)*(p_c/101325.0)^(2/3)
    return A/B
end

function pseudo_critical_properties(props, z::V) where V<:AbstractVector{T} where T
    P_pc = zero(T)
    T_pc = zero(T)
    Vc = zero(T)
    mwc = zero(T)
    @inbounds for (prop, zi) in zip(props, z)
        # Accumulate weighted properties
        mwc += zi*prop.mw
        P_pc += zi*prop.p_c
        T_pc += zi*prop.T_c
        Vc += zi*prop.V_c
    end
    return (mwc, P_pc, T_pc, Vc)
end

@inline function molar_volume(R, p, T, Z)
    return R*T*Z/p
end

function michelsen_critical_point_measure_storage(eos; T = Float64, static_size = true)
    n = number_of_components(eos)
    T∂ = ForwardDiff.Dual{nothing, T, n}
    mole_numbers_ad = Vector{T∂}(undef, n)
    for i in 1:n
        partials = ForwardDiff.Partials{n, T}(Tuple(map(j -> Float64(i == j), 1:n)))
        mole_numbers_ad[i] = T∂(1.0/n, partials)
    end
    z = mole_numbers_ad./sum(mole_numbers_ad)
    if static_size
        z = MVector{n, T∂}(z)
    end
    B = zeros(T, n, n)
    c = (p = 101325.0, T = 303.15, z = z)
    forces = force_coefficients(eos, c; static_size = static_size)
    return (B = B, forces = forces, z = z)
end

function michelsen_critical_point_measure(eos, p, T, mole_numbers; kwarg...)
    S = michelsen_critical_point_measure_storage(eos; kwarg...)
    michelsen_critical_point_measure!(S, eos, p, T, mole_numbers)
end

function michelsen_critical_point_measure!(S, eos, p, T, mole_numbers)
    # From Michelsen (1982): The Isothermal Flash Problem. Part I: Stability
    # Estimate the distance to the critical point through smallest eigenvalue.
    n = length(mole_numbers)
    z = S.z
    B = S.B
    forces = S.forces
    for i in eachindex(mole_numbers, z)
        z[i] += mole_numbers[i] - z[i].value
    end
    zt = sum(z)
    @. z /= zt
    c = (p = p, T = T, z = z)
    forces = force_coefficients!(forces, eos, c)
    scalars = force_scalars(eos, c, forces)
    Z = mixture_compressibility_factor(eos, c, forces, scalars)
    for i in eachindex(mole_numbers)
        n_i = mole_numbers[i]
        f_i = MultiComponentFlash.component_fugacity_coefficient(eos, c, i, Z, forces, scalars)
        for j in eachindex(mole_numbers)
            n_j = mole_numbers[j]
            B[i, j] = (i == j) + sqrt(n_i*n_j)*(f_i.partials[j])
        end
    end
    return eigmin(B)
end
