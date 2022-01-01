"""
    single_phase_label(mixture, cond)

Li's method for single-phase labeling of a mixture. Estimate of pure vapor/liquid.

Returns a vapor fraction that is either 1.0 (=pure vapor) or 0.0 (=pure liquid).
"""
function single_phase_label(mixture, cond)
    # Li's method for phase labeling.
    z = cond.z
    T_c = 0.0
    V_c = 0.0
    for (i, m) in enumerate(mixture.properties)
        V = m.V_c*z[i]
        T_c += m.T_c*V
        V_c += V
    end
    return Float64(cond.T < T_c)
end

"""
    lbc_viscosity(eos, p, T, ph; <keyword arguments>)

Compute the viscosity of a mixture using the [Lohrenz-Bray-Clark](https://doi.org/10.2118/915-PA) correlation.
"""
function lbc_viscosity(eos, p, temperature, ph::FlashedPhase{T}; coeff = (0.1023, 0.023364, 0.058533, -0.040758, 0.0093324), shift = -1e-4) where T
    molfactor = 1000
    Rankine = 5/9
    psia = 6.894757293168360e+03
    
    rho = 1/molar_volume(eos, p, temperature, ph)
    z = ph.mole_fractions
    
    properties = eos.mixture.properties
    P_pc = zero(T)
    T_pc = zero(T)
    Vc = zero(T)
    mwc = zero(T)

    a = zero(T)
    b = zero(T)
    @inbounds for (prop, zi) in zip(properties, z)
        mw = molfactor*prop.mw
        T_c = prop.T_c/Rankine
        p_c = prop.p_c/psia
        V_c = prop.V_c
        # Accumulate weighted properties
        mwc += zi*mw
        P_pc += zi*p_c
        T_pc += zi*T_c
        Vc += zi*V_c
        # Add contributions to atmospheric mu
        tr = temperature/T_c
        mwi = sqrt(mw)
        e_i = (5.4402*T_c^(1/6))/(mwi*p_c^(2/3)*(1e-3))
        large = tr > 1.5
        if large
            mu_i = (17.78e-5*(4.58*tr - 1.67)^0.625)/e_i
        else
            mu_i = 34e-5*tr^(0.94)/e_i
        end
        a += zi*mu_i*mwi
        b += zi*mwi
    end
    mu_atm = a/b
    e_mix = 5.4402*(T_pc)^(1/6)/((mwc)^(1/2)*(P_pc)^(2/3)*(1e-3))
    rhor = Vc*rho
    # From Jossi et al
    # coeffs = [0.1023, 0.023364, 0.058533, -0.040758, 0.0093724, 1e-4]
    # From LBC paper
    # coeffs = [0.1023, 0.023364, 0.058533, -0.040758, 0.0093324, -1e-4]
    corr = zero(T)
    for (i, c) in enumerate(coeff)
        corr += c*rhor^(i-1)
    end
    mu = mu_atm + (corr^4 + shift)/e_mix
    return mu::T
end

@inline function molar_volume(R, p, T, Z)
    return R*T*Z/p
end
