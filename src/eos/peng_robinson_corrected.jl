
function weight_ai(eos::GenericCubicEOS{PengRobinsonCorrected}, cond, i)
    mix = eos.mixture
    m = molecular_property(mix, i)
    a = acentric_factor(m)
    T_r = reduced_temperature(mix, cond, i)
    aa = a*a
    if a > 0.49
        # Use alternate expression.
        D = (0.379642 + 1.48503*a - 0.164423*aa + 0.016666*a*aa)
    else
        # Use standard expression.
        D = (0.37464 + 1.54226*a - 0.26992*aa)
    end
    tmp = 1.0 + D*(1.0-T_r^0.5)
    return eos.ω_a*(tmp*tmp);
end
