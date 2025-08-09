# SoaveRedlichKwong
function weight_ai(eos::GenericCubicEOS{SoaveRedlichKwong}, cond, i)
    mix = eos.mixture
    m = molecular_property(mix, i)
    a = acentric_factor(m)
    T_r = reduced_temperature(mix, cond, i)
    return eos.ω_a*(1 + (0.48 + 1.574*a - 0.176*a^2)*(1-T_r^(1/2)))^2;
end
