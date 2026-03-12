
# PengRobinson specialization
function static_coefficients(::AbstractPengRobinson)
    return (0.457235529, 0.077796074, 1.0 + sqrt(2), 1.0 - sqrt(2))
end

function weight_ai(eos::GenericCubicEOS{T}, cond, i) where T<:AbstractPengRobinson
    mix = eos.mixture
    m = molecular_property(mix, i)
    a = acentric_factor(m)
    T_r = reduced_temperature(mix, cond, i)
    α_half = (1.0 + (0.37464 + 1.54226*a - 0.26992*a*a)*(1-T_r^0.5))
    return eos.ω_a*(α_half*α_half)
end
