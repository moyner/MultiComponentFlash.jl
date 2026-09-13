# ZudkevitchJoffe
function weight_ai(eos::GenericCubicEOS{E}, cond, i) where E<:ZudkevitchJoffe
    zj = eos.type
    mix = eos.mixture
    T = cond.T
    T_r = reduced_temperature(mix, cond, i)
    return eos.ω_a*zj.F_a(T, i)*T_r^(-0.5)
end

function weight_bi(eos::GenericCubicEOS{E}, cond, i) where E<:ZudkevitchJoffe
    zj = eos.type
    T = cond.T
    return eos.ω_b*zj.F_b(T, i)
end
