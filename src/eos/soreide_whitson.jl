# Søreide-Whitson Peng-Robinson variant
function weight_ai(eos::GenericCubicEOS{T}, cond, i) where T<:SoreideWhitson
    sw = eos.type
    mix = eos.mixture
    m = molecular_property(mix, i)
    a = acentric_factor(m)
    T_r = reduced_temperature(mix, cond, i)
    if sw.component_types[i] == COMPONENT_H2O
        # Use the water-specific expression.
        w1, w2, w3 = sw.water_coefficients
        α_half = 1.0 + w1*(1.0 - w2*T_r) + w3*(T_r^(-3)-1.0)
    else
        α_half = 1.0 + (0.37464 + 1.54226*a - 0.26992*a*a)*(1-T_r^0.5)
    end
    return eos.ω_a*(α_half*α_half)
end

function binary_interaction(eos::GenericCubicEOS{T}, i::Int, j::Int, cond) where T<:SoreideWhitson
    phase = get_phase(cond)
    sw = eos.type
    if phase == :liquid
        cat_i = eos.type.component_types[i]
        cat_j = eos.type.component_types[j]

        # Use eq 12 from paper
        i_is_h2o = cat_i == COMPONENT_H2O
        j_is_h2o = cat_j == COMPONENT_H2O
        # in paper: i is hc, j is h2o
        if i_is_h2o || j_is_h2o
            # Special cases: H2O-H2S, H2O-N2 and H2O-CO2
            if i_is_h2o
                cat_other = cat_j
                ix = j
            else
                cat_other = cat_i
                ix = i
            end
            hc_property = molecular_property(eos.mixture, ix)
            acf = acentric_factor(hc_property)
            T_r = cond.T/critical_temperature(hc_property)
            if cat_other == COMPONENT_CO2
                bic = soreide_whitson_co2_aqueous_bic(sw::SoreideWhitson, acf, T_r)
            elseif cat_other == COMPONENT_H2S
                bic = soreide_whitson_h2s_aqueous_bic(sw::SoreideWhitson, acf, T_r)
            elseif cat_other == COMPONENT_N2
                bic = soreide_whitson_n2_aqueous_bic(sw::SoreideWhitson, acf, T_r)
            else
                bic = soreide_whitson_hc_aqueous_bic(sw, acf, T_r)
            end
        else
            # Use default.
            bic = binary_interaction(eos.mixture, i, j)
        end
    elseif phase == :vapor
        bic = binary_interaction(eos.mixture, i, j)
    else
        error("Søreide-Whitson requires the phase state to be specified for binary interactions (was: $phase).")
    end
    return bic
end

function soreide_whitson_n2_aqueous_bic(sw::SoreideWhitson, acf_i, T_ri)
    c_sw = sw.molality
    return -1.70235*(1 + 0.025587*c_sw^0.75) + 0.44338*(1 + 0.08126*c_sw^0.75)*T_ri
end

function soreide_whitson_h2s_aqueous_bic(sw::SoreideWhitson, acf_i, T_ri)
    return -0.20441 + 0.234267*T_ri
end

function soreide_whitson_co2_aqueous_bic(sw::SoreideWhitson, acf_i, T_ri)
    c_sw = sw.molality
    return -0.31092*(1 + 0.15587*c_sw^0.7505) +0.23580 * (1 + 0.17837*c_sw^0.979)*T_ri - 21.2566*exp(-6.7222*T_ri- c_sw)
end

function soreide_whitson_hc_aqueous_bic(sw::SoreideWhitson, acf_i, T_ri)
    c_sw = sw.molality
    a_0, a_1, a_2 = sw.A
    a_mw_0, a_mw_1, a_mw_2 = sw.A_mw
    α_0, α_1, α_2 = sw.alphas

    A_0 = a_0 + a_mw_0*acf_i^(-0.1)
    A_1 = a_1 + a_mw_1*acf_i
    A_2 = a_2 + a_mw_2*acf_i
    return A_0*(1 + α_0*c_sw) + A_1*T_ri*(1 + α_1*c_sw) + A_2*T_ri*(1 + α_2*c_sw)
end
