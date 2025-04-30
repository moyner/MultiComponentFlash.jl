function secant_update(T_curr, T_prev, R_curr, R_prev)        
    return T_curr - R_curr*(T_curr - T_prev)/(R_curr - R_prev)
end

function secant_intialization(T_0, R_0, F_L, F_V, Cp_L, Cp_V)        
    return T_0 - (R_0/(F_L*Cp_L + F_V*Cp_V))
end

function regular_falsi_update(T_left, T_right, R_curr, R_left, R_right)
    return T_left  - R_curr*(T_left-T_right)/(R_left - R_right)
end


function enthalpy_ideal(T,T_r, Cp)
# Input: T   -> Current Temperature, 
#        T_r -> Reference Temperature 
#        Cp  -> Matrix nc x 4. The coeficients 
#                   nc = number of components
#                   each component needs   with 4 Heat Capacity Coefficients

nc, _ = size(Cp)
h_c = 0.0;

# Output: Ideal Mixture Enthalpy
    p_exp = collect(1:4)
    for i = 0:nc
        # x_c -> fraction of c in a given phase
        h_c = h_c + (x_c) * enthalpy_ideal_component(T,T_r, Cp[i,:]);
    end
    return hc
end

function enthalpy_ideal_component(T,T_r, cp)
# Input: T   -> Current Temperature, 
#        T_r -> Reference Temperature 
#        Cp  -> Matrix 1 x 4. The coeficients 
#                   nc = number of components
#                   each component needs   with 4 Heat Capacity Coefficients
# Output: Ideal Component Enthalpy
    p_exp = collect(1:4)
    return sum((fill(T,4,).^p_exp  - fill(T_r,4,).^p_exp ).*cp')
end


function enthalpy_departure(eos, cond, i, Z, forces, scalars) 
    # Input: T   -> Current Temperature, 
    #        T_r -> Reference Temperature 
    # Output: Departure Enthalpy
    # for each component
    #  Hc -> Array n_c components x 1
    R =   8.314 # J / mol·K
    @inbounds for i in eachindex(c)
        #  d_component_fugacity_dt = d/dT(component_fugacity) with AD
        Hc[i] = d_component_fugacity_dt(eos, cond, i, Z, forces, scalars)    
    end
    return -sum(H_c)*R*(T^2) 
end

function enthalpy_phase(eos, cond, i, Z, forces, scalars, T,T_r, Cp)
    # Input: eos -> Eqution of State
    #        T   -> Current Temperature, 
    #        T_r -> Reference Temperature 
    #        Cp  -> Matrix with  4 Heat Capacity Coefficients
    # Output: Phase Enthalpy
    return enthalpy_ideal(T,T_r, Cp) + enthalpy_departure(eos, cond, i, Z, forces, scalars) 
end

