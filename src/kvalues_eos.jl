function evaluate_K(eos::KValuesEOS, cond)
    return eos.K_values_evaluator(cond)
end

function evaluate_K(eos::KValuesEOS{T, <:Any, <:Any}, cond) where T<:Tuple
    return eos.K_values_evaluator
end

function evaluate_K(eos::KValuesEOS{T, <:Any, <:Any}, cond) where T<:AbstractVector
    return eos.K_values_evaluator
end

function initial_guess_K(eos::KValuesEOS, cond)
    return evaluate_K(eos, cond)
end

function flash_storage(eos::KValuesEOS, cond = missing; kwarg...)
    return nothing
end

flash_storage(eos::KValuesEOS, cond, method, config::FlashConfig{PrintOutput, true}) where PrintOutput = nothing
flash_storage(eos::KValuesEOS, cond, method, config::FlashConfig{PrintOutput, false}) where PrintOutput = nothing

function flash_2ph!(storage, K, eos::KValuesEOS, cond, V = NaN; kwarg...)
    return solve_rachford_rice(K, cond.z, V)
end


function flash_2ph!(storage, K, eos::KValuesEOS, cond, V,
        config::FlashConfig{PrintOutput, UseDictStorage}; kwarg...) where {PrintOutput, UseDictStorage}
    return solve_rachford_rice(K, cond.z, V)
end
