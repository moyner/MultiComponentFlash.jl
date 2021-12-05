
"Abstract type for all flash types"
abstract type AbstractFlash end
"Abstract type for all flash types that use Newton in some form"
abstract type AbstractNewtonFlash <: AbstractFlash end

"""
Flash method that uses successive subtition.

Unconditionally convergent, does not require derivatives, but is very slow around critical regions.

See also: [`flash_2ph!`](@ref), [`NewtonFlash`](@ref) [`SSINewtonFlash`](@ref)
"""
struct SSIFlash <: AbstractFlash end

"""
    NewtonFlash([dMax = 0.2])

Flash using Newton's method for zero solve.

Only conditionally convergent, but has better convergence rate than SSI.

# Arguments
- `dMax`: dampening factor for the newton iteration

See also: [`flash_2ph!`](@ref), [`SSIFlash`](@ref) [`SSINewtonFlash`](@ref)
"""
@Base.kwdef struct NewtonFlash <: AbstractNewtonFlash
    dMax::Float64 = 0.2
end

"""
    SSINewtonFlash([swap_iter = 5,dMax = 0.2])

Perform a number of SSI iterations, followed by Newton until convergence.

See also: [`flash_2ph!`](@ref), [`SSIFlash`](@ref) [`NewtonFlash`](@ref)
"""
@Base.kwdef struct SSINewtonFlash <: AbstractNewtonFlash
    swap_iter::Int = 5
    dMax::Float64 = 0.2
end

"""
    flash_2ph(eos, c, [K], [V]; <keyword arguments>)

Perform two-phase flash with a given EOS under a set of specific conditions. Returns vapor fraction. Modifies K in-place.

Given a mixture with pressure, temperature and mole fractions, this routine performs a vapor-liquid
equilibrium calculation after a stability test.

Two outcomes are possible:
1. A single-phase condition prevails (returned vapor fraction is `NaN`) and a single phase (liquid or vapor) is stable.
2. A two-phase condition is possible. The routine produces K-values and vapor fraction so that the following holds:
    1. Isofugacity constraint for all components (``f_{li} = f_{vi}``)
    2. Molar balance for all components (``(1-V) x_i - V y_i - z_i``)
    3. Unity condition (``\\sum_i (x_i - y_i) = 0 ``)

# Arguments
- `eos`: the equation-of-state to be used for the flash
- `c`: conditions to flash the mixture at on the form `(p = 10e5, T = 303.15, z = [0.5, 0.3, 0.2])`
- `K`: optionally a buffer of length `number_of_components(eos)` used to hold K-values. Modified in-place.
- `V`: optionally the initial guess for V. If this value is not `NaN`, the stability check will be skipped.

# Keyword arguments
- `method = SSIFlash()`: Flash method to use. Can be `SSIFlash()`, `NewtonFlash()` or `SSINewtonFlash()`.
- `tolerance = 1e-8`: Tolerance for the convergence criterion. Relative to ``\\|1-R_i\\|_\\infty`` where ``R_i = f_{il}/f_{iv}``
- `maxiter = 10000`: Maximum nubmer of iterations for both stability tests and the main flash.
- `verbose = false`: Emit extra information during solve.
- `extra_out = false`: Return `(V, K, conv)` where `conv` contains iterations and oncergence status instead of just `V`.

See also: [`flash_2ph!`](@ref), [`single_phase_label`](@ref)
"""
function flash_2ph(eos, c::T, K = initial_guess_K(eos, c); method = SSIFlash(), kwarg...) where T
    nc = number_of_components(eos)
    @assert hasfield(T, :p)
    @assert hasfield(T, :T)
    @assert hasfield(T, :z)
    @assert length(c.z) == nc
    @assert length(K) == nc
    method::AbstractFlash

    S = flash_storage(eos, c, method = method)
    flash_2ph!(S, K, eos, c, update_forces = false; method = method, kwarg...)
end

"""
    flash_2ph!(storage, K, eos, c, [V]; <keyword arguments>)

Non-allocating version of [`flash_2ph`](@ref) where storage is pre-allocated.

Useful if you are performing many flashes of the same system with varying conditions.

# Arguments
- `storage`: Should be output from `flash_storage(eos, c, method = method)`. Preallocated storage.
Remaining arguments documented in [`flash_2ph`](@ref).

# Keyword arguments
- `update_forces = true`: Update the `p, T` dependent forces in storage initially.

"""
function flash_2ph!(storage, K, eos, c, V = NaN;
                                method = SSIFlash(), verbose = false, maxiter = 10000,
                                tolerance = 1e-8, extra_out = false, update_forces = true, kwarg...)
    z_min = MINIMUM_COMPOSITION
    z = c.z
    for i in eachindex(z)
        @inbounds z[i] = max(z[i], z_min)
    end
    @assert all(isfinite, K)
    forces = storage.forces
    if update_forces
        force_coefficients!(forces, eos, c)
    end
    single_phase_init = isnan(V) || V == 1.0 || V == 0.0
    if single_phase_init
        stable = stability_2ph!(storage, K, eos, c, maxiter = maxiter, verbose = verbose; kwarg...)
    else
        stable = false
    end
    converged = false
    if stable
        i = 0
    else
        i = 1
        V = solve_rachford_rice(K, z, V)
        while !converged
            V, ϵ = flash_update!(K, storage, method, eos, c, forces, V, i)
            converged = ϵ ≤ tolerance
            if converged
                if verbose
                    @info "Flash done in $i iterations." V K
                end
                break
            end
            if i == maxiter
                @warn "Flash did not converge in $i iterations. Final ϵ = $ϵ > $tolerance = tolerance"
            end
            i += 1
        end
    end
    if extra_out
        return (V, K, (its = i, converged = converged))
    else
        return V
    end
end

"""
    flash_storage(eos, [c]; <keyword arguments>)

Pre-allocate storage for `flash_2ph!`.

# Arguments
- `eos`: the equation-of-state to be used for the flash
- `c`: conditions for the mixture. Types in conditions must match later usage (e.g. for use with ForwardDiff).

# Keyword arguments
- `method = SSIFlash()`: Flash method to use. Can be `SSIFlash()`, `NewtonFlash()` or `SSINewtonFlash()`.
- `static_size = false`: Use `SArrays` and `MArrays` for fast flash, but slower compile times.
- `inc_jac`: Allocate storage for Newton/Jacobian. Required for Newton (and defaults to `true` for that method) or for `diff_externals`.
- `diff_externals = false`: Allocate storage for matrix inversion required to produce partial derivatives of flash using `set_partials`.

See also: [`flash_2ph!`](@ref) [`set_partials`](@ref)
"""
function flash_storage(eos, cond = (p = 10e5, T = 273.15, z = zeros(number_of_components(eos)));  method = SSIFlash(), kwarg...)
    out = Dict{Symbol,Any}()
    d = flash_storage_internal!(out, eos, cond, method; kwarg...)
    # Convert to named tuple
    return NamedTuple(pairs(d))
end

function flash_storage_internal!(out, eos, cond, method; inc_jac = isa(method, AbstractNewtonFlash), static_size = false, kwarg...)
    n = number_of_components(eos)
    out[:forces] = force_coefficients(eos, cond, static_size = static_size)
    if static_size
        alloc_vec = () -> @MVector zeros(n)
    else
        alloc_vec = () -> zeros(n)
    end
    out[:x] = alloc_vec()
    out[:y] = alloc_vec()

    out[:buffer1] = alloc_vec()
    out[:buffer2] = alloc_vec()
    if inc_jac
        flash_storage_internal_newton!(out, eos, cond, method, static_size = static_size; kwarg...)
    end
    return out
end

function flash_storage_internal_newton!(out, eos, cond, method; static_size = false, diff_externals = false, kwarg...)
    n = number_of_components(eos)
    np = 2*n + 1
    primary_ad(ix) = get_ad(0.0, np, :Flash, ix)
    V_ad = primary_ad(np)
    T = typeof(V_ad)
    if static_size
        x_ad = @MVector zeros(T, n)
        y_ad = @MVector zeros(T, n)
        r = @MVector zeros(np)
        J = @MMatrix zeros(np, np)
    else
        x_ad = zeros(T, n)
        y_ad = zeros(T, n)
        r = zeros(np)
        J = zeros(np, np)
    end
    out[:r] = r
    out[:J] = J

    for i = 1:n
        x_ad[i] = primary_ad(i)
        y_ad[i] = primary_ad(i+n)
    end
    out[:AD] = (x = x_ad, y = y_ad, V = V_ad)
    if diff_externals
        flash_storage_internal_inverse!(out, eos, cond, method, static_size = static_size; kwarg...)
    end
    return out
end

function flash_storage_internal_inverse!(out, eos, cond, method; static_size = false, npartials = nothing)
    n = number_of_components(eos)
    np = length(out[:r])
    external_partials = n + 2 # p, T, z_1, ... z_n
    secondary_ad(ix) = get_ad(0.0, external_partials, :InverseFlash, ix)
    p_ad = secondary_ad(1)
    T_ad = secondary_ad(2)
    T_cond = typeof(p_ad)
    
    if static_size
        z_ad = @MVector zeros(T_cond, n)
        J_inv = @MMatrix zeros(np, external_partials)
    else
        z_ad = zeros(T_cond, n)
        J_inv = zeros(np, external_partials)
    end
    out[:J_inv] = J_inv
    for i = 1:n
        z_ad[i] = secondary_ad(i+2)
    end
    out[:AD_cond] = (p = p_ad, T = T_ad, z = z_ad)
    if !isnothing(npartials)
        out[:buf_inv] = zeros(npartials)
    end
end

function get_ad(v::T, npartials, tag, diag_pos = nothing) where {T<:Real}
    D = ntuple(x -> T.(x == diag_pos), npartials)
    partials = ForwardDiff.Partials{npartials, T}(D)
    v = ForwardDiff.Dual{tag, T, npartials}(v, partials)
    return v
end

update_value(v::Real, newv::Real) = newv
function update_value(v::T, newv::Real) where T<:ForwardDiff.Dual
    P = v.partials
    return T(newv, P)
end

function flash_update!(K, storage, type::SSIFlash, eos, cond, forces, V, iteration)
    z = cond.z
    x, y = storage.x, storage.y
    p, T = cond.p, cond.T
    return ssi!(K, p, T, x, y, z, V, eos, forces)
end

function ssi!(K, p::F, T::F, x, y, z, V::F, eos, forces) where {F<:Real}
    # Initialize conditions for vapor and liquid phases based on K-values
    x = liquid_mole_fraction!(x, z, K, V)
    y = vapor_mole_fraction!(y, x, K)

    liquid = (p = p, T = T, z = x)
    vapor = (p = p, T = T, z = y)

    Z_l, s_l = prep(eos, liquid, forces)
    Z_v, s_v = prep(eos, vapor, forces)

    ϵ = zero(F)
    @inbounds for c in eachindex(K)
        f_l = component_fugacity(eos, liquid, c, Z_l, forces, s_l)
        f_v = component_fugacity(eos, vapor, c, Z_v, forces, s_v)
        r = f_l/f_v
        K[c] *= r
        ϵ = max(ϵ, abs(1-r))
    end
    V = solve_rachford_rice(K, z, V)
    return (V, ϵ)::Tuple{F, F}
end

cap_z(z) = min(max(z, MINIMUM_COMPOSITION), one(z))
cap_unit(v) = min(max(v, zero(z)), one(z))
cap_VL(v) = min(max(v, MINIMUM_COMPOSITION), 1 - MINIMUM_COMPOSITION)

function flash_update!(K, storage, type::NewtonFlash, eos, cond, forces, V, iteration)
    x, y = storage.x, storage.y
    z = cond.z
    x = liquid_mole_fraction!(x, z, K, V)
    y = vapor_mole_fraction!(y, x, K)
    # Newton part
    Δ = update_and_solve!(storage, eos, cond, forces, x, y, V)
    newton_dampen!(type.dMax, Δ)
    V, ϵ = update_newton_from_increment!(K, x, y, V, Δ)
    return (V, ϵ)
end

function update_newton_from_increment!(K, x, y, V, Δ)
    ϵ = zero(eltype(K))
    n = length(x)
    @inbounds for i in 1:n
        x_n = cap_z(x[i] - Δ[i])
        y_n = cap_z(y[i] - Δ[i+n])
        K_next = y_n/x_n
        r = K_next/K[i]
        ϵ = max(ϵ, abs(1-r))
        # Assign new value, overwriting old
        K[i] = K_next
    end
    V = cap_VL(V - Δ[end])
    return (V, ϵ)
end

function get_newton_ad_values(AD, x, y, V)
    xAD = AD.x
    yAD = AD.y
    @inbounds for i in eachindex(x)
        xAD[i] = update_value(xAD[i], x[i])
        yAD[i] = update_value(yAD[i], y[i])
    end
    VAD = update_value(AD.V, V)
    return (xAD, yAD, VAD)
end

function update_primary_jacobian!(storage, eos, cond, forces, x, y, V)
    J = storage.J
    r = storage.r
    # Transfer to AD buffers
    xAD, yAD, VAD = get_newton_ad_values(storage.AD, x, y, V)
    update_flash_jacobian!(J, r, eos, cond.p, cond.T, cond.z, xAD, yAD, VAD, forces)
    return (J, r)
end

function update_secondary_jacobian!(storage, eos, c, x, y, V)
    c_ad = storage.AD_cond
    J = storage.J_inv
    p, T, z = get_inverse_ad_values(c_ad, c)
    forces = force_coefficients(eos, (p = p, T = T, z = z))
    update_flash_jacobian!(J, nothing, eos, p, T, z, x, y, V, forces)
    return (J, storage.r)
end

function update_and_solve!(storage, eos, cond, forces, x, y, V)
    J, r = update_primary_jacobian!(storage, eos, cond, forces, x, y, V)
    return solve_jacobian!(J, r)
end

function solve_jacobian!(J, r)
    # Factorize in-place, overwriting Jacobian with its factorization
    F = lu!(J)
    ldiv!(F, r)
end

function solve_jacobian!(J::MMatrix, r)
    M = length(r)
    tmp = SMatrix{M, M}(J)\SVector{M}(r)
    r .= tmp
end

function get_inverse_ad_values(cond_ad, cond)
    p, T, z = cond.p, cond.T, cond.z
    z_ad = cond_ad.z
    for i in eachindex(z)
        @inbounds z_ad[i] = update_value(z_ad[i], z[i])
    end
    p_ad = update_value(cond_ad.p, p)
    T_ad = update_value(cond_ad.T, T)
    return (p_ad, T_ad, z_ad)
end

newton_dampen!(::Nothing, Δ) = nothing

function newton_dampen!(dMax, Δ)
    ω = minimum((x) -> min(abs(dMax/x), 1.0), Δ)
    @. Δ *= ω
end

function flash_update!(K, storage, type::SSINewtonFlash, eos, cond, forces, V, iteration)
    if iteration >= type.swap_iter
        flash_update!(K, storage, NewtonFlash(type.dMax), eos, cond, forces, V, iteration)
    else
        flash_update!(K, storage, SSIFlash(), eos, cond, forces, V, iteration)
    end
end

function prep(eos, cond, forces)
    s = force_scalars(eos, cond, forces)
    Z = mixture_compressibility_factor(eos, cond, forces, s)
    return (Z, s)
end

function update_flash_jacobian!(J, r, eos, p, T, z, x, y, V, forces)
    has_r = !isnothing(r)
    n = number_of_components(eos)
    liquid = (p = p, T = T, z = x)
    vapor = (p = p, T = T, z = y)

    Z_l, s_l = prep(eos, liquid, forces)
    Z_v, s_v = prep(eos, vapor, forces)

    if isa(V, ForwardDiff.Dual)
        np = length(V.partials)
    else
        np = length(p.partials)
    end
    # Isofugacity constraint
    # f_li - f_vi ∀ i
    @inbounds for c in 1:n
        f_l = component_fugacity(eos, liquid, c, Z_l, forces, s_l)
        f_v = component_fugacity(eos, vapor, c, Z_v, forces, s_v)
        Δf = f_l - f_v
        if has_r
            r[c+n] = Δf.value
        end
        for i = 1:np
            J[c+n, i] = Δf.partials[i]
        end
    end
    # x_i*(1-V) - V*y_i - z_i = 0 ∀ i
    # Σ_i x_i - y_i = 0
    T = Base.promote_type(eltype(x), eltype(y))
    L = 1 - V
    Σxy = zero(T)
    @inbounds for c in 1:n
        xc, yc = x[c], y[c]
        Σxy += (xc - yc)
        M = L*xc + V*yc - z[c]
        if has_r
            r[c] = M.value
        end
        @inbounds for i = 1:np
            J[c, i] = ∂(M, i)
        end
    end
    if has_r
        @inbounds r[end] = Σxy.value
    end
    if isa(Σxy, ForwardDiff.Dual)
        @inbounds for i = 1:np
            J[end, i] = ∂(Σxy, i)
        end
    else
        @. J[end, :] = 0
    end
end

∂(D, i) = D.partials[i]

"Compute liquid mole fraction from overall mole fraction, K-value and V vapor fraction"
@inline liquid_mole_fraction(z, K, V) = z/(1 - V + V*K)
"Compute vapor mole fraction from overall mole fraction, K-value and V vapor fraction"
@inline vapor_mole_fraction(z, K, V) = K*liquid_mole_fraction(z, K, V)
"Compute vapor mole fraction from liquid mole fraction and K-value"
@inline vapor_mole_fraction(x, K) = x*K

liquid_mole_fraction!(x, z, K, V) = begin x .= z ./ (1 .- V .+ V .* K);x end
vapor_mole_fraction!(y, x, K) = begin y .= x .* K;y end
function vapor_mole_fraction!(y,z, K, V) 
    x = y
    liquid_mole_fraction!(x, z, K, V)
    vapor_mole_fraction!(y, x, K)
end
