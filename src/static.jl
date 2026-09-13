"""
    StaticConfig()

Marker returned by `flash_storage(...; static=true)`. It selects the immutable,
stack-oriented SSI implementation used in accelerator kernels.
"""
struct StaticConfig end

@inline print_output(::StaticConfig) = false
@inline use_dict_storage(::StaticConfig) = false

function flash_storage(eos::GenericCubicEOS, cond, method, config::StaticConfig; kwarg...)
    method isa SSIFlash || throw(ArgumentError("The static flash currently supports SSIFlash only."))
    return config
end

"""
    V, K = flash_2ph_immutable(eos, c[, storage]; <keyword arguments>)

Run the immutable, accelerator-friendly two-phase flash implementation.
`c.z` must be an `SVector`; `K` is returned as an `SVector` and `V` is the
scalar vapor fraction. When `storage` is omitted, a static storage marker is
created automatically.

Set `return_stability=true` to additionally return a
[`StaticStabilityResult`](@ref). Its immutable `storage` can be supplied as
`stability_storage` on the next call to enable the Michelsen bypass. The
default `bypass_tolerance=10` uses the same conservative pressure, temperature
and composition bounds as the mutable simulator integration.

The immutable path currently supports `SSIFlash` and generic cubic EOS values
converted with [`make_eos_immutable`](@ref).
"""
@inline function flash_2ph_immutable(eos, c;
        method = SSIFlash(),
        stability_storage = nothing,
        stability_bypass::Bool = !isnothing(stability_storage),
        return_stability::Bool = false,
        kwarg...)
    return flash_2ph_immutable(eos, c,
        flash_storage(eos, c; method = method, static = true);
        method = method,
        stability_storage = stability_storage,
        stability_bypass = stability_bypass,
        return_stability = return_stability,
        kwarg...)
end

@inline function flash_2ph_immutable(eos, c, storage::StaticConfig;
        method = SSIFlash(),
        stability_storage = nothing,
        stability_bypass::Bool = !isnothing(stability_storage),
        return_stability::Bool = false,
        kwarg...)
    c.z isa SVector || throw(ArgumentError(
        "flash_2ph_immutable requires c.z to be an SVector"))
    V, K, report = flash_2ph!(storage, initial_guess_K(eos, c, storage), eos, c,
        NaN; method = method, extra_out = true,
        stability_storage = stability_storage,
        stability_bypass = stability_bypass || return_stability,
        kwarg...)
    return immutable_flash_output(V, K, report.stability_result,
        Val(return_stability))
end

@inline immutable_flash_output(V, K, stability, ::Val{false}) = (V, K)
@inline immutable_flash_output(V, K, stability, ::Val{true}) = (V, K, stability)

"""Return an isbits representation of a mixture for accelerator kernels."""
function static_mixture(mixture::MultiComponentMixture{R, N}) where {R, N}
    names = ntuple(_ -> nothing, Val(N))
    bic = mixture.binary_interaction
    if !isnothing(bic)
        bic = SMatrix{N, N, R}(bic)
    end
    return MultiComponentMixture(mixture.properties; A_ij = bic, names = names, name = nothing)
end

"""
    make_eos_immutable(eos)

Convert a generic cubic EOS to an isbits representation for accelerator kernels.
"""
function make_eos_immutable(eos::GenericCubicEOS{T, R, N}) where {T, R, N}
    mixture = static_mixture(eos.mixture)
    volume_shift = eos.volume_shift
    if !isnothing(volume_shift)
        volume_shift = SVector{N, eltype(volume_shift)}(volume_shift)
    end
    return GenericCubicEOS(
        eos.type,
        mixture,
        eos.m_1,
        eos.m_2,
        eos.ω_a,
        eos.ω_b,
        volume_shift
    )
end

"""Return immutable Wilson K-values for static storage."""
@inline function initial_guess_K(eos::GenericCubicEOS{E, R, N}, cond,
        ::StaticConfig) where {E, R, N}
    T = Base.promote_eltype(cond.p, cond.T, cond.z[1])
    properties = eos.mixture.properties
    return SVector{N, T}(ntuple(i -> wilson_estimate(properties[i], cond.p, cond.T), Val(N)))
end

# Val phase tags keep Symbol construction and dynamic dispatch out of kernels.
@inline phase_symbol(::Val{phase}) where phase = phase

@inline function pick_root(eos, roots, cond, forces, scalars, ::Val{:liquid})
    min_root, _ = root_bounds(roots, minimum_allowable_root(eos, forces, scalars))
    return min_root
end

@inline function pick_root(eos, roots, cond, forces, scalars, ::Val{:vapor})
    _, max_root = root_bounds(roots, minimum_allowable_root(eos, forces, scalars))
    return max_root
end

"""Immutable force coefficients for accelerator kernels."""
@inline function static_force_coefficients(eos::GenericCubicEOS{E, R, N}, cond,
        ::Type{T}) where {E, R, N, T}
    A_i_static = SVector{N, T}(ntuple(i -> A_i(eos, cond, i), Val(N)))
    B_i_static = SVector{N, T}(ntuple(i -> B_i(eos, cond, i), Val(N)))
    A_ij_static = SMatrix{N, N, T}(ntuple(Val(N*N)) do index
        i = mod1(index, N)
        j = (index - 1) ÷ N + 1
        sqrt(A_i_static[i]*A_i_static[j]) *
            (one(T) - binary_interaction(eos, i, j, cond))
    end)
    return (A_ij = A_ij_static, A_i = A_i_static, B_i = B_i_static)
end

function make_eos_immutable(eos::KValuesEOS{T, R, N}) where {T, R, N}
    mixture = static_mixture(eos.mixture)
    evaluator = eos.K_values_evaluator
    if evaluator isa AbstractVector
        evaluator = SVector{N, eltype(evaluator)}(evaluator)
    end
    volume_shift = eos.volume_shift
    if !isnothing(volume_shift)
        volume_shift = SVector{N, eltype(volume_shift)}(volume_shift)
    end
    return KValuesEOS(evaluator, mixture; volume_shift = volume_shift)
end

@inline function static_condition(c, ::Type{F}, ::Val{N};
        z_min = nothing) where {F, N}
    z = SVector{N, F}(ntuple(Val(N)) do i
        isnothing(z_min) ? c.z[i] : max(c.z[i], z_min)
    end)
    return (p = convert(F, c.p), T = convert(F, c.T), z = z)
end

@inline invalid_stability_storage(cond) =
    StaticStabilityStorage(cond, convert(typeof(cond.p), NaN))

@inline function stability_bypass_available(storage::StaticStabilityStorage,
        cond; tolerance::Real = 10.0)
    tolerance > zero(tolerance) || throw(ArgumentError(
        "bypass_tolerance must be positive"))
    b = storage.critical_distance
    if !(isfinite(b) && b > zero(b))
        return false
    end
    reference = storage.reference
    return maximum(abs, reference.z - cond.z) < b/tolerance &&
        abs(reference.p - cond.p) < b*abs(cond.p)/tolerance &&
        abs(reference.T - cond.T) < b*tolerance
end

@inline static_minimum_eigenvalue(B::SMatrix) =
    minimum(eigvals(Symmetric(B)))

"""Immutable Michelsen critical-point distance used by the stability bypass."""
@generated function static_fugacity_coefficients(
        eos::GenericCubicEOS{E, R, N}, cond, Z, forces, scalars,
        ::Type{D}) where {E, R, N, D}
    values = [:(component_fugacity_coefficient(
        eos, cond, $i, Z, forces, scalars)) for i in 1:N]
    return :(SVector{N, D}(($(values...),)))
end

@inline function static_michelsen_critical_point_measure(
        eos::GenericCubicEOS{E, R, N}, p, temperature,
        mole_numbers::SVector{N, F}) where {E, R, N, F}
    D = ForwardDiff.Dual{Nothing, F, N}
    mole_numbers_ad = SVector{N, D}(ntuple(Val(N)) do i
        partials = ForwardDiff.single_seed(ForwardDiff.Partials{N, F}, Val(i))
        D(mole_numbers[i], partials)
    end)
    z = mole_numbers_ad/sum(mole_numbers_ad)
    cond = (p = convert(F, p), T = convert(F, temperature), z = z)
    forces = static_force_coefficients(eos, cond, D)
    scalars = force_scalars(eos, cond, forces)
    Z = mixture_compressibility_factor(eos, cond, forces, scalars)
    coefficients = static_fugacity_coefficients(
        eos, cond, Z, forces, scalars, D)
    B = SMatrix{N, N, F}(ntuple(Val(N*N)) do index
        i = mod1(index, N)
        j = (index - 1) ÷ N + 1
        F(i == j) + sqrt(mole_numbers[i]*mole_numbers[j])*coefficients[i].partials[j]
    end)
    # Roundoff in the AD construction can make the theoretically symmetric
    # matrix differ by a few ulps. Symmetrize before finding its eigenvalues.
    B = (B + transpose(B))/2
    return static_minimum_eigenvalue(B)
end

@inline function solve_rachford_rice(K::StaticVector{2}, z::StaticVector{2}, V = NaN)
    z1, z2 = z
    k1, k2 = K
    b1, b2 = inv(1 - k1), inv(1 - k2)
    return (z1*b2 + z2*b1)/(z1 + z2)
end

@inline function solve_rachford_rice(K::StaticVector{3}, z::StaticVector{3}, V = NaN)
    z1, z2, z3 = z
    k1, k2, k3 = K
    b1, b2, b3 = inv(1-k1), inv(1-k2), inv(1-k3)
    a2 = z1 + z2 + z3
    a1 = -b1*(z2 + z3) - b2*(z1 + z3) - b3*(z1 + z2)
    a0 = b1*b2*z3 + b1*b3*z2 + b2*b3*z1
    discriminant = a1*a1 - 4*a0*a2
    if discriminant >= zero(discriminant)
        inv_2a2 = inv(2*a2)
        root_offset = sqrt(discriminant)*inv_2a2
        root_center = -a1*inv_2a2
        root1 = root_center - root_offset
        root2 = root_center + root_offset
        if zero(root1) < root1 < one(root1)
            return root1
        elseif zero(root2) < root2 < one(root2)
            return root2
        elseif isfinite(root1 + root2)
            kmin = min(k1, k2, k3)
            kmax = max(k1, k2, k3)
            kmin > one(kmin) && return max(root1, root2)
            kmax < one(kmax) && return min(root1, root2)
        end
    end
    return solve_rachford_rice_static_iterative(K, z, V)
end

@inline solve_rachford_rice(K::StaticVector, z::StaticVector, V = NaN) =
    solve_rachford_rice_static_iterative(K, z, V)

@inline function solve_rachford_rice_static_iterative(K, z, V;
        tol = 1e-12, maxiter = 1000)
    V_lo = inv(1 - maximum(K))
    V_hi = inv(1 - minimum(K))
    if V_hi < V_lo
        V_lo, V_hi = V_hi, V_lo
    end
    if isnan(V)
        V = (V_lo + V_hi)/2
    end
    for _ in 1:maxiter
        residual = zero(V)
        denominator = zero(V)
        @inbounds for i in eachindex(K)
            delta_K = K[i] - one(K[i])
            term_denominator = one(V) + V*delta_K
            residual += z[i]*delta_K/term_denominator
            denominator += z[i]*delta_K^2/term_denominator^2
        end
        abs(residual) < tol && break
        if residual > zero(residual)
            V_lo = V
        else
            V_hi = V
        end
        V_next = V + residual/denominator
        if !(V_lo < V_next < V_hi) || !isfinite(V_next)
            V_next = (V_lo + V_hi)/2
        end
        V = V_next
    end
    return V
end

@generated function static_fugacities(eos::GenericCubicEOS{E, R, N}, cond, forces,
        ::Type{F}) where {E, R, N, F}
    values = [:(component_fugacity(
        eos, cond, $i, Z, forces, scalars)) for i in 1:N]
    return quote
        Z, scalars = prep(eos, cond, forces)
        SVector{N, F}(($(values...),))
    end
end

@inline function static_ssi(K::SVector{N, F}, p::F, T::F, z, V::F,
        eos, forces) where {N, F<:Real}
    x = SVector{N, F}(ntuple(i -> liquid_mole_fraction(z[i], K[i], V), Val(N)))
    y = SVector{N, F}(ntuple(i -> vapor_mole_fraction(x[i], K[i]), Val(N)))
    liquid = (p = p, T = T, z = x, phase = Val(:liquid))
    vapor = (p = p, T = T, z = y, phase = Val(:vapor))
    f_l = static_fugacities(eos, liquid, forces, F)
    f_v = static_fugacities(eos, vapor, forces, F)
    ratios = SVector{N, F}(ntuple(i -> f_l[i]/f_v[i], Val(N)))
    residual = zero(F)
    @inbounds for i in 1:N
        residual = max(residual, abs(one(F) - ratios[i]))
    end
    K_next = SVector{N, F}(ntuple(i -> K[i]*ratios[i], Val(N)))
    V_next = solve_rachford_rice(K_next, z, V)
    V_next = clamp(V_next, zero(V_next), one(V_next))
    return V_next, K_next, residual
end

@inline function flash_2ph(eos::GenericCubicEOS, c, K, V,
        config::StaticConfig; kwarg...)
    return flash_2ph!(config, K, eos, c, V; kwarg...)
end

@inline function flash_2ph!(config::StaticConfig, K, eos::GenericCubicEOS, c,
        V = NaN; extra_out::Bool = false, kwarg...)
    out = flash_2ph_impl!(config, K, eos, c, V; kwarg...)
    return static_flash_output(out, Val(extra_out))
end

@inline static_flash_output(out, ::Val{true}) = out
@inline static_flash_output(out, ::Val{false}) = out[1]

@inline function flash_2ph_impl!(::StaticConfig, K,
        eos::GenericCubicEOS{E, R, N}, c, V;
        method::SSIFlash = SSIFlash(),
        maxiter::Int = 25000,
        tolerance::Float64 = 1e-8,
        verbose::Bool = false,
        check::Bool = true,
        update_forces::Bool = true,
        z_min = MINIMUM_COMPOSITION,
        stability_storage = nothing,
        stability_bypass::Bool = !isnothing(stability_storage),
        bypass_tolerance::Real = 10.0,
        kwarg...
    ) where {E, R, N}
    F = Base.promote_eltype(c.p, c.T, c.z[1], K[1])
    K = SVector{N, F}(K)
    cond = static_condition(c, F, Val(N); z_min = z_min)
    z = cond.z
    forces = static_force_coefficients(eos, cond, F)
    V = convert(F, V)
    single_phase_init = isnan(V) || V == one(F) || V == zero(F)
    if single_phase_init
        stability_result = static_stability_2ph(K, eos, cond, forces;
            storage = stability_storage_value(stability_storage),
            update_bypass = stability_bypass,
            bypass_tolerance = bypass_tolerance,
            maxiter = maxiter,
            kwarg...)
        stable = stability_result.stable
        stability_report = stability_result.report
        K = stability_result.K
    else
        stable = false
        stability_report = StabilityReport(false, false, false, false)
        stability_result = StaticStabilityResult(stable, stability_report, K,
            invalid_stability_storage(cond), false)
    end
    converged = false
    if stable
        iteration = 0
    else
        iteration = 1
        if isnan(V)
            V = solve_rachford_rice(K, z, V)
        end
        while true
            V, K, residual = static_ssi(K, cond.p, cond.T, z, V, eos, forces)
            converged = residual <= tolerance
            (converged || iteration == maxiter) && break
            iteration += 1
        end
    end
    report = (its = iteration, converged = converged,
        stability = stability_report, stability_result = stability_result)
    return V, K, report
end

stability_phase(::Val{true}) = Val(:vapor)
stability_phase(::Val{false}) = Val(:liquid)

@generated function stability_xy(z::SVector{N, F}, K::SVector{N, F}, phase) where {N, F}
    values = [:(xy_value(z[$i], K[$i], phase)) for i in 1:N]
    return :(SVector{N, F}(($(values...),)))
end

@generated function static_scale(v::SVector{N, F}, scale::F) where {N, F}
    values = [:(v[$i]/scale) for i in 1:N]
    return :(SVector{N, F}(($(values...),)))
end

@generated function stability_ratios(f_z::SVector{N, F}, f_xy::SVector{N, F},
        scale::F, phase) where {N, F}
    values = [:(f_ratio(f_z[$i], scale*f_xy[$i], phase)) for i in 1:N]
    return :(SVector{N, F}(($(values...),)))
end

@generated function static_multiply(a::SVector{N, F}, b::SVector{N, F}) where {N, F}
    values = [:(a[$i]*b[$i]) for i in 1:N]
    return :(SVector{N, F}(($(values...),)))
end

@generated function static_divide(a::SVector{N, F}, b::SVector{N, F}) where {N, F}
    values = [:(a[$i]/b[$i]) for i in 1:N]
    return :(SVector{N, F}(($(values...),)))
end

@inline function static_michelsen_test(f_z, z::SVector{N, F}, K::SVector{N, F},
        eos, cond, forces, inside_is_vapor;
        tol_equil = 1e-10,
        tol_trivial = tol_equil,
        tol_sat = tol_trivial,
        maxiter = 1000
    ) where {N, F}
    trivial = false
    S = one(F)
    iter = 0
    xy = zero(SVector{N, F})
    while true
        iter += 1
        unnormalized = stability_xy(z, K, inside_is_vapor)
        S = zero(F)
        @inbounds for i in 1:N
            S += unnormalized[i]
        end
        xy = static_scale(unnormalized, S)
        inside = (p = cond.p, T = cond.T, z = xy,
            phase = stability_phase(inside_is_vapor))
        f_xy = static_fugacities(eos, inside, forces, F)
        ratios = stability_ratios(f_z, f_xy, S, inside_is_vapor)
        K = static_multiply(K, ratios)
        R_norm = zero(F)
        K_norm = zero(F)
        @inbounds for i in 1:N
            R_norm += (ratios[i] - one(F))^2
            K_norm += log(K[i])^2
        end
        trivial = K_norm < tol_trivial
        converged = R_norm < tol_equil
        if trivial || converged
            break
        elseif iter == maxiter
            trivial = true
            break
        end
    end
    stable = trivial || S <= one(F) + tol_sat
    return stable, trivial, iter, K, xy
end

@inline function static_stability_2ph(K::SVector{N, F}, eos, cond, forces;
        check_vapor::Bool = true,
        check_liquid::Bool = true,
        storage = nothing,
        update_bypass::Bool = false,
        bypass_tolerance::Real = 10.0,
        kwarg...
    ) where {N, F}
    if update_bypass && !isnothing(storage) &&
            stability_bypass_available(storage, cond;
                tolerance = bypass_tolerance)
        report = StabilityReport(true, true, true, true)
        return StaticStabilityResult(true, report, K, storage, true)
    end
    vapor_phase = (p = cond.p, T = cond.T, z = cond.z, phase = Val(:vapor))
    f_z_vapor = static_fugacities(eos, vapor_phase, forces, F)
    K_wilson = initial_guess_K(eos, cond, StaticConfig())
    if check_vapor
        stable_vapor, trivial_vapor, i_v, K_vapor, y = static_michelsen_test(
            f_z_vapor, cond.z, K_wilson, eos, cond, forces, Val(true); kwarg...)
    else
        stable_vapor, trivial_vapor, i_v, K_vapor, y = true, true, 0, K_wilson, cond.z
    end
    if check_liquid
        liquid_phase = (p = cond.p, T = cond.T, z = cond.z, phase = Val(:liquid))
        f_z_liquid = forces_per_phase(eos) ?
            static_fugacities(eos, liquid_phase, forces, F) : f_z_vapor
        stable_liquid, trivial_liquid, i_l, K_liquid, x = static_michelsen_test(
            f_z_liquid, cond.z, K_wilson, eos, cond, forces, Val(false); kwarg...)
    else
        stable_liquid, trivial_liquid, i_l, K_liquid, x = true, true, 0, K_wilson, cond.z
    end
    report = StabilityReport(stable_liquid, trivial_liquid,
        stable_vapor, trivial_vapor)
    K_out = report.stable ? K_liquid : static_divide(y, x)
    if update_bypass && report.stable && report.liquid.trivial &&
            report.vapor.trivial
        critical_distance = static_michelsen_critical_point_measure(
            eos, cond.p, cond.T, cond.z)
        next_storage = StaticStabilityStorage(cond, critical_distance)
    else
        next_storage = invalid_stability_storage(cond)
    end
    return StaticStabilityResult(report.stable, report, K_out,
        next_storage, false)
end

@inline function stability_2ph(eos::GenericCubicEOS{E, R, N}, c, K,
        config::StaticConfig; extra_out::Bool = false, kwarg...) where {E, R, N}
    F = Base.promote_eltype(c.p, c.T, c.z[1], K[1])
    cond = static_condition(c, F, Val(N))
    K = SVector{N, F}(K)
    forces = static_force_coefficients(eos, cond, F)
    result = static_stability_2ph(K, eos, cond, forces; kwarg...)
    return extra_out ? (result.stable, result.report) : result.stable
end

@inline stability_2ph!(::StaticConfig, K, eos::GenericCubicEOS, c; kwarg...) =
    stability_2ph(eos, c, K, StaticConfig(); kwarg...)

"""
    result = stability_2ph_immutable(eos, c[, storage]; <keyword arguments>)

Run the immutable stability test independently of a flash. The returned
[`StaticStabilityResult`](@ref) contains the stability report, K-values and
updated [`StaticStabilityStorage`](@ref). Pass `result.storage` to a later call
to enable Michelsen's stability bypass for nearby conditions.

The bypass is only armed after both trial phases converge to trivial stable
solutions. Calls inside the shadow region retain the full stability test.
"""
@inline function stability_2ph_immutable(eos::GenericCubicEOS{E, R, N}, c,
        storage = nothing;
        K = initial_guess_K(eos, c, StaticConfig()),
        bypass_tolerance::Real = 10.0,
        z_min = MINIMUM_COMPOSITION,
        kwarg...) where {E, R, N}
    c.z isa SVector || throw(ArgumentError(
        "stability_2ph_immutable requires c.z to be an SVector"))
    F = Base.promote_eltype(c.p, c.T, c.z[1], K[1])
    cond = static_condition(c, F, Val(N); z_min = z_min)
    K = SVector{N, F}(K)
    forces = static_force_coefficients(eos, cond, F)
    return static_stability_2ph(K, eos, cond, forces;
        storage = stability_storage_value(storage),
        update_bypass = true,
        bypass_tolerance = bypass_tolerance,
        kwarg...)
end

@inline stability_storage_value(::Nothing) = nothing
@inline stability_storage_value(storage::StaticStabilityStorage) = storage
@inline stability_storage_value(result::StaticStabilityResult) = result.storage
