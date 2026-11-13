# Bound generated code size as well as the GPU register footprint. Larger
# matrices retain the ordinary CPU eigensolver; GPU callers check this limit.
const MAX_STATIC_EIGENVALUE_SIZE = 12

@generated function eigenvalue_pack_upper(A::SMatrix{N,N,T}) where {N,T}
    entries = [:(A[$i, $j]) for j in 1:N for i in 1:j]
    return :(SVector{$(length(entries)),T}(($(entries...),)))
end

@generated function eigenvalue_trailing_matrix(A, ::Val{N}) where N
    entries = [:(A[$(j*(j-1)÷2+i)]) for j in 2:N for i in 2:j]
    return :(SVector{$(length(entries))}(($(entries...),)))
end

@generated function eigenvalue_packed_matvec(A, v::SVector{N,T}) where {N,T}
    entries = Any[]
    for i in 1:N
        terms = [:(A[$(max(i,j)*(max(i,j)-1)÷2+min(i,j))]*v[$j]) for j in 1:N]
        push!(entries, foldl((a,b) -> :($a + $b), terms))
    end
    return :(SVector{$N,T}(($(entries...),)))
end

@generated function eigenvalue_rank_two_update(A, v::SVector{N,T}, w) where {N,T}
    entries = [:(A[$(j*(j-1)÷2+i)] - (v[$i]*w[$j] + w[$i]*v[$j]))
        for j in 1:N for i in 1:j]
    return :(SVector{$(length(entries)),T}(($(entries...),)))
end

@generated function eigenvalue_packed_norm(A, ::Val{N}) where N
    entries = Any[]
    for i in 1:N
        terms = [:(abs(A[$(max(i,j)*(max(i,j)-1)÷2+min(i,j))])) for j in 1:N]
        push!(entries, foldl((a,b) -> :($a + $b), terms))
    end
    return :(maximum(SVector{$N}(($(entries...),))))
end

@inline function eigenvalue_householder_step(A::SVector{L,T}, ::Val{N}) where {L,T,N}
    M = N - 1
    x = SVector{M,T}(ntuple(i -> A[i*(i+1)÷2+1], Val(M)))
    trailing = eigenvalue_trailing_matrix(A, Val(N))
    xscale = maximum(abs, x)
    xscale == zero(T) && return trailing, zero(T)

    # Scale the reflector separately, including when its column is subnormal.
    u = x/xscale
    beta = -copysign(sqrt(sum(abs2, u)), u[1])
    divisor = u[1] - beta
    v = SVector{M,T}(ntuple(i -> i == 1 ? one(T) : u[i]/divisor, Val(M)))
    tau = (beta - u[1])/beta
    work = tau*eigenvalue_packed_matvec(trailing, v)
    w = work - (tau*dot(v, work)/2)*v
    return eigenvalue_rank_two_update(trailing, v, w), beta*xscale
end

@inline eigenvalue_tridiagonal(A::SVector{1,T}, ::Val{1}) where T =
    (SVector{1,T}(A[1]), SVector{0,T}())
@inline eigenvalue_tridiagonal(A::SVector{3,T}, ::Val{2}) where T =
    (SVector{2,T}(A[1], A[3]), SVector{1,T}(A[2]))

@generated function eigenvalue_tridiagonal(A::SVector{L,T}, ::Val{N}) where {L,T,N}
    # Keep the common small systems inline, but split larger reductions into
    # shorter compiled functions. Fully expanding N=12 Float32 triggers LLVM
    # optimizer/scheduler failures with Julia 1.13's LLVM 20 backend.
    diagonal = [i == 1 ? :(first_d) : :(rest_d[$(i-1)]) for i in 1:N]
    off_diagonal = [i == 1 ? :(first_off) : :(rest_e[$(i-1)]) for i in 1:N-1]
    body = quote
        trailing, first_off = eigenvalue_householder_step(A, Val(N))
        rest_d, rest_e = eigenvalue_tridiagonal(trailing, Val(N-1))
        first_d = A[1]
        d = SVector{$N,T}(($(diagonal...),))
        e = SVector{$(N-1),T}(($(off_diagonal...),))
        return d, e
    end
    pushfirst!(body.args, Expr(:meta, N <= 8 ? :inline : :noinline))
    return body
end

# For lambda_min we need only determine whether T-xI is positive definite.
# Stop at the first nonpositive LDL' pivot instead of counting the full inertia.
# The positive pivot floor also bounds e²/pivot and prevents overflow.
@inline function eigenvalue_below_spectrum(d::SVector{N,T}, e2, x, pivmin) where {N,T}
    pivot = d[1] - x
    pivot <= pivmin && return false
    @inbounds for i in 2:N
        pivot = (d[i] - x) - e2[i-1]/pivot
        pivot <= pivmin && return false
    end
    return true
end

@inline function eigenvalue_bisect(d::SVector{N,T}, e) where {N,T}
    e2 = map(abs2, e)
    maximum(e2) == zero(T) && return minimum(d)
    lower = T(Inf)
    upper = T(Inf)
    tnorm = zero(T)
    @inbounds for i in 1:N
        radius = (i > 1 ? abs(e[i-1]) : zero(T)) + (i < N ? abs(e[i]) : zero(T))
        lower = min(lower, d[i] - radius)
        upper = min(upper, d[i])
        tnorm = max(tnorm, abs(d[i]) + radius)
    end
    padding = T(4N)*eps(T)*max(tnorm, one(T))
    lower -= padding
    upper += padding
    tolerance = T(2)*eps(T)*max(tnorm, one(T))
    pivmin = floatmin(T)*max(one(T), maximum(e2))
    # Scaling gives an O(1) interval. 128 steps exceed the Float64 precision
    # needed for every supported size; stagnation is checked explicitly.
    for _ in 1:128
        mid = lower + (upper - lower)/2
        if upper - lower <= tolerance || mid == lower || mid == upper
            break
        end
        if eigenvalue_below_spectrum(d, e2, mid, pivmin)
            lower = mid
        else
            upper = mid
        end
    end
    return lower
end

"""
    static_minimum_eigenvalue(B::SMatrix)

Estimate the smallest eigenvalue of a real symmetric matrix, using its upper
triangle. Float32/Float64 matrices up to `MAX_STATIC_EIGENVALUE_SIZE` use packed
Householder reduction and a bracketed Sturm search without allocation or LAPACK.

The lower bracket endpoint is reduced by a norm-scaled roundoff allowance.
This deliberately favors disabling the flash bypass near a singularity over
overestimating its critical distance. It is not an interval-arithmetic proof.
On this static path, nonfinite input returns NaN. Larger matrices retain the
ordinary CPU fallback.
"""
@inline function static_minimum_eigenvalue(B::SMatrix{N,N,T}) where {N,T<:Union{Float32,Float64}}
    if N > MAX_STATIC_EIGENVALUE_SIZE
        return minimum(eigvals(Symmetric(B)))
    end
    packed = eigenvalue_pack_upper(B)
    scale = maximum(abs, packed)
    isfinite(scale) || return T(NaN)
    scale == zero(T) && return zero(T)
    N == 1 && return B[1,1]
    normalized = packed/scale
    norm = eigenvalue_packed_norm(normalized, Val(N))
    d, e = eigenvalue_tridiagonal(normalized, Val(N))
    value = eigenvalue_bisect(d, e)
    allowance = T(16N)*eps(T)*norm
    # De-normalization can round upward, particularly for subnormal outputs
    # whose spacing exceeds the normalized allowance. Take one representable
    # step down so that this final rounding cannot enlarge the bypass region.
    return prevfloat((value - allowance)*scale)
end

@inline static_minimum_eigenvalue(B::SMatrix) = minimum(eigvals(Symmetric(B)))
