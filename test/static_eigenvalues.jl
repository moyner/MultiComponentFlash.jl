using Test, StaticArrays, LinearAlgebra

function smallest_eigenvalue_allocations(A)
    return @allocated MultiComponentFlash.static_minimum_eigenvalue(A)
end

function check_static_smallest_eigenvalues(::Val{N}, ::Type{T}) where {N,T}
    solve = MultiComponentFlash.static_minimum_eigenvalue
    identity = SMatrix{N,N,T}(I)
    matrices = SMatrix{N,N,T,N*N}[zero(identity), identity, -identity]
    diagonal = T.(collect(1:N))
    push!(matrices, SMatrix{N,N,T}(Diagonal(diagonal)))
    for seed in 1:12
        X = T[sin(i*j + seed/3) + cos(i + 2j + seed) for i in 1:N, j in 1:N]
        A = (X + transpose(X))/2
        push!(matrices, SMatrix{N,N,T}(A), SMatrix{N,N,T}(X*transpose(X)))
    end
    v = T.(collect(1:N))
    singular = identity - SVector{N,T}(v)*transpose(SVector{N,T}(v))/sum(abs2,v)
    for shift in (zero(T), -8eps(T), 8eps(T))
        push!(matrices, singular + shift*identity)
    end
    for scale in (T(1), sqrt(floatmin(T)), sqrt(floatmax(T))/T(N),
                  T === Float64 ? T(1e-310) : T(1e-40))
        for A in matrices
            B = scale*A
            magnitude = maximum(abs,B)
            reference = magnitude == 0 ? zero(T) : eigmin(Symmetric(Matrix(B/magnitude)))
            value = @inferred solve(B)
            normalized = magnitude == 0 ? value : value/magnitude
            @test isfinite(value)
            norm = magnitude == 0 ? one(T) : opnorm(Matrix(B/magnitude),Inf)
            spacing = magnitude == 0 ? zero(T) : eps(value)/magnitude
            @test abs(normalized-reference) <= 64N*eps(T)*norm + 2spacing
            @test normalized <= reference + 4N*eps(T)*norm
        end
    end
    solve(singular)
    @test smallest_eigenvalue_allocations(singular) == 0
end

@testset "Static smallest symmetric eigenvalue" begin
    for N in (1,2,3,4,5,8,10,12), T in (Float32,Float64)
        check_static_smallest_eigenvalues(Val(N),T)
    end

    # A repeated zero eigenvalue and a negative eigenvalue hidden behind a zero
    # leading pivot must never become positive bypass distances.
    for T in (Float32,Float64)
        A = @SMatrix T[0 1 0 0; 1 0 0 0; 0 0 1 0; 0 0 0 1]
        @test MultiComponentFlash.static_minimum_eigenvalue(A) <= -one(T)
        @test MultiComponentFlash.static_minimum_eigenvalue(zero(A)) == zero(T)
        invalid = setindex(A,T(NaN),1,1)
        @test isnan(MultiComponentFlash.static_minimum_eigenvalue(invalid))
        invalid = setindex(A,T(Inf),1,1)
        @test isnan(MultiComponentFlash.static_minimum_eigenvalue(invalid))

        # Compare cancellation-prone 2x2 eigenvalues to an independent exact
        # input, high precision closed form rather than another Float eigensolve.
        for offset in (-2eps(T),zero(T),2eps(T))
            B = @SMatrix [one(T) one(T); one(T) one(T)+offset]
            setprecision(BigFloat,256) do
                a,b,d = BigFloat(B[1,1]),BigFloat(B[1,2]),BigFloat(B[2,2])
                reference = (a+d)/2 - sqrt(((a-d)/2)^2+b^2)
                @test BigFloat(MultiComponentFlash.static_minimum_eigenvalue(B)) <= reference
            end
        end
    end

    # Respect the upper-triangle convention of the previous Symmetric wrapper.
    A = @SMatrix [2.0 0.25; 999.0 1.0]
    @test MultiComponentFlash.static_minimum_eigenvalue(A) ≈ eigmin(Symmetric(Matrix(A))) atol=1e-13

    # Matrices outside the generated size range retain the CPU LAPACK fallback.
    n = MultiComponentFlash.MAX_STATIC_EIGENVALUE_SIZE + 1
    A = SMatrix{n,n,Float64}(Diagonal(collect(1.0:n)))
    @test MultiComponentFlash.static_minimum_eigenvalue(A) == 1.0
end
