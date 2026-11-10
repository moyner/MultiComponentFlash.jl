using Test, StaticArrays, LinearAlgebra

function static_eigenvalue_allocations(B)
    MultiComponentFlash.static_minimum_eigenvalue(B)
    return @allocated MultiComponentFlash.static_minimum_eigenvalue(B)
end

@testset "Static minimum eigenvalue without LAPACK" begin
    for T in (Float32, Float64), n in (4, 10, 20)
        matrices = [zeros(T, n, n), Matrix{T}(I, n, n),
            Matrix(Diagonal(T.(range(-1, 2; length = n))))]
        for shift in 1:4
            A = T[sin(i*j + shift) for i in 1:n, j in 1:n]
            push!(matrices, A)
        end
        # Include clustered and nearly singular positive eigenvalues.
        Q = Matrix(qr(T[sin(i + n*j) for i in 1:n, j in 1:n]).Q)
        push!(matrices, Q*Diagonal(T.(range(eps(T), 1; length = n)))*Q')
        for A in matrices, scale in (T(1), sqrt(floatmin(T)), sqrt(floatmax(T)))
            B = SMatrix{n, n}(A*scale)
            distance = MultiComponentFlash.static_minimum_eigenvalue(B)
            expected = eigmin(Symmetric(A))*scale
            tolerance = 100*n*eps(T)
            @test isfinite(distance)
            @test isapprox(distance, expected;
                atol = tolerance*maximum(abs, B), rtol = tolerance)
            if VERSION >= v"1.12"
                @test static_eigenvalue_allocations(B) == 0
            end
        end
        invalid = SMatrix{n, n}(fill(T(NaN), n, n))
        @test isnan(MultiComponentFlash.static_minimum_eigenvalue(invalid))
    end
end
