using Test, MultiComponentFlash, StaticArrays, ForwardDiff, LinearAlgebra
using KernelAbstractions, JLArrays

@kernel function static_newton_kernel!(out, converged, eos, cond, K0, V0, tolerance)
    i = @index(Global)
    V, _, report = flash_2ph!(
        MultiComponentFlash.StaticConfig(), K0, eos, cond, V0;
        method = NewtonFlash(dMax = 1.0), maxiter = 3, extra_out = true,
        tolerance = tolerance, z_min = nothing, check = false
    )
    @inbounds out[i] = V
    @inbounds converged[i] = report.converged
end

static_newton_allocations(eos, cond, K, V) = @allocated flash_2ph!(
    MultiComponentFlash.StaticConfig(), K, eos, cond, V;
    method = NewtonFlash(dMax = 1.0), maxiter = 3, extra_out = true,
    z_min = nothing, check = false
)

@testset "Static Gibbs Newton flash" begin
    for F in (Float64, Float32)
        eos = make_eos_immutable(get_test_eos(); float_type = F)
        cond = (p = F(1.0e6), T = F(300), z = SVector{3, F}(0.5, 0.3, 0.2))
        reference_eos = make_eos_immutable(get_test_eos())
        V_reference, K_reference = flash_2ph_immutable(
            reference_eos,
            (p = Float64(cond.p), T = Float64(cond.T), z = SVector{3, Float64}(cond.z))
        )
        V0, K0 = F(V_reference), SVector{3, F}(K_reference)
        nearby = (
            p = F(1.0001e6), T = F(300.01),
            z = SVector{3, F}(0.50001, 0.29999, 0.2),
        )
        tolerance = F == Float32 ? 1.0e-4 : 1.0e-8
        V, K, report = flash_2ph!(
            MultiComponentFlash.StaticConfig(),
            K0, eos, nearby, V0; method = NewtonFlash(dMax = 1.0),
            maxiter = 3, extra_out = true, tolerance = tolerance,
            z_min = nothing, check = false
        )
        expected_V, expected_K = flash_2ph_immutable(
            reference_eos,
            (p = Float64(nearby.p), T = Float64(nearby.T), z = SVector{3, Float64}(nearby.z))
        )
        @test report.converged
        @test report.its <= 3
        @test V isa F
        @test K isa SVector{3, F}
        @test V ≈ expected_V atol = 2tolerance
        @test K ≈ expected_K rtol = 10tolerance
        x = nearby.z ./ (1 .+ V .* (K .- 1))
        y = K .* x
        @test sum(x) ≈ one(F) atol = 10eps(F)
        @test sum(y) ≈ one(F) atol = 10eps(F)
        @test (1 - V) * x + V * y ≈ nearby.z
        if isdefined(JLArrays, :JLBackend)
            out = JLArray(zeros(F, 8))
            converged = JLArray(falses(8))
            backend = JLArrays.JLBackend()
            static_newton_kernel!(backend, 8)(
                out, converged, eos, nearby,
                K0, V0, tolerance; ndrange = 8
            )
            @test all(Array(converged))
            @test all(isapprox.(Array(out), V; atol = tolerance, rtol = 0))
        end
        if F == Float64
            static_newton_allocations(eos, nearby, K0, V0)
            @test static_newton_allocations(eos, nearby, K0, V0) == 0
        end

        # A finite but stale V is normalized before evaluating fugacities.
        other_V, other_K, other_report = flash_2ph!(
            MultiComponentFlash.StaticConfig(),
            K0, eos, nearby, F(0.1); method = NewtonFlash(dMax = 1.0),
            maxiter = 3, extra_out = true, tolerance = tolerance,
            z_min = nothing, check = false
        )
        @test other_report.converged
        @test other_V ≈ V
        @test other_K ≈ K

        _, _, unfinished = flash_2ph!(
            MultiComponentFlash.StaticConfig(),
            K0, eos, nearby, V0; method = NewtonFlash(dMax = 1.0),
            maxiter = 0, extra_out = true, tolerance = 1.0e-12,
            z_min = nothing, check = false
        )
        @test !unfinished.converged
        _, _, trivial = flash_2ph!(
            MultiComponentFlash.StaticConfig(),
            one.(K0), eos, nearby, V0; method = NewtonFlash(),
            maxiter = 3, extra_out = true, z_min = nothing, check = false
        )
        @test !trivial.converged
    end
end

@testset "Gibbs Hessian from the two-phase derivative Jacobian" begin
    MCF = MultiComponentFlash
    eos = make_eos_immutable(get_test_eos())
    cond = (p = 1.0e6, T = 300.0, z = @SVector [0.5, 0.3, 0.2])
    V, K = flash_2ph_immutable(eos, cond)
    # Check away from the stationary solution so the residual and its
    # derivatives are both nonzero. The component amounts still sum to z.
    vapor = V .* K .* (cond.z ./ (1 .+ V .* (K .- 1)))
    vapor += @SVector [0.001, -0.0001, 0.0002]
    liquid = cond.z - vapor
    gradient, hessian, scale = MCF.static_newton_system(eos, cond, liquid, vapor)
    forces = MCF.static_force_coefficients(eos, cond, Float64)
    potentials(v) = MCF.static_newton_phase_potentials(eos, cond, v, Val(:vapor), forces) -
        MCF.static_newton_phase_potentials(eos, cond, cond.z - v, Val(:liquid), forces)
    independent_hessian = ForwardDiff.jacobian(potentials, vapor)
    @test gradient ≈ potentials(vapor) rtol = 1.0e-10 atol = 1.0e-12
    @test hessian ≈ Diagonal(scale) * independent_hessian * Diagonal(scale) rtol = 1.0e-10
    @test isposdef(Symmetric(hessian))
    @test scale .* (hessian \ (-scale .* gradient)) ≈
        independent_hessian \ (-gradient) rtol = 1.0e-10
end
