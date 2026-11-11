using Test, MultiComponentFlash, StaticArrays, LinearAlgebra

@testset "SW aqueous/nonaqueous stability references" begin
    mixture = MultiComponentMixture(["Water", "Methane"])
    eos = GenericCubicEOS(mixture, SoreideWhitson(mixture))
    immutable_eos = make_eos_immutable(eos)
    # Water/methane tie line at 9 MPa and 333.15 K. Changing the overall
    # composition changes the phase amounts, not these phase compositions.
    water_liquid = 0.99890792980097
    water_vapor = 0.00475898719416
    for water in (0.001, 0.01, 0.8, 1 - 1e-10), static in (false, true)
        z = SVector(water, 1 - water)
        cond = (p = 9e6, T = 333.15, z = z)
        if static
            V, K, report = flash_2ph!(MultiComponentFlash.StaticConfig(),
                initial_guess_K(immutable_eos, cond, MultiComponentFlash.StaticConfig()),
                immutable_eos, cond, NaN; extra_out = true, z_min = nothing)
        else
            V, K, report = flash_2ph(eos,
                (p = cond.p, T = cond.T, z = collect(z));
                extra_out = true, z_min = nothing)
        end
        two_phase = water_vapor < water < water_liquid
        @test report.stability.stable == !two_phase
        if two_phase
            @test report.converged
            @test V ≈ (water_liquid - water)/(water_liquid - water_vapor) atol = 1e-9
            x = z ./ (1 .+ V .* (K .- 1))
            y = K .* x
            @test x[1] ≈ water_liquid atol = 1e-9
            @test y[1] ≈ water_vapor atol = 1e-9
            @test sum(x) ≈ 1 atol = 1e-12
            @test sum(y) ≈ 1 atol = 1e-12
            @test (1 - V)*x + V*y ≈ z atol = 1e-14
            forces = MultiComponentFlash.static_force_coefficients(immutable_eos, cond, Float64)
            fl = MultiComponentFlash.static_fugacities(immutable_eos,
                (; cond..., z = x, phase = Val(:liquid)), forces, Float64)
            fv = MultiComponentFlash.static_fugacities(immutable_eos,
                (; cond..., z = y, phase = Val(:vapor)), forces, Float64)
            @test maximum(abs, log.(fl ./ fv)) < 1e-8
        else
            @test isnan(V)
        end
    end
end

@testset "SW does not reverse the aqueous and gas phases" begin
    mixture = MultiComponentMixture(["Water", "Methane", "Hydrogen"])
    eos = GenericCubicEOS(mixture, SoreideWhitson(mixture))
    immutable_eos = make_eos_immutable(eos)
    cond = (p = 9e6, T = 320.0, z = @SVector [0.8, 0.18, 0.02])
    for static in (false, true)
        if static
            V, K, report = flash_2ph!(MultiComponentFlash.StaticConfig(),
                initial_guess_K(immutable_eos, cond, MultiComponentFlash.StaticConfig()),
                immutable_eos, cond, NaN; extra_out = true, z_min = nothing)
        else
            V, K, report = flash_2ph(eos,
                (p = cond.p, T = cond.T, z = collect(cond.z));
                extra_out = true, z_min = nothing)
        end
        x = cond.z ./ (1 .+ V .* (K .- 1))
        y = K .* x
        @test report.converged
        @test V ≈ 0.19966310869043 atol = 1e-9
        @test x[1] > 0.99
        @test y[1] < 0.01
    end
end

@testset "An unfinished mutable stability trial is not stable" begin
    eos = get_test_eos()
    cond = test_conditions()
    stable, report = stability_2ph(
        eos, cond, initial_guess_K(eos, cond),
        MultiComponentFlash.FlashConfig(print_output = false);
        maxiter = 1, extra_out = true
    )
    @test !stable
    @test !report.liquid.trivial
    @test !report.vapor.trivial
end
