using Test, MultiComponentFlash, StaticArrays, LinearAlgebra

@testset "One feed reference for both stability trials" begin
    # A phase-dependent model isolates the reference error from any
    # hydrogen correlation or reservoir input. Both self-trials used to report
    # trivial stability, despite a negative cross-model tangent-plane distance.
    mixture = MultiComponentMixture(
        ["Water", "Methane", "Ethane"];
        A_ij = [0.0 0.8 -0.8; 0.8 0.0 0.0; -0.8 0.0 0.0]
    )
    sw = SoreideWhitson(
        mixture; A = (0.0, 0.0, 0.0),
        A_mw = (0.0, 0.0, 0.0), alphas = (0.0, 0.0, 0.0)
    )
    eos = GenericCubicEOS(mixture, sw)
    immutable_eos = make_eos_immutable(eos)
    cond = (p = 7.5e6, T = 333.15, z = @SVector [1 - 2.0e-6, 1.0e-6, 1.0e-6])
    static_result = stability_2ph_immutable(immutable_eos, cond; z_min = nothing)
    mutable_cond = (p = cond.p, T = cond.T, z = collect(cond.z))
    mutable_stable, mutable_report = stability_2ph(
        eos, mutable_cond;
        extra_out = true
    )
    @test !static_result.stable
    @test !mutable_stable
    @test mutable_report.liquid.stable == static_result.report.liquid.stable
    @test mutable_report.vapor.stable == static_result.report.vapor.stable
    @test isnan(static_result.storage.critical_distance)

    # Independently evaluate the resulting trial against one feed tangent.
    MCF = MultiComponentFlash
    forces = MCF.static_force_coefficients(immutable_eos, cond, Float64)
    vapor_feed = (; cond..., phase = Val(:vapor))
    liquid_feed = (; cond..., phase = Val(:liquid))
    f_v = MCF.static_fugacities(immutable_eos, vapor_feed, forces, Float64)
    f_l = MCF.static_fugacities(immutable_eos, liquid_feed, forces, Float64)
    reference = dot(cond.z, log.(f_l ./ f_v)) > 0 ? f_v : f_l
    K0 = initial_guess_K(immutable_eos, cond, MCF.StaticConfig())
    distances = map((Val(true), Val(false))) do inside_is_vapor
        _, _, _, _, trial = MCF.static_michelsen_test(
            reference, cond.z, K0,
            immutable_eos, cond, forces, inside_is_vapor
        )
        trial_cond = (
            p = cond.p, T = cond.T, z = trial,
            phase = MCF.stability_phase(inside_is_vapor),
        )
        trial_f = MCF.static_fugacities(immutable_eos, trial_cond, forces, Float64)
        dot(trial, log.(trial_f ./ reference))
    end
    @test minimum(distances) < -1.0e-10
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
