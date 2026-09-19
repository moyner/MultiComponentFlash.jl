include("test_setup.jl")
include("eos_equations.jl")

@testset "Rachford-Rice" begin
    test_rachford_rice()
end

flash_methods = [SSIFlash(), NewtonFlash(), SSINewtonFlash()]
@testset "Peng-Robinson" begin
    for m in flash_methods
        name = typeof(m)
        @testset "$name" begin
            test_flash_pr(m)
        end
    end
end
@testset "Soave-Redlich-Kwong" begin
    for m in flash_methods
        name = typeof(m)
        @testset "$name" begin
            test_flash_srk(m)
        end
    end
end
@testset "Redlich-Kwong" begin
    for m in flash_methods
        name = typeof(m)
        @testset "$name" begin
            test_flash_rk(m)
        end
    end
end
@testset "Zudkevitch-Joffe" begin
    for m in flash_methods
        name = typeof(m)
        @testset "$name" begin
            test_flash_zj(m)
        end
    end
end
@testset "Zero allocating flash" begin
    test_flash_inplace(SSIFlash())
end

@testset "Negative flash from V=Inf" begin
    host_eos = get_test_eos()
    static_eos = make_eos_immutable(host_eos)
    z = @SVector [0.5, 0.3, 0.2]
    static_config = MultiComponentFlash.StaticConfig()

    # These states are stable according to the ordinary flash, but their
    # extrapolated equilibrium has V below zero or above one, respectively.
    for (p, T, below_zero) in ((1e7, 250.0, true), (1e5, 400.0, false))
        cond = (p = p, T = T, z = collect(z))
        V_stable, _, stable_report = flash_2ph(host_eos, cond;
            extra_out = true)
        @test isnan(V_stable)
        @test stable_report.stability.stable

        static_cond = (p = p, T = T, z = z)
        K0 = initial_guess_K(static_eos, static_cond, static_config)
        V_static, K_static, static_report = flash_2ph!(static_config,
            K0, static_eos, static_cond, Inf; extra_out = true)
        @test static_report.converged
        @test !static_report.stability_result.stable
        @test (below_zero ? V_static < 0 : V_static > 1)
        V_immutable, K_immutable = flash_2ph_immutable(
            static_eos, static_cond, Inf)
        @test V_immutable ≈ V_static
        @test K_immutable ≈ K_static

        for method in (SSIFlash(), NewtonFlash(), SSINewtonFlash())
            host_cond = (p = p, T = T, z = collect(z))
            V, K, report = flash_2ph(host_eos, host_cond,
                initial_guess_K(host_eos, host_cond), Inf;
                method = method, extra_out = true)
            @test report.converged
            @test !report.stability.stable
            @test V ≈ V_static rtol = 1e-6
            @test K ≈ K_static rtol = 1e-6
            x = liquid_mole_fraction.(z, K, V)
            y = vapor_mole_fraction.(x, K)
            @test all(>(0), x)
            @test all(>(0), y)
            @test sum(x) ≈ 1 atol = 1e-7
            @test sum(y) ≈ 1 atol = 1e-7
        end
    end

    # A trivial K ≈ 1 fixed point is not a converged negative flash.
    trivial_cond = (p = 5e7, T = 250.0, z = collect(z))
    _, _, trivial_report = flash_2ph(host_eos, trivial_cond,
        initial_guess_K(host_eos, trivial_cond), Inf;
        extra_out = true, check = false, maxiter = 100)
    @test !trivial_report.converged

    K4 = @SVector [0.2, 0.4, 2.0, 4.0]
    z4 = @SVector [0.25, 0.25, 0.25, 0.25]
    @test solve_rachford_rice(K4, z4, Inf) ≈
        solve_rachford_rice(K4, z4, NaN)
    @test solve_rachford_rice(collect(K4), collect(z4), Inf) ≈
        solve_rachford_rice(collect(K4), collect(z4), NaN)

    # A guess at a pole must be reinitialized inside the positive-composition
    # window, even when the correct negative-flash root lies outside [0, 1].
    for (z_negative, below_zero) in
            ((@SVector([0.7, 0.2, 0.05, 0.05]), true),
             (@SVector([0.05, 0.05, 0.2, 0.7]), false))
        for (K_test, z_test) in ((K4, z_negative),
                (collect(K4), collect(z_negative)))
            V_expected = solve_rachford_rice(K_test, z_test)
            @test below_zero ? V_expected < 0 : V_expected > 1
            V_pole = 1/(1 - maximum(K_test))
            @test solve_rachford_rice(K_test, z_test, V_pole) ≈ V_expected
        end
    end

    # Regression for a phase-diagram failure: all K-values were above one,
    # yet the old RR solve returned a root between two negative poles. That
    # root gave negative/huge phase compositions and crashed the EOS.
    K_runaway = @SVector [1.0496438050320456, 1.8228323947565352,
        3.6255787603685725, 8.673724738551629,
        23.930790867112876, 63.42734611565452]
    z_runaway = @SVector [0.635, 0.115, 0.05, 0.1, 0.075, 0.025]
    V_pole = -0.016018621040647378
    @test isnan(solve_rachford_rice(K_runaway, z_runaway, V_pole))
    @test isnan(solve_rachford_rice(collect(K_runaway), collect(z_runaway), V_pole))
    @test isnan(solve_rachford_rice(K_runaway, z_runaway))

    # With K almost equal to one, a valid negative flash can have |V| >> 1.
    # The RR stopping test and composition formula must both preserve the
    # normalization of the extrapolated liquid and vapor phases.
    K_near = @SVector [0.1, 0.5, 1.0 + 1e-14, 1.0 + 2e-14]
    z_near = @SVector [0.01, 0.01, 0.49, 0.49]
    for (K_test, z_test) in ((K_near, z_near),
            (collect(K_near), collect(z_near)))
        V_near = solve_rachford_rice(K_test, z_test)
        x_near = liquid_mole_fraction.(z_test, K_test, V_near)
        y_near = vapor_mole_fraction.(x_near, K_test)
        @test V_near < 0
        @test all(>(0), x_near)
        @test sum(x_near) ≈ 1 atol = 1e-10
        @test sum(y_near) ≈ 1 atol = 1e-10
    end

    K_no_split = @SVector [1.1, 2.0, 3.0]
    V_bad, _, bad_report = flash_2ph!(static_config, K_no_split,
        static_eos, (p = 1e7, T = 250.0, z = z), Inf; extra_out = true)
    @test isnan(V_bad)
    @test !bad_report.converged
    V_bad_host, _, bad_host_report = flash_2ph(host_eos,
        (p = 1e7, T = 250.0, z = collect(z)), collect(K_no_split), Inf;
        extra_out = true)
    @test isnan(V_bad_host)
    @test !bad_host_report.converged
end

@testset "Static accelerator path" begin
    host_eos = get_test_eos()
    eos = make_eos_immutable(host_eos)
    c = (p = 1e6, T = 300.0, z = @SVector [0.5, 0.3, 0.2])
    storage = flash_storage(eos, c; method = SSIFlash(), static = true)
    K = initial_guess_K(eos, c, storage)
    V, K, report = flash_2ph!(storage, K, eos, c, 0.5; extra_out = true)

    @test isbitstype(typeof(eos))
    @test isbitstype(typeof(storage))
    @test storage isa MultiComponentFlash.StaticConfig
    @test report.converged
    @test report.its == 6
    @test V ≈ 0.7632068334421974
    @test K ≈ @SVector [4.553402802323027, 17.73895830809456, 0.0004031451448211194]

    normal_config = MultiComponentFlash.FlashConfig(print_output=false)
    @test typeof(normal_config) == MultiComponentFlash.FlashConfig
    @test flash_storage(eos, c, SSIFlash(), normal_config).x isa Vector
    deprecated_storage = @test_deprecated flash_storage(eos, c; static_size = true)
    @test deprecated_storage isa MultiComponentFlash.StaticConfig
    dynamic_storage = @test_deprecated flash_storage(eos, c; static_size = false)
    @test dynamic_storage.x isa Vector

    config = MultiComponentFlash.FlashConfig(print_output=false, use_dict_storage=false)
    @test typeof(config) == MultiComponentFlash.FlashConfig
    @test flash_storage(eos, c, SSIFlash(), config) isa MultiComponentFlash.StaticConfig
    K_config = initial_guess_K(eos, c, storage)
    V_config, K_config, config_report = flash_2ph(eos, c, K_config, NaN, config;
        method=SSIFlash(), extra_out=true, z_min=nothing)
    @test config_report.stability isa MultiComponentFlash.StabilityReport
    @test config_report.converged
    @test V_config ≈ V
    @test K_config ≈ K

    V_immutable, K_immutable = flash_2ph_immutable(eos, c)
    @test V_immutable ≈ V
    @test K_immutable ≈ K
    @test K_immutable isa SVector{3, Float64}
    @test flash_2ph_immutable(eos, c, storage) == (V_immutable, K_immutable)
    @test_throws ArgumentError flash_2ph_immutable(eos,
        (p = c.p, T = c.T, z = collect(c.z)))

    @testset "Standalone stability and bypass" begin
        stable_cond = (p = 1e5, T = 800.0, z = c.z)
        stability = stability_2ph_immutable(eos, stable_cond)
        @test stability isa MultiComponentFlash.StaticStabilityResult
        @test stability.storage isa MultiComponentFlash.StaticStabilityStorage
        @test isbitstype(typeof(stability))
        @test stability.stable
        @test stability.report.liquid.trivial
        @test stability.report.vapor.trivial
        @test isfinite(stability.storage.critical_distance)
        @test !stability.bypassed
        @test test_static_stability_allocs(eos, stable_cond) == 0

        nearby = (p = 1.001e5, T = 800.01,
            z = @SVector [0.50001, 0.29999, 0.2])
        bypassed = stability_2ph_immutable(eos, nearby, stability)
        @test bypassed.stable
        @test bypassed.bypassed
        @test bypassed.storage == stability.storage

        far_away = (p = 1e6, T = 800.0, z = c.z)
        retested = stability_2ph_immutable(eos, far_away, stability.storage)
        @test !retested.bypassed
        @test retested.storage.reference == far_away

        # A stable state inside the shadow region is deliberately not armed.
        shadow = stability_2ph_immutable(eos,
            (p = 1e5, T = 500.0, z = c.z))
        @test shadow.stable
        @test !shadow.report.liquid.trivial
        @test isnan(shadow.storage.critical_distance)

        V_stable, K_stable, flash_stability = flash_2ph_immutable(
            eos, stable_cond; return_stability = true)
        @test isnan(V_stable)
        @test all(isfinite, K_stable)
        @test flash_stability.stable
        V_nearby, K_nearby, nearby_stability = flash_2ph_immutable(
            eos, nearby;
            stability_storage = flash_stability,
            return_stability = true)
        @test isnan(V_nearby)
        @test all(isfinite, K_nearby)
        @test nearby_stability.bypassed
        @test_throws ArgumentError stability_2ph_immutable(eos, nearby,
            stability.storage; bypass_tolerance = 0.0)
    end

    @testset "Nearly absent components" begin
        tiny = 1e-16
        compositions = (
            SVector(tiny, 0.3, 0.7 - tiny),
            SVector(0.5, tiny, 0.5 - tiny),
            SVector(0.3, 0.7 - tiny, tiny),
            SVector(tiny, tiny, 1.0 - 2tiny)
        )
        for z_tiny in compositions
            cond_tiny = (p = 1e5, T = 250.0, z = z_tiny)
            result = stability_2ph_immutable(eos, cond_tiny)
            @test isbitstype(typeof(result))
            @test all(isfinite, result.K)

            V_tiny, K_tiny, report_tiny = flash_2ph!(storage,
                initial_guess_K(eos, cond_tiny, storage), eos, cond_tiny, NaN;
                extra_out = true)
            @test report_tiny.converged || result.stable
            @test isnan(V_tiny) || 0.0 <= V_tiny <= 1.0
            @test all(isfinite, K_tiny)

            host_cond_tiny = (p = cond_tiny.p, T = cond_tiny.T,
                z = collect(z_tiny))
            V_host, K_host, host_report = flash_2ph(host_eos,
                host_cond_tiny; extra_out = true)
            @test host_report.converged || host_report.stability.stable
            @test isnan(V_host) || 0.0 <= V_host <= 1.0
            @test all(isfinite, K_host)
            @test isapprox(V_tiny, V_host; nans = true, atol = 1e-12)

            # Exercise the AD/eigenvalue part of the bypass without flooring
            # away the trace component.
            distance = MultiComponentFlash.michelsen_critical_point_measure(
                eos, cond_tiny.p, cond_tiny.T, z_tiny)
            @test isfinite(distance)

            stable_tiny = (p = 5e7, T = 800.0, z = z_tiny)
            tiny_stability = stability_2ph_immutable(eos, stable_tiny;
                z_min = nothing)
            @test tiny_stability.stable
            @test isfinite(tiny_stability.storage.critical_distance)
            tiny_nearby = (p = 5.001e7, T = 800.01, z = z_tiny)
            tiny_bypass = stability_2ph_immutable(eos, tiny_nearby,
                tiny_stability; z_min = nothing)
            @test tiny_bypass.bypassed
        end
    end
end

@testset "Partial derivatives" begin
    test_flash_partials()
end

@testset "Constructors" begin
    m1 = MolecularProperty(0.0440, 7.38e6, 304.1, 9.412e-5, 0.224)
    m2 = MolecularProperty(mw = 0.0440, p_c = 7.38e6, T_c = 304.1, V_c = 9.412e-5, acentric_factor = 0.224)
    @test m1 == m2
end

@testset "K-value EOS" begin
    mixture = MultiComponentMixture(["CarbonDioxide", "Water"])
    eos = KValuesEOS([0.01, 2.0], mixture)
    cond = (p = 1e5, T = 273.15, z = (0.1, 0.9))
    @test round(flash_2ph(eos, cond), digits = 4) ≈ 0.8091
    @test number_of_components(eos) == 2
    eos2 = KValuesEOS(cond -> [0.01, 2.0], mixture)
    @test round(flash_2ph(eos2, cond), digits = 4) ≈ 0.8091

    static_eos = make_eos_immutable(eos)
    static_cond = (p = cond.p, T = cond.T, z = @SVector [0.1, 0.9])
    @test isbitstype(typeof(static_eos))
    @test static_eos.K_values_evaluator isa SVector{2, Float64}
    @test round(flash_2ph(static_eos, static_cond), digits = 4) ≈ 0.8091
end

@testset "Static flashed mixture storage" begin
    x = @SVector [0.8, 0.2]
    y = @SVector [0.1, 0.9]
    flashed = FlashedMixture2Phase(
        MultiComponentFlash.two_phase_lv, SVector(0.125, 4.5),
        0.4, x, y, 0.9, 1.1)
    @test isbitstype(typeof(flashed))
    @test phase_data(flashed, Val(:liquid)).mole_fractions === x
    @test phase_data(flashed, Val(:vapor)).mole_fractions === y
end

using ForwardDiff
@testset "Static flashed mixture promotion" begin
    x = @SVector [0.8, 0.2]
    y = @SVector [0.1, 0.9]
    dZ = ForwardDiff.derivative(0.4) do V
        flashed = FlashedMixture2Phase(
            MultiComponentFlash.two_phase_lv, SVector(0.125, 4.5),
            V, x, y, V + 0.5, V + 0.7)
        @test eltype(flashed.liquid.mole_fractions) === typeof(V)
        @test eltype(flashed.vapor.mole_fractions) === typeof(V)
        @test isbitstype(typeof(flashed))
        return flashed.liquid.Z
    end
    @test dZ == 1.0
end

@testset "Rachford-Rice derivatives" begin
    N = 25
    for z_light in range(0.0, 1.0, length = N)
        z = [z_light, 1.0 - z_light]
        for K1 in range(0.001, 1000.0, length = N)
            for K2 in range(0.001, 1000.0, length = N)
                K = [K1, K2]
                for V in range(0, 1, length = N)
                    f_v(V) = MultiComponentFlash.objectiveRR(V, K, z)
                    dv_ad = ForwardDiff.derivative(f_v, V)
                    dv_a = MultiComponentFlash.objectiveRR_dV(V, K, z)
                    @test dv_a ≈ dv_ad
                    f_K(K) = MultiComponentFlash.objectiveRR(V, K, z)
                    dK_ad = ForwardDiff.gradient(f_K, K)
                    for i in eachindex(K)
                        dK_a_i = MultiComponentFlash.objectiveRR_dK(V, K, z, i)
                        @test dK_a_i ≈ dK_ad[i]
                    end
                    f_z(z) = MultiComponentFlash.objectiveRR(V, K, z)
                    dz_ad = ForwardDiff.gradient(f_z, z)
                    for i in eachindex(z)
                        dz_a_i = MultiComponentFlash.objectiveRR_dz(V, K, z, i)
                        @test dz_a_i ≈ dz_ad[i]
                    end
                end
            end
        end
    end
end

@testset "Critical point measure" begin
    decane = MolecularProperty("n-Decane")
    # The light component is given with explicit properties
    mw = 0.0160428  # Molar mass (kg/mole)
    P_c = 4.5992e6  # Critical pressure (Pa)
    T_c = 190.564   # Critical temperature (°K)
    V_c = 9.4118e-5 # Critical volume (m^3/mole)
    ω = 0.22394     # Acentric factor
    methane = MolecularProperty(mw, P_c, T_c, V_c, ω)
    mixture = MultiComponentMixture((methane, decane))
    equation_of_state = GenericCubicEOS(mixture, PengRobinson())
    z = [0.4, 0.6]
    @test MultiComponentFlash.michelsen_critical_point_measure(equation_of_state, 5e6, 303.15, z) ≈ 0.776435 atol = 1e-4
    @test MultiComponentFlash.michelsen_critical_point_measure(equation_of_state, 5e6, 303.15, z, static_size = false) ≈ 0.776435 atol = 1e-4
    equation_of_state_static = make_eos_immutable(equation_of_state)
    @test MultiComponentFlash.michelsen_critical_point_measure(
        equation_of_state_static, 5e6, 303.15, SVector{2}(z)) ≈ 0.776435 atol = 1e-4

    ethane = MolecularProperty(0.03007, 4.872e6, 305.32, 1.455e-4, 0.099)
    carbon_dioxide = MolecularProperty(0.0440, 7.38e6, 304.1, 9.412e-5, 0.224)
    mixture_4 = MultiComponentMixture((methane, ethane, carbon_dioxide, decane))
    eos_4 = GenericCubicEOS(mixture_4, PengRobinson())
    eos_4_static = make_eos_immutable(eos_4)
    z_4 = @SVector [1e-16, 0.2, 0.3, 0.5 - 1e-16]
    dynamic_distance = MultiComponentFlash.michelsen_critical_point_measure(
        eos_4, 5e7, 700.0, collect(z_4); static_size = false)
    static_distance = MultiComponentFlash.michelsen_critical_point_measure(
        eos_4_static, 5e7, 700.0, z_4)
    @test static_distance ≈ dynamic_distance rtol = 1e-12
end

using StaticArrays, KernelAbstractions, JLArrays
@testset "Static flash with KernelAbstractions/JLArrays" begin
    @kernel function static_flash_kernel!(out, pressure, temperature, z, eos, storage)
        i = @index(Global)
        if i <= length(out)
            @inbounds cond = (p = pressure[i], T = temperature[i], z = z)
            K = initial_guess_K(eos, cond, storage)
            V = flash_2ph!(storage, K, eos, cond, NaN;
                method = SSIFlash(), check = false, verbose = false, z_min = nothing)
            @inbounds out[i] = V
        end
    end
    @kernel function static_bypass_kernel!(distance, bypassed, eos, z)
        i = @index(Global)
        if i <= length(distance)
            initial = (p = 5e7, T = 800.0, z = z)
            stability = stability_2ph_immutable(eos, initial;
                z_min = nothing)
            nearby = (p = 5.001e7, T = 800.01, z = z)
            next_stability = stability_2ph_immutable(eos, nearby,
                stability; z_min = nothing)
            @inbounds distance[i] = stability.storage.critical_distance
            @inbounds bypassed[i] = next_stability.bypassed
        end
    end
    if isdefined(JLArrays, :JLBackend)
        host_eos = get_test_eos()
        eos = make_eos_immutable(host_eos)
        z = @SVector [0.5, 0.3, 0.2]
        storage = flash_storage(eos, (p = 1e5, T = 300.0, z = z);
            method = SSIFlash(), static = true)
        n = 16
        pressure_host = collect(range(1e5, 4e6, length = n))
        temperature_host = collect(range(280.0, 320.0, length = n))
        expected = map(pressure_host, temperature_host) do p, T
            flash_2ph(host_eos, (p = p, T = T, z = collect(z));
                method = SSIFlash(), check = false)
        end

        pressure = JLArray(pressure_host)
        temperature = JLArray(temperature_host)
        out = JLArray(zeros(n))
        backend = JLArrays.JLBackend()
        kernel! = static_flash_kernel!(backend, 8)
        kernel!(out, pressure, temperature, z, eos, storage; ndrange = n)

        @test Array(out) ≈ expected rtol = 1e-11

        trace_z = @SVector [1e-16, 0.3, 0.7 - 1e-16]
        distance = JLArray(zeros(1))
        bypassed = JLArray(falses(1))
        bypass_kernel! = static_bypass_kernel!(backend, 1)
        bypass_kernel!(distance, bypassed, eos, trace_z; ndrange = 1)
        @test isfinite(only(Array(distance)))
        @test only(Array(bypassed))
    else
        # JLArrays 0.1 supports Julia 1.6 but predates the KernelAbstractions backend.
        @test_skip false
    end
end

