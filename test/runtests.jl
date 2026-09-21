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

@testset "Float32 flashed mixture storage" begin
    x = SVector(0.8f0, 0.2f0)
    y = SVector(0.1f0, 0.9f0)
    K = SVector(0.125, 4.5)
    flashed = FlashedMixture2Phase(
        MultiComponentFlash.two_phase_lv, K, 0.4f0,
        x, y, 0.9, 1.1)
    @test flashed.V isa Float32
    @test flashed.liquid.Z isa Float32
    @test eltype(flashed.liquid.mole_fractions) === Float32

    widened = FlashedMixture2Phase(
        MultiComponentFlash.two_phase_lv, K, 0.4,
        SVector(0.8, 0.2), SVector(0.1, 0.9), 0.9, 1.1)
    @test isbitstype(typeof(widened))
    target = typeof(flashed)
    converted = convert(target, widened)
    @test converted.V isa Float32
    @test converted.liquid.mole_fractions == x
    @test isequal(converted.flash_cond.z, widened.flash_cond.z)

    narrow_K = SVector(0.125f0, 4.5f0)
    narrow = FlashedMixture2Phase(
        MultiComponentFlash.two_phase_lv, narrow_K, 0.4f0,
        x, y, 0.9f0, 1.1f0, Float32(NaN),
        (p = 1.0f6, T = 300.0f0, z = narrow_K))
    @test narrow.critical_distance isa Float32
    @test isbitstype(typeof(narrow))
    @test narrow.flash_cond.p isa Float32
    @test narrow.flash_cond.T isa Float32
    @test eltype(narrow.flash_cond.z) === Float32
    converted_narrow = convert(typeof(narrow), widened)
    @test converted_narrow.K == narrow_K
    @test converted_narrow.flash_cond.p isa Float32
    partial_target = FlashedMixture2Phase{Float32, typeof(x), typeof(narrow_K)}
    @test convert(partial_target, widened) isa typeof(narrow)
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

@testset "Static flashed mixture allocation" begin
    function allocated_static_flash(V)
        x = SVector(0.8*one(V), 0.2*one(V))
        y = SVector(0.1*one(V), 0.9*one(V))
        K = SVector(0.125, 4.5)
        cond = (p = 1.0, T = 273.15, z = SVector(0.5, 0.5))
        allocated_bytes = @allocated flashed = FlashedMixture2Phase(
            MultiComponentFlash.two_phase_lv, K, V, x, y,
            one(V), one(V), NaN, cond)
        return allocated_bytes, flashed
    end

    # The flash is stored once per cell, so even a small per-result allocation
    # produces substantial GC traffic on reservoir grids.
    observed_bytes = Ref{Int}(0)
    for _ in 1:2
        ForwardDiff.derivative(0.4) do V
            observed_bytes[], flashed = allocated_static_flash(V)
            @test flashed.V == V
            return flashed.V
        end
    end
    if VERSION >= v"1.12"
        @test observed_bytes[] == 0
    end
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

