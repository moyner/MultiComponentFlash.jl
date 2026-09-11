include("test_setup.jl")

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
    for m in flash_methods
        name = typeof(m)
        @testset "$name - Arrays" begin
            test_flash_inplace(m, static_size = true)
        end
        @testset "$name - StaticArrays" begin
            test_flash_inplace(m, static_size = true)
        end
    end
end

@testset "Static accelerator path" begin
    eos = static_eos(get_test_eos())
    c = (p = 1e6, T = 300.0, z = @SVector [0.5, 0.3, 0.2])
    K = zero(MVector{3, Float64})
    initial_guess_K!(K, eos, c)
    V, K, iterations, converged = flash_2ph_static(eos, c, K, 0.5)

    @test isbitstype(typeof(eos))
    @test converged
    @test iterations == 6
    @test V ≈ 0.7632068334421974
    @test K ≈ @SVector [4.553402802323027, 17.73895830809456, 0.0004031451448211194]

    K_static = initial_guess_K_static(eos, c)
    V_static, K_static, iterations_static, converged_static = flash_2ph_static(eos, c, K_static, 0.5)
    @test converged_static
    @test iterations_static == iterations
    @test V_static ≈ V
    @test K_static ≈ K

    config = FlashConfig(print_output=false, use_dict_storage=false)
    storage = flash_storage(eos, c, SSIFlash(), config)
    K_config = initial_guess_K_static(eos, c)
    V_config, K_config, report = flash_2ph!(storage, K_config, eos, c, NaN, config;
        method=SSIFlash(), extra_out=true, z_min=nothing)
    @test report.stability isa MultiComponentFlash.StabilityReport
    @test report.converged
    @test V_config ≈ V
    @test K_config ≈ K
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
end

using ForwardDiff
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
end

using StaticArrays, KernelAbstractions, JLArrays
@testset "Static flash with KernelAbstractions/JLArrays" begin
    @kernel function static_flash_kernel!(out, pressure, temperature, z, eos)
        i = @index(Global)
        if i <= length(out)
            @inbounds cond = (p = pressure[i], T = temperature[i], z = z)
            config = FlashConfig{false, false}()
            K = initial_guess_K_static(eos, cond)
            V = flash_2ph(eos, cond, K, NaN, config;
                method = SSIFlash(), check = false, verbose = false, z_min = nothing)
            @inbounds out[i] = V
        end
    end
    if isdefined(JLArrays, :JLBackend)
        host_eos = get_test_eos()
        eos = static_eos(host_eos)
        z = @SVector [0.5, 0.3, 0.2]
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
        kernel!(out, pressure, temperature, z, eos; ndrange = n)

        @test Array(out) ≈ expected rtol = 1e-11
    else
        # JLArrays 0.1 supports Julia 1.6 but predates the KernelAbstractions backend.
        @test_skip false
    end
end

