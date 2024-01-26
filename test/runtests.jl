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
    for z_light in range(0.0, 1.0, N)
        z = [z_light, 1.0 - z_light]
        for K1 in range(0.001, 1000.0, N)
            for K2 in range(0.001, 1000.0, N)
                K = [K1, K2]
                for V in range(0, 1, N)
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
