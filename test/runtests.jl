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
