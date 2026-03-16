using Test
using MultiComponentFlash
import MultiComponentFlash.PVTExperiments as PVTExp

@testset "PVTExperiments" begin
    # Set up a simple oil mixture (methane + n-decane)
    mixture = MultiComponentMixture(["Methane", "n-Decane"])
    eos = GenericCubicEOS(mixture, PengRobinson())
    z_oil = [0.5, 0.5]
    T_res = 373.15  # 100°C

    # Gas condensate mixture (methane-heavy)
    z_gas = [0.9, 0.1]

    @testset "Saturation Pressure" begin
        p_sat, is_bp = PVTExp.find_saturation_pressure(eos, z_oil, T_res)
        @test p_sat > 0.0
        @test is_bp  # Should be bubble point for oil mixture

        p_sat_gas, is_bp_gas = PVTExp.find_saturation_pressure(eos, z_gas, T_res)
        @test p_sat_gas > 0.0
    end

    @testset "Flash and Properties" begin
        props = PVTExp.flash_and_properties(eos, 5e6, T_res, z_oil)
        @test haskey(props, :V)
        @test haskey(props, :ρ_l)
        @test haskey(props, :ρ_v)
        @test haskey(props, :μ_l)
        @test haskey(props, :μ_v)
        @test props.ρ_l > 0
        @test props.ρ_v > 0
        @test props.μ_l > 0
        @test props.μ_v > 0
    end

    @testset "CCE" begin
        result = PVTExp.cce(eos, z_oil, T_res; n_points = 10)
        @test result isa PVTExp.CCEResult
        @test result.T == T_res
        @test length(result.pressures) == 10
        @test all(result.density_oil .> 0)
        @test all(result.density_gas .> 0)
        @test all(result.viscosity_oil .> 0)
        @test result.p_sat > 0
        # Relative volume should be 1.0 near saturation pressure
        idx_sat = argmin(abs.(result.pressures .- result.p_sat))
        @test isapprox(result.relative_volume[idx_sat], 1.0, atol = 0.1)
        # Test printing
        io = IOBuffer()
        show(io, result)
        output = String(take!(io))
        @test contains(output, "CCE")
    end

    @testset "DLE" begin
        result = PVTExp.dle(eos, z_oil, T_res; n_points = 10)
        @test result isa PVTExp.DLEResult
        @test result.T == T_res
        @test length(result.pressures) == 11  # n_points + 1
        @test result.p_sat > 0
        # Bo at bubble point should be > 1
        @test result.Bo[1] >= 1.0
        # Rs should decrease with pressure
        @test result.Rs[1] >= result.Rs[end]
        # All densities should be positive
        @test all(result.density_oil .> 0)
        # Test printing
        io = IOBuffer()
        show(io, result)
        output = String(take!(io))
        @test contains(output, "DLE")
    end

    @testset "CVD" begin
        result = PVTExp.cvd(eos, z_gas, T_res; n_points = 10)
        @test result isa PVTExp.CVDResult
        @test result.T == T_res
        @test length(result.pressures) == 11  # n_points + 1
        @test result.p_sat > 0
        # Liquid dropout starts at 0 (dew point)
        @test result.liquid_dropout[1] ≈ 0.0
        # Gas Z-factors should be positive
        @test all(result.Z_gas .> 0)
        # Test printing
        io = IOBuffer()
        show(io, result)
        output = String(take!(io))
        @test contains(output, "CVD")
    end

    @testset "MSS" begin
        stages = [
            PVTExp.SeparatorStage(3e6, 330.0),   # High pressure separator
            PVTExp.SeparatorStage(101325.0, 288.706)  # Stock tank
        ]
        result = PVTExp.mss(eos, z_oil, T_res, stages)
        @test result isa PVTExp.MSSResult
        @test result.Bo > 0
        @test result.Rs >= 0
        @test result.density_oil_st > 0
        @test result.density_gas_st > 0
        @test length(result.gas_composition) == 2
        @test length(result.oil_composition) == 2
        # Test printing
        io = IOBuffer()
        show(io, result)
        output = String(take!(io))
        @test contains(output, "Separator")
    end

    @testset "PVTO Table" begin
        table = PVTExp.pvto_table(eos, z_oil, T_res; n_rs = 5, n_undersaturated = 2)
        @test table isa PVTExp.PVTOTable
        @test length(table.Rs) == length(table.p_bub)
        @test length(table.Rs) > 0
        # Rs should be non-negative and sorted
        @test all(table.Rs .>= 0)
        @test issorted(table.Rs)
        # Bo should be positive
        for bo_vec in table.Bo
            @test all(bo_vec .> 0)
        end
        # Viscosity should be positive
        for mu_vec in table.mu_o
            @test all(mu_vec .> 0)
        end
        # Test printing
        io = IOBuffer()
        show(io, table)
        output = String(take!(io))
        @test contains(output, "PVTO")
    end

    @testset "PVDG Table" begin
        table = PVTExp.pvdg_table(eos, z_gas, T_res; n_points = 10)
        @test table isa PVTExp.PVDGTable
        @test length(table.p) == 10
        @test all(table.Bg .> 0)
        @test all(table.mu_g .> 0)
        # Bg should decrease with pressure
        @test table.Bg[1] > table.Bg[end]
        # Test printing
        io = IOBuffer()
        show(io, table)
        output = String(take!(io))
        @test contains(output, "PVDG")
    end

    @testset "PVTG Table" begin
        table = PVTExp.pvtg_table(eos, z_gas, T_res; n_rv = 5, n_undersaturated = 2)
        @test table isa PVTExp.PVTGTable
        @test length(table.p) > 0
        # Bg should be positive
        for bg_vec in table.Bg
            @test all(bg_vec .> 0)
        end
        # Test printing
        io = IOBuffer()
        show(io, table)
        output = String(take!(io))
        @test contains(output, "PVTG")
    end

    @testset "Surface Densities" begin
        sd = PVTExp.surface_densities(eos, z_oil, T_res)
        @test sd isa PVTExp.SurfaceDensities
        @test sd.oil > 0
        @test sd.gas > 0
        # Oil should be denser than gas at surface
        @test sd.oil > sd.gas
        # Test printing
        io = IOBuffer()
        show(io, sd)
        output = String(take!(io))
        @test contains(output, "Surface")
    end

    @testset "High-Level Interface - Oil" begin
        tables = generate_pvt_tables(eos, z_oil, T_res;
            n_pvto = 5, n_pvdg = 10, n_undersaturated = 2)
        @test tables.pvto !== nothing
        @test tables.pvdg !== nothing
        @test tables.surface_densities isa PVTExp.SurfaceDensities
        @test tables.saturation_pressure > 0
        @test tables.is_bubblepoint == true
    end

    @testset "High-Level Interface - Gas" begin
        tables = generate_pvt_tables(eos, z_gas, T_res;
            n_pvtg = 5, n_pvdg = 10)
        @test tables.pvdg !== nothing
        @test tables.surface_densities isa PVTExp.SurfaceDensities
        @test tables.saturation_pressure > 0
    end

    @testset "High-Level Interface with Separator" begin
        stages = [
            PVTExp.SeparatorStage(3e6, 330.0),
            PVTExp.SeparatorStage(101325.0, 288.706)
        ]
        tables = generate_pvt_tables(eos, z_oil, T_res;
            separator_stages = stages,
            n_pvto = 5, n_pvdg = 10, n_undersaturated = 2)
        @test tables.pvto !== nothing
        @test tables.pvdg !== nothing
        @test tables.surface_densities isa PVTExp.SurfaceDensities
    end
end
