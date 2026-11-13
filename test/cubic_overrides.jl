using MultiComponentFlash, StaticArrays, Test

@testset "Cubic coefficient overrides" begin
    MCF = MultiComponentFlash
    mix = MultiComponentMixture(["Water", "Methane", "Hydrogen"])
    c = (p = 9e6, T = 320.0, z = SVector(0.8, 0.18, 0.02))
    shift = (0.01, 0.02, 0.03)
    sw = SoreideWhitson(mix)
    base = GenericCubicEOS(mix, sw)
    uniform = GenericCubicEOS(mix, sw;
        omega_a = ntuple(_ -> base.ω_a, 3), omega_b = ntuple(_ -> base.ω_b, 3))
    @test typeof(uniform) === typeof(base)
    @test uniform.ω_a === base.ω_a
    @test uniform.ω_b === base.ω_b
    legacy = GenericCubicEOS{typeof(sw), Float64, 3, Nothing, typeof(mix)}(
        sw, mix, base.m_1, base.m_2, base.ω_a, base.ω_b, nothing)
    @test typeof(legacy) === typeof(base)
    scalar = GenericCubicEOS(mix, sw; omega_a = 0.46f0, omega_b = 0.08, m1 = 1, m2 = 2)
    @test (scalar.ω_a, scalar.ω_b, scalar.m_1, scalar.m_2) === (Float64(0.46f0), 0.08, 1.0, 2.0)
    @test_throws ArgumentError GenericCubicEOS(mix, sw; omega_a = (0.46, 0.45))

    for eos_type in (sw, ZudkevitchJoffe(F_a = (T, i) -> 1 + 0.01i,
                                         F_b = (T, i) -> 1 - 0.01i))
        base = GenericCubicEOS(mix, eos_type)
        a = (0.45, 0.46, 0.47)
        b = (0.076, 0.079, 0.08)
        eos = GenericCubicEOS(mix, eos_type; omega_a = a, omega_b = b, volume_shift = shift)
        @test eos.type === eos_type
        for phase in (:liquid, :vapor), i in 1:3
            cond = (; c..., phase)
            @test MCF.A_i(eos, cond, i) ≈ MCF.A_i(base, cond, i)*(a[i]/base.ω_a)
            @test MCF.B_i(eos, cond, i) ≈ MCF.B_i(base, cond, i)*(b[i]/base.ω_b)
        end
        # The volume shift uses the same component-specific covolume.
        ref_shift = sum(i -> begin
            prop = mix.properties[i]
            weight = MCF.weight_bi(base, c, i)*(b[i]/base.ω_b)
            shift[i]*c.z[i]*weight*MCF.IDEAL_GAS_CONSTANT*MCF.critical_temperature(prop)/MCF.critical_pressure(prop)
        end, 1:3)
        ph = FlashedPhase(c.z, 1.0)
        @test MCF.molar_volume(base, c.p, c.T, ph) - MCF.molar_volume(eos, c.p, c.T, ph) ≈ ref_shift
        for F in (Float32, Float64)
            static = make_eos_immutable(eos; float_type = F)
            @test isbitstype(typeof(static))
            @test static.ω_a === map(F, a)
            @test static.ω_b === map(F, b)
        end
        if eos_type isa SoreideWhitson
            V0, K0, report = flash_2ph(eos, (; c..., z = collect(c.z)); extra_out = true, z_min = nothing)
            V, K = flash_2ph_immutable(make_eos_immutable(eos), c; z_min = nothing)
            @test report.converged
            @test V ≈ V0
            @test K ≈ K0
            xw = liquid_mole_fraction(c.z[1], K[1], V)
            @test xw > 0.99 && K[1]*xw < 0.01
        end
    end
end
