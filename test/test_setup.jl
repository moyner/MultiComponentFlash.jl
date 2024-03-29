using Test, MultiComponentFlash, StaticArrays

function get_test_eos(type = PengRobinson())
    co2 = MolecularProperty(0.0440, 7.38e6, 304.1, 9.412e-5, 0.224)
    c1 = MolecularProperty(0.0160, 4.60e6, 190.6, 9.863e-5, 0.011)
    c10 = MolecularProperty(0.0142, 2.10e6, 617.7, 6.098e-4, 0.488)

    mixture = MultiComponentMixture([co2, c1, c10], names = ["CO2", "C1", "C10"])
    return GenericCubicEOS(mixture, type)
end

test_conditions() = (p = 10e5, T = 300.0, z = [0.5, 0.3, 0.2])

test_allocs(S, K, eos, c, m) = @allocated flash_2ph!(S, K, eos, c, NaN, method = m)

function test_flash_inplace(m, do_test = true; static_size = false)
    eos = get_test_eos()
    c = test_conditions()
    S = flash_storage(eos, c, method = m, static_size = static_size)
    K = initial_guess_K(eos, c)
    if static_size
        n = number_of_components(eos)
        K = MVector{n}(K)
    end
    # Just in case of compilation
    test_allocs(S, K, eos, c, m)
    # Then evaluation
    allocs = test_allocs(S, K, eos, c, m)
    if do_test
        if allocs > 0
            println("Flash resulted in allocations: $allocs")
        end
        @test allocs == 0
    end
end

function test_flash_pr(m)
    eos = get_test_eos(PengRobinson())
    c = test_conditions()
    V, K, rep = flash_2ph(eos, c, method = m, extra_out = true)
    @test V ≈ 0.7632068334421974
    @test K ≈ [4.553402802323027, 17.73895830809456, 0.0004031451448211194]
    @test rep.its == 6
    @test rep.converged
end

function test_flash_srk(m)
    eos = get_test_eos(SoaveRedlichKwong())
    c = test_conditions()
    V, K, rep = flash_2ph(eos, c, method = m, extra_out = true)
    @test V ≈ 0.76245810
    @test K ≈ [4.447493042737206, 18.274788690856127, 0.0003414157433042269]
    @test rep.its == 6
    @test rep.converged
end

function test_flash_zj(m)
    eos = get_test_eos(ZudkevitchJoffe())
    c = test_conditions()
    V, K, rep = flash_2ph(eos, c, method = m, extra_out = true)
    @test V ≈ 0.7619590164
    @test K ≈ [4.254067899595903, 13.13041424605875, 0.004978904096288983]
    @test rep.its == 6
    @test rep.converged
end

function test_flash_rk(m)
    eos = get_test_eos(RedlichKwong())
    c = test_conditions()
    V, K, rep = flash_2ph(eos, c, method = m, extra_out = true)
    @test V ≈ 0.7619590164
    @test K ≈ [4.254067899595903, 13.13041424605875, 0.004978904096288983]
    @test rep.its == 6
    @test rep.converged
end


function test_flash_partials()
    eos = get_test_eos()
    c = test_conditions()

    n = number_of_components(eos)
    # Set up typical primary variable set for a simulation (p, T, z_1, ... z_(n-1))
    np = n + 1 # p, T, z1, ..., zn-1

    get_ad_scalar(v, i) = MultiComponentFlash.get_ad(v, np, :FlashTest, i)
    p = get_ad_scalar(c.p, 1)
    T = get_ad_scalar(c.T, 2)
    ∂T = typeof(p)
    z = zeros(∂T, n)

    zt = 0.0
    for i = 1:n
        if i < n
            z[i] = get_ad_scalar(c.z[i], i + 2)
            zt += z[i]
        else
            z[i] = 1 - zt
        end
    end

    S = flash_storage(eos, c, inc_jac = true, diff_externals = true, npartials = np)
    K = initial_guess_K(eos, c)
    V, K, = flash_2ph!(S, K, eos, c, NaN, extra_out = true)
    inverse_flash_update!(S, eos, c, V)

    ∂c = (p = p, T = T, z = z)
    x = liquid_mole_fraction.(z, K, V)
    set_partials_phase_mole_fractions!(x, S, eos, ∂c, :liquid)

    y = vapor_mole_fraction.(x, K)
    set_partials_phase_mole_fractions!(y, S, eos, ∂c, :vapor)

    V = set_partials_vapor_fraction(convert(∂T, V), S, eos, ∂c)
    @testset "Liquid mole fractions" begin
        @test x[1].value ≈ 0.1346989544117768
        @test x[1].partials[1] ≈ 1.228287897818972e-7
        @test x[1].partials[2] ≈ -0.0020147184844332847
        @test x[1].partials[3] ≈ 0.1197359143763492
        @test x[1].partials[4] ≈ -0.1572114180863043

        @test x[2].value ≈ 0.02177813007039174
        @test x[2].partials[1] ≈ 2.1912542920313383e-8
        @test x[2].partials[2] ≈ -0.0001541487185624217
        @test x[2].partials[3] ≈ -0.03058244779079394
        @test x[2].partials[4] ≈ 0.040154284625321765

        @test x[3].value ≈ 0.8435229155178313
        @test x[3].partials[1] ≈ -1.447413327022105e-7
        @test x[3].partials[2] ≈  0.0021688672029956967
        @test x[3].partials[3] ≈ -0.08915346658555479
        @test x[3].partials[4] ≈ 0.1170571334609829
    end

    @testset "Vapor mole fractions" begin
        @test y[1].value ≈ 0.6133385964871853
        @test y[1].partials[1] ≈ -1.249597008210954e-8
        @test y[1].partials[2] ≈ 0.00023074314715075357
        @test y[1].partials[3] ≈ 0.5450354755216841
        @test y[1].partials[4] ≈ -0.7156232151431435

        @test y[2].value ≈ 0.38632134134487833
        @test y[2].partials[1] ≈  1.2708851844895392e-8
        @test y[2].partials[2] ≈ -0.00025251614754867813
        @test y[2].partials[3] ≈ -0.5450312504283377
        @test y[2].partials[4] ≈ 0.7156176676604193

        @test y[3].value ≈ 0.00034006216793640385
        @test y[3].partials[1] ≈ -2.1288176278579448e-10
        @test y[3].partials[2] ≈ 2.177300039790918e-5
        @test y[3].partials[3] ≈ -4.22509334621591e-6
        @test y[3].partials[4] ≈ 5.547482724505394e-6
    end

    @testset "Vapor fraction" begin
        @test V.value ≈ 0.76320683344208
        @test V.partials[1] ≈ -4.084076328578355e-8
        @test V.partials[2] ≈ 0.0006287962728977443
        @test V.partials[3] ≈ 1.160941354394336
        @test V.partials[4] ≈ 1.2188608427571124
    end
end

function test_rachford_rice()
    @test solve_rachford_rice([0.1, 10.0], [0.5, 0.5]) ≈ 0.5
    @test solve_rachford_rice([0.0001, 10.0], [0.1, 0.9]) ≈ 0.88897889788978
end
