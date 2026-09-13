const MCF = MultiComponentFlash

@testset "Cubic EOS defining equations" begin
    property = MolecularProperty(0.050, 5.0e6, 400.0, 1.0e-4, 0.80)
    mixture = MultiComponentMixture((property,); names = ["test component"])
    cond = (p = 2.0e6, T = 320.0, z = [1.0])
    Tr = cond.T/property.T_c

    @testset "Static and temperature coefficients" begin
        pr = GenericCubicEOS(mixture, PengRobinson())
        prc = GenericCubicEOS(mixture, PengRobinsonCorrected())
        srk = GenericCubicEOS(mixture, SoaveRedlichKwong())
        rk = GenericCubicEOS(mixture, RedlichKwong())

        @test collect(MCF.static_coefficients(PengRobinson())) ≈
            [0.457235529, 0.077796074, 1 + sqrt(2), 1 - sqrt(2)]
        @test collect(MCF.static_coefficients(SoaveRedlichKwong())) ≈
            [0.4274802327, 0.08664035, 0.0, 1.0]

        kappa_pr = 0.37464 + 1.54226property.ω - 0.26992property.ω^2
        alpha_pr = (1 + kappa_pr*(1 - sqrt(Tr)))^2
        @test MCF.weight_ai(pr, cond, 1) ≈ 0.457235529*alpha_pr

        kappa_pr78 = 0.379642 + 1.48503property.ω -
            0.164423property.ω^2 + 0.016666property.ω^3
        alpha_pr78 = (1 + kappa_pr78*(1 - sqrt(Tr)))^2
        @test MCF.weight_ai(prc, cond, 1) ≈ 0.457235529*alpha_pr78

        kappa_srk = 0.48 + 1.574property.ω - 0.176property.ω^2
        alpha_srk = (1 + kappa_srk*(1 - sqrt(Tr)))^2
        @test MCF.weight_ai(srk, cond, 1) ≈ 0.4274802327*alpha_srk
        @test MCF.weight_ai(rk, cond, 1) ≈ 0.4274802327/sqrt(Tr)
    end

    @testset "Zudkevitch-Joffe modifiers" begin
        Fa(T, i) = 1 + 1.0e-3T + 0.1i
        Fb(T, i) = 0.8 + 5.0e-4T + 0.05i
        zj = GenericCubicEOS(mixture, ZudkevitchJoffe(; F_a = Fa, F_b = Fb))
        @test MCF.weight_ai(zj, cond, 1) ≈
            0.4274802327*Fa(cond.T, 1)/sqrt(Tr)
        @test MCF.weight_bi(zj, cond, 1) ≈ 0.08664035*Fb(cond.T, 1)

        rk = GenericCubicEOS(mixture, RedlichKwong())
        zj_default = GenericCubicEOS(mixture, ZudkevitchJoffe())
        @test MCF.weight_ai(zj_default, cond, 1) == MCF.weight_ai(rk, cond, 1)
        @test MCF.weight_bi(zj_default, cond, 1) == MCF.weight_bi(rk, cond, 1)
    end

    @testset "Generalized cubic polynomial" begin
        for eos_type in (PengRobinson(), SoaveRedlichKwong(), RedlichKwong())
            eos = GenericCubicEOS(mixture, eos_type)
            A = 0.31
            B = 0.047
            polynomial = MCF.cubic_polynomial(eos, A, B)
            for Z in (0.09, 0.8, 1.3)
                cubic_residual = Z^3 + polynomial[1]*Z^2 +
                    polynomial[2]*Z + polynomial[3]
                eos_residual = 1 - 1/(Z - B) +
                    A/((Z + eos.m_1*B)*(Z + eos.m_2*B))
                denominator = (Z - B)*(Z + eos.m_1*B)*(Z + eos.m_2*B)
                @test cubic_residual ≈ eos_residual*denominator
            end
        end
    end
end

@testset "Soreide-Whitson correlations" begin
    molality = 2.5
    sw = SoreideWhitson(["Water", "generic hydrocarbon"]; molality = molality)
    acentric = 0.35
    Tr = 1.20

    A0 = sw.A[1] + sw.A_mw[1]*sign(acentric)*abs(acentric)^(-0.1)
    A1 = sw.A[2] + sw.A_mw[2]*acentric
    A2 = sw.A[3] + sw.A_mw[3]*acentric
    expected_hc_bic = A0*(1 + sw.alphas[1]*molality) +
        A1*Tr*(1 + sw.alphas[2]*molality) +
        A2*Tr^2*(1 + sw.alphas[3]*molality)
    @test MCF.soreide_whitson_hc_aqueous_bic(sw, acentric, Tr) ≈
        expected_hc_bic

    expected_n2_bic = -1.70235*(1 + 0.025587molality^0.75) +
        0.44338*(1 + 0.08126molality^0.75)*Tr
    expected_co2_bic = -0.31092*(1 + 0.15587molality^0.7505) +
        0.23580*(1 + 0.17837molality^0.979)*Tr -
        21.2566exp(-6.7222Tr - molality)
    @test MCF.soreide_whitson_n2_aqueous_bic(sw, acentric, Tr) ≈ expected_n2_bic
    @test MCF.soreide_whitson_h2s_aqueous_bic(sw, acentric, Tr) ≈
        -0.20441 + 0.234267Tr
    @test MCF.soreide_whitson_co2_aqueous_bic(sw, acentric, Tr) ≈
        expected_co2_bic

    water = MolecularProperty("Water")
    hydrocarbon = MolecularProperty(0.050, 5.0e6, 400.0, 1.0e-4, acentric)
    interaction = [0.0 0.12; 0.12 0.0]
    mixture = MultiComponentMixture((water, hydrocarbon);
        names = ["Water", "generic hydrocarbon"], A_ij = interaction)
    sw = SoreideWhitson(mixture; molality = molality)
    eos = GenericCubicEOS(mixture, sw)
    cond = (p = 1.0e6, T = 350.0, z = [0.4, 0.6])
    water_Tr = cond.T/water.T_c
    alpha_half = 1 + 0.4530*(1 - (1 - 0.0103molality^1.1)*water_Tr) +
        0.0034*(water_Tr^(-3) - 1)
    @test MCF.weight_ai(eos, cond, 1) ≈ 0.457235529*alpha_half^2

    liquid = MCF.set_phase(cond, :liquid)
    vapor = MCF.set_phase(cond, :vapor)
    expected = MCF.soreide_whitson_hc_aqueous_bic(
        sw, acentric, cond.T/hydrocarbon.T_c)
    @test MCF.binary_interaction(eos, 1, 2, liquid) ≈ expected
    @test MCF.binary_interaction(eos, 2, 1, liquid) ≈ expected
    @test MCF.binary_interaction(eos, 1, 2, vapor) == 0.12

    storage = flash_storage(eos, cond)
    @test MCF.get_force_coefficients(storage.forces, eos, liquid) ===
        storage.forces.liquid
    @test MCF.get_force_coefficients(storage.forces, eos, vapor) ===
        storage.forces.vapor
    vapor_fraction, K, report = flash_2ph(eos, cond; extra_out = true)
    @test report.converged
    @test isfinite(vapor_fraction)
    @test all(isfinite, K)
end

@testset "Cubic root degeneracies" begin
    @test MCF.solve_cubic_positive_roots(-3.0, 3.0, -1.0) ≈ 1.0
    @test MCF.solve_cubic_positive_roots(0.0, 0.0, 0.0) ≈ 0.0
    @test sort(collect(MCF.solve_cubic_positive_roots(-4.0, 5.0, -2.0))) ≈
        [1.0, 1.0, 2.0]

    root = MCF.solve_cubic_positive_roots(0.0, 1.0, 1.0)
    @test root^3 + root + 1 ≈ 0.0 atol = 1.0e-14
end
