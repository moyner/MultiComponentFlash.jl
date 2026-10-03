using ForwardDiff

primal(x) = x
primal(x::ForwardDiff.Dual) = ForwardDiff.value(x)

function central_difference(f, u, i, h)
    plus = copy(u)
    minus = copy(u)
    plus[i] += h
    minus[i] -= h
    return (f(plus) - f(minus))/(2h)
end

@testset "Implicit static flash derivatives" begin
    cubic = make_eos_immutable(get_test_eos())
    function cubic_result(u)
        z = SVector(u[3], u[4], 1-u[3]-u[4])
        ad_cond = (p = 1e6*u[1], T = u[2], z = z)
        numeric_cond = (p = primal(ad_cond.p), T = primal(ad_cond.T),
            z = map(primal, z))
        V, K = flash_2ph_immutable(cubic, numeric_cond; check = false)
        return implicit_flash_derivatives(cubic, numeric_cond, ad_cond, V, K)
    end
    u = [1.0, 300.0, 0.6, 0.1]
    @test 0 < cubic_result(u)[1] < 1
    for f in (v -> cubic_result(v)[1],
            v -> cubic_result(v)[2][1])
        for (i, h) in ((2, 1e-2), (3, 1e-5))
            derivative = ForwardDiff.derivative(t -> begin
                v = [j == i ? t : u[j] for j in eachindex(u)]
                f(v)
            end, u[i])
            @test derivative ≈ central_difference(f, u, i, h) rtol=1e-4
        end
    end

    mixture = MultiComponentMixture(["CarbonDioxide", "Water"])
    kvalue = KValuesEOS(cond -> SVector(
        0.05*(cond.p/1e6)*(cond.T/300)^0.2,
        5.0*(1e6/cond.p)*(cond.T/300)^(-0.1)), mixture)
    function kvalue_result(u)
        z = SVector(u[3], 1-u[3])
        ad_cond = (p = 1e6*u[1], T = u[2], z = z)
        K_ad = SVector(initial_guess_K(kvalue, ad_cond))
        K_numeric = map(primal, K_ad)
        numeric_cond = (p = primal(ad_cond.p), T = primal(ad_cond.T),
            z = map(primal, z))
        V = solve_rachford_rice(K_numeric, numeric_cond.z)
        return implicit_flash_derivatives(kvalue, numeric_cond, ad_cond,
            V, K_numeric, K_ad)
    end
    u_k = [1.0, 300.0, 0.4]
    @test 0 < kvalue_result(u_k)[1] < 1
    gradient = ForwardDiff.gradient(v -> kvalue_result(v)[1], u_k)
    for i in eachindex(u_k)
        f = v -> kvalue_result(v)[1]
        @test gradient[i] ≈ central_difference(f, u_k, i, 1e-4) rtol=1e-4
    end

    # A constant K evaluator can leave K_ad numeric even when the condition
    # carries AD seeds. The returned vapor fraction and K-values must still
    # use the condition's promoted type.
    constant_K = @SVector [0.05, 5.0]
    numeric_cond = (p = 1e6, T = 300.0, z = @SVector [0.4, 0.6])
    V_numeric = solve_rachford_rice(constant_K, numeric_cond.z)
    function constant_k_result(; p = numeric_cond.p, T = numeric_cond.T,
            z = numeric_cond.z)
        ad_cond = (; p, T, z)
        return implicit_flash_derivatives(kvalue, numeric_cond, ad_cond,
            V_numeric, constant_K, constant_K)
    end
    for variable in (:p, :T)
        value = getproperty(numeric_cond, variable)
        sensitivity = ForwardDiff.derivative(value) do seeded
            cond = merge(numeric_cond, NamedTuple{(variable,)}((seeded,)))
            constant_k_result(; cond...)[1]
        end
        @test sensitivity == 0
    end
    composition_sensitivity = ForwardDiff.derivative(0.4) do z1
        constant_k_result(z = SVector(z1, 1-z1))[1]
    end
    composition_finite_difference = (solve_rachford_rice(constant_K,
        SVector(0.4+1e-5, 0.6-1e-5)) - solve_rachford_rice(constant_K,
        SVector(0.4-1e-5, 0.6+1e-5)))/(2e-5)
    @test composition_sensitivity ≈ composition_finite_difference rtol=1e-6
end
