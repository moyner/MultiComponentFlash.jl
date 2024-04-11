# Basic usage

This section covers the basics of how to perform vapor-liquid flashes. For examples of how to get better performance, please see [Advanced usage](@ref).

## Vapor-liquid equilibrium for constants

Many vapor-liquid problems can be solved under the assumption of equilibrium constants (K-values). If the K-values are independent of the phase mole fractions, the vapor fraction ``V`` can be determined by a solution of the Rachford-Rice equations. We define the standard relations for molar balance:
`` z_i = V y_i + (1-V) x_i ``
where `` z_i `` is the overall mole fraction of component ``i``, ``y_i`` the vapor mole fraction of that same component and ``x_i`` the liquid fraction:
``x_i \frac{z_i}{(1 - V + V*K_i)}, \quad y_i = K_i x_i``.

The Rachford-Rice reformulation of these equations is the natural choice for a numerical solution. For more details, the [Wikipedia page on flash evaporation](https://en.wikipedia.org/wiki/Flash_evaporation) is a good starting point.

We can demonstrate this by a binary system where the first component is light and easy to vaporize (it is found in the vapor phase 9 times out of 10) and the second is heavy (being found in the liquid phase 9 times out of 10). Let us define this mixturem, and take a 70-30 mixture in moles and perform a flash to find the vapor fraction ``V``:

```jldoctest
using MultiComponentFlash
K = [0.1, 9.0] # K-values
z = [0.7, 0.3] # Mole fractions
solve_rachford_rice(K, z)

# output

0.24583333333333332
```

The result indicates that we can expect to have about 3 moles of liquid per mole of vapor.

## Two-phase multicomponent flash

For more complex mixtures, the assumption of constant K-values is not very accurate. We need to perform a full flash by defining a mixture together with an equation-of-state:

```jldoctest
using MultiComponentFlash
# Define two species: One heavy and one light.
# The heavy component uses table lookup:
decane = MolecularProperty("n-Decane")
# The light component is given with explicit properties
mw = 0.0160428  # Molar mass (kg/mole)
P_c = 4.5992e6  # Critical pressure (Pa)
T_c = 190.564   # Critical temperature (°K)
V_c = 9.4118e-5 # Critical volume (m^3/mole)
ω = 0.22394     # Acentric factor
methane = MolecularProperty(mw, P_c, T_c, V_c, ω)
# or, equivialently,
# methane = MolecularProperty("Methane")
# Create a mixture
mixture = MultiComponentMixture((methane, decane))
eos = GenericCubicEOS(mixture, PengRobinson())
# Define conditions to flash at
p = 5e6        # 5 000 000 Pa, or 50 bar
T = 303.15     # 30 °C = 303.15 °K
z = [0.4, 0.6] # 1 mole methane per 9 moles of decane
conditions = (p = p, T = T, z = z)
# Perform a flash to get the vapor fraction
V = flash_2ph(eos, conditions)
round(V, digits = 5)

# output

0.2661
```

## K-values and fractions

```@meta
DocTestSetup = quote
    using MultiComponentFlash
    decane = MolecularProperty("n-Decane")
    methane = MolecularProperty("Methane")
    mixture = MultiComponentMixture((methane, decane))
    eos = GenericCubicEOS(mixture, PengRobinson())
    m = SSIFlash()
    # Define conditions to flash at
    p = 5e6        # 5 000 000 Pa, or 50 bar
    T = 303.15     # 30 °C = 303.15 °K
    z = [0.4, 0.6] # 1 mole methane per 9 moles of decane
    conditions = (p = p, T = T, z = z)
    # Perform a flash to get the vapor fraction
    V, K, report = flash_2ph(eos, conditions, extra_out = true, method = m)
end
```

We can also get more output by turning on the `extra_out` flag. We can use this to examine the K-values (ratio between vapor and liquid mole fractions for each component):

```jldoctest
V, K, report = flash_2ph(eos, conditions, extra_out = true);
round.(K, digits = 5)

# output

2-element Vector{Float64}:
 4.12272
 0.00047
```

From the chosen overall mole fractions `z`, and the flashed `K`-values together with the vapor fraction `V` we can get the phase mole fractions in the liquid phase:

```jldoctest
julia> liquid_mole_fraction.(z, K, V)
2-element Vector{Float64}:
 0.24247146623483776
 0.7575285337651623
```

As expected, the liquid phase has more of the heavy component than in the overall mole fractions (0.75 relative to 0.6). If we compute the vapor fractions,

```jldoctest
julia> vapor_mole_fraction.(z, K, V)
2-element Vector{Float64}:
 0.9996430842149211
 0.000356915785078925
```

we see that the vapor phase is almost entirely made up of the lighter methane at the chosen conditions.

## Changing the solution algorithm to Newton

If we examine the third output, we can see output about number of iterations, the result of the phase stability test and a verification that the flash converged within the default tolerance:

```jldoctest
julia> report
(its = 8, converged = true, stability = StabilityReport (two-phase, liquid = stable, vapor = unstable))
```

Here, we can see that forming a second vapor-like phase led to the solver to conclude that the mixture is two-phase.

The default flash method is the `SSIFlash()` method. The name is an abbreviation for successive-substitution, a simple but very robust method.

We could alternatively switch to `NewtonFlash()` to use Newton's method with AD instead to reduce the number of iterations:

```jldoctest
julia> V, K, report = flash_2ph(eos, conditions, extra_out = true, method = NewtonFlash()); report
(its = 5, converged = true, stability = StabilityReport (two-phase, liquid = stable, vapor = unstable))
```

Note that Newton's method is not unconditionally stable: It is therefore also possible to use `SSINewtonFlash()` that switches from SSI to Newton after a prescribed number of iterations, which is effective around the critical region where SSI has slow convergence.
