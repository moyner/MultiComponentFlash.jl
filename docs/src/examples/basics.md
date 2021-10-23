
# Introduction

## Vapor-liquid equilibrium for constants


## Two-phase multicomponent flash
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

# output

0.265849860967563
```

### K-values and fractions
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
K

# output

2-element Vector{Float64}:
 4.119472435769978
 0.00048241080032170367
```
From the chosen overall mole fractions `z`, and the flashed `K`-values together with the vapor fraction `V` we can get the phase mole fractions in the liquid phase:
```jldoctest
julia> liquid_mole_fraction.(z, K, V)
2-element Vector{Float64}:
 0.24266084237653415
 0.7573391576234659
```
As expected, the liquid phase has more of the heavy component than in the overall mole fractions (0.75 relative to 0.6). If we compute the vapor fractions,
```jldoctest
julia> vapor_mole_fraction.(z, K, V)
2-element Vector{Float64}:
 0.999634651410856
 0.00036534858914410106
```
we see that the vapor phase is almost entirely made up of the lighter methane at the chosen conditions.

### Switching algorithms
If we examine the third output, we can see output about number of iterations and a verification that the flash converged within the default tolerance:
```jldoctest
julia> report
(its = 8, converged = true)
```
The default method is the `SSIFlash()` method. The name is an abbreviation for successive-substitution, a simple but very robust method.

We could alternatively switch to `NewtonFlash()` to use Newton's method with AD instead to reduce the number of iterations:
```jldoctest
julia> V, K, report = flash_2ph(eos, conditions, extra_out = true, method = NewtonFlash()); report
(its = 5, converged = true)
```
It is also possible to use `SSINewtonFlash()` that switches from SSI to Newton after a prescribed number of iterations, which is effective around the critical region where SSI has slow convergence.

## Avoiding allocations
If many flashes of the same mixture are to be performed at different conditions, you may want to pre-allocate the storage buffers for the flash:
```jldoctest
m = SSIFlash()
K = zeros(number_of_components(eos))
S = flash_storage(eos, conditions, method = m)
@allocated V = flash_2ph!(S, K, eos, conditions, method = m)

# output

16
```
See the unit tests for examples where the flash can use `StaticArrays` to avoid allocations entirely.