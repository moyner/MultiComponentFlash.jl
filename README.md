# MultiComponentFlash

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://moyner.github.io/MultiComponentFlash.jl/dev)
[![Build Status](https://github.com/moyner/MultiComponentFlash.jl/workflows/CI/badge.svg)](https://github.com/moyner/MultiComponentFlash.jl/actions)
[![Coverage](https://codecov.io/gh/moyner/MultiComponentFlash.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/moyner/MultiComponentFlash.jl)


# Introduction
This package implements several equations of state for multicomponent vapor-liquid equilibrium, also called flash, for mixtures. These can be used to determine vapor fractions, molar partition between the phases and predict properties such as density and volume.

 The following equations of state (EOS) are implemented as a class of generic cubics:

* [Peng-Robinson](https://doi.org/10.1021/i160057a011)
* [Redlich-Kwong](https://doi.org/10.1021/cr60137a013)
* [Soave-Redlich-Kwong](https://doi.org/10.1016/0009-2509(72)80096-4)
* [Zudkevitch-Joffe](https://doi.org/10.1002/aic.690160122)

The code is fully type stable, easy to use and fairly performant, with additional options to avoid allocations if you need to perform many flashes. The main implementation goal is to have a compact, performant and easy to use code suitable for integration in simulators of multiphase flow.

# Highlights
We highlight a few of the features. For more details, please see the API and the [Basic usage](@ref) examples.

* Support for different solution methods (successive substition (SSI), Newton, SSI+Newton) and type dispatch makes it easy for users to add their own.
* User friendly interface - with options to pre-allocate for increase performance.
* Compatible with AD. Mostly tested with `ForwardDiff`.
* Phase stability test.
* Rachford-Rice for constant K-values.
* Useful utilities for molar density, molar volume, single-phase liquid-vapor estimation and Lohrenz-Bray-Clark viscosity correlations that are compatible with custom types and AD.

For more details, please see the [documentation](https://moyner.github.io/MultiComponentFlash.jl/stable)

# Quick start
## Vapor fraction for constant K-values
We solve the Rachford-Rice equations:
```julia
using MultiComponentFlash
K = [0.1, 9.0] # K-values
z = [0.7, 0.3] # Mole fractions
V = solve_rachford_rice(K, z)
```
`V` now holds the vapor fraction:
```
0.24583333333333332
```
## Flash with cubic equations of state
Solve vapor-liquid equilibrium with a two-component mixture and the Peng-Robinson equation of state:
```julia
props = MolecularProperty.(["Methane", "n-Decane"])
mixture = MultiComponentMixture(props)
eos = GenericCubicEOS(mixture, PengRobinson())
# Define conditions to flash at
p = 5e6        # 5 000 000 Pa, or 50 bar
T = 303.15     # 30 °C = 303.15 °K
z = [0.4, 0.6] # 1 mole methane per 9 moles of decane
conditions = (p = p, T = T, z = z)
# Perform a flash to get the vapor fraction
V = flash_2ph(eos, conditions)
```

```
julia> V
0.20785284212697513
```
## Get K-values and calculate liquid and vapor mole fractions
If we also want to know how the components are partitioned in the two phases, we can turn on the `extra_out` flag.
```
V, K, = flash_2ph(eos, conditions, extra_out = true)
x = liquid_mole_fraction.(z, K, V)
y = vapor_mole_fraction.(z, K, V)
```

```
julia> x
2-element Vector{Float64}:
 0.24266084237653415
 0.7573391576234659

julia> y
2-element Vector{Float64}:
 0.999634651410856
 0.00036534858914410106

julia> K
2-element Vector{Float64}:
 4.119472435769978
 0.00048241080032170367
 ```

## Generation of phase diagrams
More examples, and details are found in the [documentation](https://moyner.github.io/MultiComponentFlash.jl/stable). Here is a p-T phase diagram for methane, n-decane and carbon dioxide found in the advanced examples:
![Phase diagram](docs/src/assets/phase_diagram_simple.png)


