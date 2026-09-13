# Advanced usage

These are intended as more advanced examples for users who may want to use `MultiComponentFlash` as a part of another code, or get better performance by pre-allocating buffers. Please read [Basic usage](@ref) first.

## Avoiding allocations

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

    S = flash_storage(eos, conditions, method = m)
    @allocated V = flash_2ph!(S, K, eos, conditions, method = m)
end
```

If many flashes of the same mixture are to be performed at different conditions, you may want to pre-allocate the storage buffers for the flash:

```jldoctest
m = SSIFlash()
K = zeros(number_of_components(eos))
S = flash_storage(eos, conditions, method = m)
@allocated flash_2ph!(S, K, eos, conditions, method = m)

# output

16
```

## Immutable performance and GPU use

[`flash_2ph_immutable`](@ref) is the public interface to the fully static SSI path.
It is useful when the component count is small and fixed, particularly inside CPU
or GPU kernels. Convert the EOS once with [`make_eos_immutable`](@ref), and provide the
overall composition as an `SVector`:

```julia
using BenchmarkTools, MultiComponentFlash, StaticArrays

eos_static = make_eos_immutable(eos)
conditions_static = (p = p, T = T, z = SVector{length(z)}(z))

V, K = flash_2ph_immutable(eos_static, conditions_static)
@btime flash_2ph_immutable($eos_static, $conditions_static)
```

`V` is the scalar vapor fraction and `K` is an `SVector`. All working vectors are
immutable values local to the call; the input `conditions_static.z` must also be an
`SVector`. The implementation currently supports `GenericCubicEOS` with
`SSIFlash`, and compilation is specialized on the number of components.

The two-argument form creates the static storage marker automatically. It can also
be constructed once and passed as the final positional argument:

```julia
storage = flash_storage(eos_static, conditions_static; static = true)
V, K = flash_2ph_immutable(eos_static, conditions_static, storage)
@btime flash_2ph_immutable($eos_static, $conditions_static, $storage)
```

Static storage is a zero-size immutable marker rather than a mutable work buffer,
so constructing it inline normally compiles away and does not allocate. Passing it
explicitly can still be convenient when setting up a kernel. For example, each
kernel work item can construct its conditions and call:

```julia
conditions_i = (p = pressure[i], T = temperature[i], z = z_static)
V, K = flash_2ph_immutable(eos_static, conditions_i, storage)
```

For a sequence of nearby states, the immutable Michelsen stability bypass can
reuse the last fully tested single-phase condition:

```julia
V, K, stability = flash_2ph_immutable(eos_static, conditions_static;
    return_stability = true)

next_conditions = (p = 1.001p, T = T + 0.01,
    z = SVector{length(z)}(z))
V, K, stability = flash_2ph_immutable(eos_static, next_conditions;
    stability_storage = stability,
    return_stability = true)
```

The stability calculation is also available on its own with
`stability_2ph_immutable`. Its result and nested storage are isbits values and
can be passed through accelerator kernels. `stability.bypassed` indicates
whether the full Michelsen test was skipped. Storage is only armed when both
trial phases converged to trivial stable solutions; the shadow region is
always retested.

Do not share ordinary mutable storage from `flash_storage(...; static = false)`
between kernel work items.

## Generate and plot a phase diagram

We create a three-component mixture and flash for a range of pressure and temperature conditions:

```@example phase-diagram
using MultiComponentFlash, CairoMakie
CairoMakie.activate!(type = "svg")

ns = 100
ubar = 1e5
# Pressure range
p0 = 1*ubar
p1 = 120*ubar
# Temperature range
T0 = 273.15 + 1
T1 = 263.15 + 350
# Define mixture + eos
names = ["Methane", "CarbonDioxide", "n-Decane"]
props = MolecularProperty.(names)
mixture = MultiComponentMixture(props)
eos = GenericCubicEOS(mixture)
# Constant mole fractions, vary p-T
z = [0.3, 0.1, 0.6]
p  = range(p0, p1, length = ns)
T = range(T0, T1, length = ns)
cond = (p = p0, T = T0, z = z)

m = SSIFlash()
S = flash_storage(eos, cond, method = m)
K = initial_guess_K(eos, cond)
data = zeros(length(T), length(p))
for (iT, temperature) in pairs(T)
    for (ip, pressure) in pairs(p)
        c = (p = pressure, T = temperature, z = z)
        V = flash_2ph!(S, K, eos, c, NaN, method = m)
        data[iT, ip] = V
    end
end

fig = Figure(size = (760, 480))
ax = Axis(fig[1, 1];
    xlabel = "Temperature [°C]",
    ylabel = "Pressure [bar]",
    title = "Vapor fraction")
contours = contourf!(ax, T .- 273.15, p./ubar, data;
    levels = range(0.0, 1.0, length = 11), colormap = :hot)
Colorbar(fig[1, 2], contours; label = "Vapor mole fraction")
fig
```

The colored region is the two-phase envelope. Unfilled states are stable
single-phase conditions, for which the two-phase solver reports no intermediate
vapor fraction.

### PVT table generation

There is experimental support for generating simulator input blackoil tables (e.g. PVTG/PVDG and PVTO/PVDO).

```@docs
generate_pvt_tables
```

### Coupling to simulators and other utilities

```@docs
cubic_benchmark
FlashedMixture2Phase
FlashedPhase
```
