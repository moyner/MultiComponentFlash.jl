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
@allocated V = flash_2ph!(S, K, eos, conditions, method = m)

# output

16
```

See the unit tests for examples where the flash can use `StaticArrays` to avoid allocations entirely.

## Performance example

The default interface is designed for ease-of-use with standard Julia types, but the module  also supports further by using `StaticArrays`:

```julia
using MultiComponentFlash, BenchmarkTools, StaticArrays
function bench(m, static_size = false)
    p = 6e6
    T = 480.0
    # Take the SPE5 benchmark
    eos, data = cubic_benchmark("spe5")
    n = number_of_components(eos)
    z = repeat([1/n], n)
    conditions = (p = p, T = T, z = z)
    S = flash_storage(eos, conditions, method = m, static_size = static_size)
    K = initial_guess_K(eos, conditions)
    if static_size
        N = number_of_components(eos)
        K = MVector{N}(K)
    end
    V, K, status = flash_2ph!(S, K, eos, conditions, NaN, method = m, extra_out = true)
    println("V = $V (Completed in $(status.its) iterations)")
    @btime flash_2ph!($S, $K, $eos, $conditions, NaN, method = $m)
    return nothing
end
println("SSI:")
bench(SSIFlash())
println("SSI (static arrays):")
bench(SSIFlash(), true)
##
println("Newton:")
bench(NewtonFlash())
println("Newton (static arrays):")
bench(NewtonFlash(), true)
```

The output will be a bit different on other CPUs, but this flash generally takes around 20 microseconds to complete, including both stability test and flash.

```julia
SSI:
V = 0.03279769425318795 (Completed in 14 iterations)
  18.500 μs (0 allocations: 0 bytes)
SSI (static arrays):
V = 0.03279769425318795 (Completed in 14 iterations)
  16.500 μs (0 allocations: 0 bytes)

Newton:
V = 0.032797694260046494 (Completed in 4 iterations)
  20.100 μs (0 allocations: 0 bytes)
Newton (static arrays):
V = 0.032797694260046494 (Completed in 4 iterations)
  19.900 μs (0 allocations: 0 bytes)
```

!!! note "Use of `StaticArrays`"
    Switching to statically sized arrays can improve the speed, at the cost of longer compilation times. Please note that for `StaticArrays` there will be compilation that is dependent on the number of components in your mixture. For example, switching from a five to six component mixture will trigger a full recompilation of your chosen flash.

## Generate and plot a phase diagram

We create a three-component mixture and flash for a range of pressure and temperature conditions:

```julia
using MultiComponentFlash, Plots
ns = 1000
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
data = zeros(ns, ns)
for ip = 1:ns
    for iT = 1:ns
        c = (p = p[ip], T = T[iT], z = z)
        data[ip, iT] = flash_2ph!(S, K, eos, c, NaN, method = m)
    end
end

contour(p./ubar, T .- 273.15, data, levels = 10, fill=(true,cgrad(:hot)))
ylabel!("Pressure [Bar]")
xlabel!("T [°Celsius]")
```

![Phase diagram](../assets/phase_diagram_simple.png)
