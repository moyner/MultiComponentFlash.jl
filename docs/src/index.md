# MultiComponentFlash Manual

```@meta
CurrentModule = MultiComponentFlash
```

Welcome to the documentation for [MultiComponentFlash](https://github.com/moyner/MultiComponentFlash.jl).

This package implements several equations of state for multicomponent vapor-liquid equilibrium, also called flash, for mixtures. The following equations of state (EOS) are implemented as a class of generic cubics:

* Peng-Robinson
* Redlich-Kwong
* Soave-Redlich-Kwong
* Zudkevitch-Joffe

The code is fully type stable, easy to use and fairly performant, with additional options to avoid allocations if you need to perform many flashes. A small example that benchmarks the performance:
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

```
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