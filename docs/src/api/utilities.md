# Utilities
```@index
Pages = ["utilities.md"]
```

## Partial derivatives
!!! note "Experimental features" 
    Functions for obtaining partial derivatives of the flashed results. Please note that this is an experimental feature. Examples of usage are found in the unit tests.

```@autodocs
Modules = [MultiComponentFlash]
Pages   = ["derivatives.jl"]
Order   = [:type, :function]
Private = false
```

## Coupling utilities
!!! note "Experimental features" 
    Utilities for coupling flash to simulation codes. Please note that this is an experimental feature.

```@autodocs
Modules = [MultiComponentFlash]
Pages   = ["flow_coupler.jl"]
Order   = [:type, :function]
Private = false
```


## Various
```@autodocs
Modules = [MultiComponentFlash]
Pages   = ["utils.jl"]
Order   = [:type, :function]
Private = false
```
