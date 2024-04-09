"Abstract type for all flash types"
abstract type AbstractFlash end
"Abstract type for all flash types that use Newton in some form"
abstract type AbstractNewtonFlash <: AbstractFlash end

"""
Flash method that uses successive subtition.

Unconditionally convergent, does not require derivatives, but is very slow around critical regions.

See also: [`flash_2ph!`](@ref), [`NewtonFlash`](@ref) [`SSINewtonFlash`](@ref)
"""
struct SSIFlash <: AbstractFlash end

"""
    NewtonFlash([dMax = 0.2])

Flash using Newton's method for zero solve.

Only conditionally convergent, but has better convergence rate than SSI.

# Arguments
- `dMax`: dampening factor for the newton iteration

See also: [`flash_2ph!`](@ref), [`SSIFlash`](@ref) [`SSINewtonFlash`](@ref)
"""
@Base.kwdef struct NewtonFlash <: AbstractNewtonFlash
    dMax::Float64 = 0.2
end

"""
    SSINewtonFlash([swap_iter = 5,dMax = 0.2])

Perform a number of SSI iterations, followed by Newton until convergence.

See also: [`flash_2ph!`](@ref), [`SSIFlash`](@ref) [`NewtonFlash`](@ref)
"""
@Base.kwdef struct SSINewtonFlash <: AbstractNewtonFlash
    swap_iter::Int = 5
    dMax::Float64 = 0.2
end


struct PhaseStabilityStatus
    stable::Bool
    trivial::Bool
    function PhaseStabilityStatus(stable = false; trivial = stable)
        return new(stable, trivial && stable)
    end
end
struct StabilityReport
    stable::Bool
    liquid::PhaseStabilityStatus
    vapor::PhaseStabilityStatus
    function StabilityReport(;
            stable_liquid::Bool = false,
            trivial_liquid::Bool = stable_liquid,
            stable_vapor::Bool = false,
            trivial_vapor::Bool = stable_vapor,
        )
        new(stable_liquid && stable_vapor,
            PhaseStabilityStatus(stable_liquid, trivial = trivial_liquid),
            PhaseStabilityStatus(stable_vapor, trivial = trivial_vapor)
        )
    end
end
