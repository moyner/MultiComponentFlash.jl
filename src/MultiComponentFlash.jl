
module MultiComponentFlash
    using LinearAlgebra, ForwardDiff, StaticArrays, Roots
    # Constants.
    const MINIMUM_COMPOSITION = 1e-10
    const IDEAL_GAS_CONSTANT = 8.3144598
    # Load types first
    include("mixture_types.jl")
    include("eos_types.jl")
    include("flash_types.jl")
    include("flow_coupler_types.jl")
    # Types that define specific cubic equations of state
    export SoaveRedlichKwong, RedlichKwong, PengRobinson, PengRobinsonCorrected, ZudkevitchJoffe
    # The generic cubic form that supports the above
    export GenericCubicEOS
    export number_of_components
    # Flash interfaces
    export flash_2ph, flash_2ph!, flash_storage
    export stability_2ph, stability_2ph!
    # Algorithms for flash
    export SSIFlash, NewtonFlash, SSINewtonFlash
    # Mixtures and their molecular makeup
    export MolecularProperty, MultiComponentMixture
    # K-values
    export wilson_estimate, wilson_estimate!, initial_guess_K, initial_guess_K!
    # Vapor-liquid equilibrium
    export solve_rachford_rice

    export liquid_mole_fraction, vapor_mole_fraction
    export component_fugacity, mixture_fugacities, mixture_compressibility_factor

    export force_scalars, force_coefficients, force_coefficients!
    export critical_pressure, critical_temperature, critical_volume, acentric_factor, molar_weight
    export cubic_benchmark
    export single_phase_label

    # Types and utilities for coupling flash to other codes.
    export SinglePhaseLiquid, SinglePhaseVapor, TwoPhaseLiquidVapor, FlashedPhase, UnknownPhaseState, FlashedMixture2Phase
    export phase_saturations, mass_density, mass_densities, molar_volume, phase_is_present
    export inverse_flash_update!
    export lbc_viscosity, lbc_viscosities
    export set_partials, set_partials_phase_mole_fractions!, set_partials_vapor_fraction

    include("mixtures.jl")
    include("kvalues.jl")
    include("rachford_rice.jl")
    include("eos.jl")
    include("flash.jl")
    include("derivatives.jl")
    include("stability.jl")
    include("tables.jl")

    include("flow_coupler.jl")
    include("utils.jl")
end
