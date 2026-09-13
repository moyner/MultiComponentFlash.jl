module PVTExperiments
    # Notes: This module was in part generated with the help of Claude Opus. We
    # export a single function, `generate_pvt_tables`, which is the main
    # interface for generating all the tables at once. Additional experiments
    # are available, but their interfaces are not exported and can be changed
    # without breaking versions.
    using MultiComponentFlash
    using Printf

    import MultiComponentFlash: flashed_mixture_2ph, number_of_components

    const WATER_DENSITY_SC = 999.0  # kg/m³ at standard conditions
    const R_GAS = MultiComponentFlash.IDEAL_GAS_CONSTANT  # 8.3144598 J/(mol·K)

    include("types.jl")
    include("utils.jl")
    include("cce.jl")
    include("dle.jl")
    include("cvd.jl")
    include("mss.jl")
    include("tables.jl")
    include("interface.jl")

    export generate_pvt_tables
end
