var documenterSearchIndex = {"docs":
[{"location":"","page":"Home","title":"Home","text":"CurrentModule = MultiComponentFlash","category":"page"},{"location":"#MultiComponentFlash","page":"Home","title":"MultiComponentFlash","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Documentation for MultiComponentFlash.","category":"page"},{"location":"","page":"Home","title":"Home","text":"","category":"page"},{"location":"","page":"Home","title":"Home","text":"Modules = [MultiComponentFlash]","category":"page"},{"location":"#MultiComponentFlash.AbstractFlash","page":"Home","title":"MultiComponentFlash.AbstractFlash","text":"Abstract type for all flash types\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.AbstractGeneralizedCubic","page":"Home","title":"MultiComponentFlash.AbstractGeneralizedCubic","text":"Abstract base type for generalized cubics.\n\nAny subtype may override the any of following internal functions to alter the EOS:\n\n* static_coefficients\n* weight_ai\n* weight_bi\n\nIn addition, the GenericCubicEOS is parametric on the specific cubic type for further dispatch modifications.\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.AbstractNewtonFlash","page":"Home","title":"MultiComponentFlash.AbstractNewtonFlash","text":"Abstract type for all flash types that use Newton in some form\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.FlashedMixture2Phase","page":"Home","title":"MultiComponentFlash.FlashedMixture2Phase","text":"Type that holds liquid and vapor phase states together with their state\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.FlashedPhase","page":"Home","title":"MultiComponentFlash.FlashedPhase","text":"Type that holds values for a flashed phase (mole fractions + compressibility factor)\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.GenericCubicEOS","page":"Home","title":"MultiComponentFlash.GenericCubicEOS","text":"GenericCubicEOS is an implementation of generalized cubic equations of state.\n\nMany popular cubic equations can be written in a single form with a few changes in definitions for the terms (they are, after all, all cubic in form). References:\n\nCubic Equations of State-Which? by J.J. Martin\nSimulation of Gas Condensate Reservoir Performance  by K.H. Coats\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.GenericCubicEOS-2","page":"Home","title":"MultiComponentFlash.GenericCubicEOS","text":"GenericCubicEOS(mixture, [type = PengRobinson()])\n\nInstantiate a generic cubic equation-of-state for a MultiComponentMixture and  a specified EOS.\n\nCurrently supported choices for type:\n\n1. `PengRobinson` (default)\n2. `ZudkevitchJoffe`\n3. `RedlichKwong`\n4. `SoaveRedlichKwong`\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.MolecularProperty","page":"Home","title":"MultiComponentFlash.MolecularProperty","text":"MolecularProperty(molar_mass, p_crit, T_crit, V_crit, acentric_factor = 0.0)\n\nType that defines the static properties of a molecular species.\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.MolecularProperty-Tuple{String}","page":"Home","title":"MultiComponentFlash.MolecularProperty","text":"MolecularProperty(\"Name\")\n\nConvenience constructor that looks up molecular properties from a table in tabulated_properties.\n\nThe properties are taken from the wonderful MIT-licensed CoolProp. Please note that the equations of state included in this module may not be approprioate for all the available fluids, especially for mixtures!\n\nSee list of species at CoolProp website.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.MultiComponentMixture","page":"Home","title":"MultiComponentFlash.MultiComponentMixture","text":"MultiComponentMixture(properties::NTuple{N, MolecularProperty}; A_ij = nothing, names = [\"C1\", \"C2\", ...], name = \"UnnamedMixture\")\n\nCreate a multicomponent mixture with an optional binary interaction coefficient matrix A_ij.\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.NewtonFlash","page":"Home","title":"MultiComponentFlash.NewtonFlash","text":"Flash using Newton's method for zero solve.\n\nOnly conditionally convergent, but has better convergence rate than SSI.\n\nSee also: flash_2ph!, SSIFlash SSINewtonFlash\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.PengRobinson","page":"Home","title":"MultiComponentFlash.PengRobinson","text":"Specializes the GenericCubicEOS to the Peng-Robinson cubic equation of state.\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.RedlichKwong","page":"Home","title":"MultiComponentFlash.RedlichKwong","text":"Specializes the GenericCubicEOS to the Redlich-Kwong cubic equation of state.\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.SSIFlash","page":"Home","title":"MultiComponentFlash.SSIFlash","text":"Flash method that uses successive subtition.\n\nUnconditionally convergent, does not require derivatives, but is very slow around critical regions.\n\nSee also: flash_2ph!, NewtonFlash SSINewtonFlash\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.SSINewtonFlash","page":"Home","title":"MultiComponentFlash.SSINewtonFlash","text":"SSINewtonFlash([swap_iter = 5])\n\nPerform a number of SSI iterations, followed by Newton until convergence.\n\nSee also: flash_2ph!, SSIFlash NewtonFlash\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.SinglePhaseLiquid","page":"Home","title":"MultiComponentFlash.SinglePhaseLiquid","text":"Single-phase liquid state for dispatch\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.SinglePhaseVapor","page":"Home","title":"MultiComponentFlash.SinglePhaseVapor","text":"Single-phase vapor state for dispatch\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.SoaveRedlichKwong","page":"Home","title":"MultiComponentFlash.SoaveRedlichKwong","text":"Specializes the GenericCubicEOS to the Soave-Redlich-Kwong cubic equation of state.\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.TwoPhaseLiquidVapor","page":"Home","title":"MultiComponentFlash.TwoPhaseLiquidVapor","text":"Two-phase liquid-vapor state for dispatch\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.UnknownPhaseState","page":"Home","title":"MultiComponentFlash.UnknownPhaseState","text":"Unknown phase state (not initialized)\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.ZudkevitchJoffe","page":"Home","title":"MultiComponentFlash.ZudkevitchJoffe","text":"Specializes the GenericCubicEOS to the Zudkevitch-Joffe cubic equation of state.\n\nThe Zudkevitch-Joffe equations of state allows for per-component functions of  temperature that modify the weight_ai and weight_bi functions. These additional fitting parameters allows for more flexibility when matching complex mixtures.\n\n\n\n\n\n","category":"type"},{"location":"#MultiComponentFlash.acentric_factor-Union{Tuple{MolecularProperty{R}}, Tuple{R}} where R","page":"Home","title":"MultiComponentFlash.acentric_factor","text":"Get the acentric factorfor a species.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.component_fugacity-Tuple{GenericCubicEOS, Any, Any, Any, Any, Any}","page":"Home","title":"MultiComponentFlash.component_fugacity","text":"component_fugacity(eos, cond, i, Z, forces, scalars)\n\nGet fugacity of component i in a phase with compressibility Z and EOS constants scalars.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.component_fugacity_coefficient-Tuple{MultiComponentFlash.AbstractCubicEOS, Any, Any, Any, Any, Any}","page":"Home","title":"MultiComponentFlash.component_fugacity_coefficient","text":"Fugacity\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.critical_pressure-Union{Tuple{MolecularProperty{R}}, Tuple{R}} where R","page":"Home","title":"MultiComponentFlash.critical_pressure","text":"Get the critical pressure for a species.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.critical_temperature-Union{Tuple{MolecularProperty{R}}, Tuple{R}} where R","page":"Home","title":"MultiComponentFlash.critical_temperature","text":"Get the critical temperature for a species.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.critical_volume-Union{Tuple{MolecularProperty{R}}, Tuple{R}} where R","page":"Home","title":"MultiComponentFlash.critical_volume","text":"Get the critical volume for a species.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.cubic_benchmark-Tuple{Any}","page":"Home","title":"MultiComponentFlash.cubic_benchmark","text":"cubic_benchmark(name)\n\nGet a benchmark equation-of-state instance together with some data for testing.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.cubic_polynomial-Tuple{GenericCubicEOS, Real, Real}","page":"Home","title":"MultiComponentFlash.cubic_polynomial","text":"(Precompute version)\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.eos_polynomial-Tuple{GenericCubicEOS, Any, Any}","page":"Home","title":"MultiComponentFlash.eos_polynomial","text":"Generic version (only depend on scalars, ignores forces)\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.eos_polynomial-Tuple{GenericCubicEOS, Any}","page":"Home","title":"MultiComponentFlash.eos_polynomial","text":"(Expensive version)\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.flash_2ph!","page":"Home","title":"MultiComponentFlash.flash_2ph!","text":"flash_2ph!(storage, eos, c, K; <keyword arguments>)\n\nNon-allocating version of flash_2ph where storage is pre-allocated.\n\nUseful if you are performing many flashes of the same system with varying conditions.\n\nArguments\n\nstorage: Should be output from flash_storage(eos, c, method = method). Preallocated storage.\n\nRemaining arguments documented in flash_2ph.\n\nKeyword arguments\n\nupdate_forces = true: Update the p, T dependent forces in storage initially.\n\n\n\n\n\n","category":"function"},{"location":"#MultiComponentFlash.flash_2ph-Union{Tuple{T}, Tuple{Any, T}, Tuple{Any, T, Any}} where T","page":"Home","title":"MultiComponentFlash.flash_2ph","text":"flash_2ph(eos, c, [K]; <keyword arguments>)\n\nPerform two-phase flash with a given EOS under a set of specific conditions. Returns vapor fraction. Modifies K in-place.\n\nGiven a mixture with pressure, temperature and mole fractions, this routine performs a vapor-liquid equilibrium calculation after a stability test.\n\nTwo outcomes are possible:\n\nA single-phase condition prevails (returned vapor fraction is NaN) and a single phase (liquid or vapor) is stable.\nA two-phase condition is possible. The routine produces K-values and vapor fraction so that the following holds:\nIsofugacity constraint for all components (f_li = f_vi)\nMolar balance for all components ((1-V) x_i - V y_i - z_i)\nUnity condition (sum_i (x_i - y_i) = 0)\n\nArguments\n\neos: the equation-of-state to be used for the flash\nc: conditions to flash the mixture at on the form (p = 10e5, T = 303.15, z = [0.5, 0.3, 0.2])\nK: optionally a buffer of length number_of_components(eos) used to hold K-values. Modified in-place.\n\nKeyword arguments\n\nmethod = SSIFlash(): Flash method to use. Can be SSIFlash(), NewtonFlash() or SSINewtonFlash().\ntolerance = 1e-8: Tolerance for the convergence criterion. Relative to 1-R_i_infty where R_i = f_ilf_iv\nmaxiter = 10000: Maximum nubmer of iterations for both stability tests and the main flash.\nverbose = false: Emit extra information during solve.\nextra_out = false: Return (V, K, conv) where conv contains iterations and oncergence status instead of just V.\n\nSee also: flash_2ph!, single_phase_label\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.flash_storage","page":"Home","title":"MultiComponentFlash.flash_storage","text":"flash_storage(eos, [c]; <keyword arguments>)\n\nPre-allocate storage for flash_2ph!.\n\nArguments\n\neos: the equation-of-state to be used for the flash\nc: conditions for the mixture. Types in conditions must match later usage (e.g. for use with ForwardDiff).\n\nKeyword arguments\n\nmethod = SSIFlash(): Flash method to use. Can be SSIFlash(), NewtonFlash() or SSINewtonFlash().\nstatic_size = false: Use SArrays and MArrays for fast flash, but slower compile times.\ninc_jac: Allocate storage for Newton/Jacobian. Required for Newton (and defaults to true for that method) or for diff_externals.\ndiff_externals = false: Allocate storage for matrix inversion required to produce partial derivatives of flash using set_partials.\n\nSee also: flash_2ph! set_partials\n\n\n\n\n\n","category":"function"},{"location":"#MultiComponentFlash.force_coefficients-Tuple{MultiComponentFlash.AbstractCubicEOS, Any}","page":"Home","title":"MultiComponentFlash.force_coefficients","text":"force_coefficients(eos, cond)\n\nGet coefficients for forces for a specific EOS (component interactions). For most cubics, these are a set of attractive (linear and quadratic) forces and a set of linear repulsive forces.\n\nNote that the current implementation of flash assumes that these are independent of the  compositions themselves.\n\nSee also force_coefficients!\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.force_scalars-Tuple{MultiComponentFlash.AbstractCubicEOS, Any, Any}","page":"Home","title":"MultiComponentFlash.force_scalars","text":"force_scalars(eos, cond, forces)\n\nCompute EOS specific scalars for the current conditions based on the forces.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.initial_guess_K!-Tuple{Any, Any, Any}","page":"Home","title":"MultiComponentFlash.initial_guess_K!","text":"initial_guess_K(eos, cond)\n\nIn-place version of initial_guess_K.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.initial_guess_K-Tuple{Any, Any}","page":"Home","title":"MultiComponentFlash.initial_guess_K","text":"initial_guess_K(eos, cond)\n\nProduce a plausible initial guess for K values for eos under current cond.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.inverse_flash_update!-NTuple{4, Any}","page":"Home","title":"MultiComponentFlash.inverse_flash_update!","text":"inverse_flash_update!(storage, eos, c, V)\n\nUpdate internal matrix of partial derivatives for a converged flash result.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.lbc_viscosities-Tuple{Any, Any, Any, FlashedMixture2Phase}","page":"Home","title":"MultiComponentFlash.lbc_viscosities","text":"lbc_viscosities(eos, p, T, flashed_mixture)\n\nCompute phase viscosities for a flashed two-phase mixture using the LBC correlation.\n\nAlways returns a named tuple of (μl, μv), even if the mixture is single-phase.\n\nThe value in the absent phase will be zero.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.lbc_viscosity-NTuple{4, Any}","page":"Home","title":"MultiComponentFlash.lbc_viscosity","text":"lbc_viscosity(eos, p, T, ph; <keyword arguments>)\n\nCompute the viscosity of a mixture using the Lohrenz-Bray-Clark correlation.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.liquid_mole_fraction-Tuple{Any, Any, Any}","page":"Home","title":"MultiComponentFlash.liquid_mole_fraction","text":"Compute liquid mole fraction from overall mole fraction, K-value and V vapor fraction\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.mass_densities-Tuple{Any, Any, Any, FlashedMixture2Phase}","page":"Home","title":"MultiComponentFlash.mass_densities","text":"mass_densities(eos, p, T, flashed_mixture)\n\nCompute mass densities for a flashed two-phase mixture.\n\nAlways returns a named tuple of (ρl, ρv), even if the mixture is single-phase.\n\nThe value in the absent phase will be zero.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.mass_density-Tuple{Any, Any, Any, FlashedPhase}","page":"Home","title":"MultiComponentFlash.mass_density","text":"Compute mass density of a flashed phase\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.mixture_fugacities!","page":"Home","title":"MultiComponentFlash.mixture_fugacities!","text":"Compute fugacities for all components.\n\n\n\n\n\n","category":"function"},{"location":"#MultiComponentFlash.mixture_fugacities-Tuple{Any, Vararg{Any, N} where N}","page":"Home","title":"MultiComponentFlash.mixture_fugacities","text":"Allocating version of mixture_fugacities!\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.molar_volume-Tuple{Any, Any, Any, FlashedPhase}","page":"Home","title":"MultiComponentFlash.molar_volume","text":"Compute molar volume of a flashed phase\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.molar_weight-Union{Tuple{MolecularProperty{R}}, Tuple{R}} where R","page":"Home","title":"MultiComponentFlash.molar_weight","text":"Get the molar weight for a species.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.number_of_components-Tuple{MultiComponentFlash.AbstractEOS}","page":"Home","title":"MultiComponentFlash.number_of_components","text":"number_of_components(eos)\n\nReturn number of components for the underlying mixture of the EOS.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.number_of_components-Union{Tuple{MultiComponentMixture{R, N}}, Tuple{N}, Tuple{R}} where {R, N}","page":"Home","title":"MultiComponentFlash.number_of_components","text":"number_of_components(mixture)\n\nReturn number of components in the MultiComponentMixture.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.phase_is_present-Tuple{Any, Any}","page":"Home","title":"MultiComponentFlash.phase_is_present","text":"phase_is_present(label, phase_state)\n\nCheck if a phase (symbol :liquid/:vapor) is present with the provided phase state.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.phase_saturations-Tuple{FlashedMixture2Phase}","page":"Home","title":"MultiComponentFlash.phase_saturations","text":"phase_saturations(flashed_mixture)\n\nCompute phase saturations for a flashed two-phase mixture.\n\nAlways returns a named tuple of (Sl, Sv), even if the mixture is single-phase.\n\nThe value in the absent phase will be zero.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.single_phase_label-Tuple{Any, Any}","page":"Home","title":"MultiComponentFlash.single_phase_label","text":"single_phase_label(mixture, cond)\n\nLi's method for single-phase labeling of a mixture. Estimate of pure vapor/liquid.\n\nReturns a vapor fraction that is either 1.0 (=pure vapor) or 0.0 (=pure liquid).\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.solve_rachford_rice","page":"Home","title":"MultiComponentFlash.solve_rachford_rice","text":"solve_rachford_rice(K, z, [V]; <keyword arguments>)\n\nCompute vapor mole fraction V for given equilibrium constants K and mole fractions z.\n\nArguments\n\nK - Equal length to z, containing the equilibrium constants for each component. z - Mole fractions. Should sum up to unity.\n\nKeyword arguments\n\ntol = 1e-12 - Tolerance for solve. maxiter - Maximum number of iterations ad - Use automatic differentiation (ForwardDiff) instead of analytical gradient.\n\nExamples\n\njulia> solve_rachford_rice([0.5, 1.5], [0.3, 0.7])\n0.8000000000000002\n\n\n\n\n\n","category":"function"},{"location":"#MultiComponentFlash.solve_roots-Tuple{MultiComponentFlash.AbstractCubicEOS, Any}","page":"Home","title":"MultiComponentFlash.solve_roots","text":"Specialization of solve_rotos for cubics\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.update_attractive_linear!-Tuple{Any, GenericCubicEOS, Any}","page":"Home","title":"MultiComponentFlash.update_attractive_linear!","text":"Update linear part of attractive term\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.update_attractive_quadratic!-Tuple{Any, Any, MultiComponentFlash.AbstractCubicEOS, Any}","page":"Home","title":"MultiComponentFlash.update_attractive_quadratic!","text":"Update quadratic part of attractive term\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.vapor_mole_fraction-Tuple{Any, Any, Any}","page":"Home","title":"MultiComponentFlash.vapor_mole_fraction","text":"Compute vapor mole fraction from overall mole fraction, K-value and V vapor fraction\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.vapor_mole_fraction-Tuple{Any, Any}","page":"Home","title":"MultiComponentFlash.vapor_mole_fraction","text":"Compute vapor mole fraction from liquid mole fraction and K-value\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.wilson_estimate!-Union{Tuple{R}, Tuple{AbstractVector{R}, Any, R, R}} where R<:Real","page":"Home","title":"MultiComponentFlash.wilson_estimate!","text":"wilson_estimate!(K, properties, p, T)\n\nUpdate a vector K in-place with K-values from wilson_estimate.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.wilson_estimate-Tuple{MultiComponentMixture, Any, Any}","page":"Home","title":"MultiComponentFlash.wilson_estimate","text":"wilson_estimate(properties, p, T)\n\nCreate vector of K-values that holds the wilson_estimate for each species.\n\n\n\n\n\n","category":"method"},{"location":"#MultiComponentFlash.wilson_estimate-Union{Tuple{R}, NTuple{5, R}} where R<:Real","page":"Home","title":"MultiComponentFlash.wilson_estimate","text":"wilson_estimate(p, T, ω, p_c, T_c)\n\nEstimate K-values for a given acentric factor ω and pressure and temperature at current and critical conditions.\n\nReference Vapor-Liquid Equilibrium. XI. A New Expression for the Excess Free Energy of Mixing by GM Wilson\n\n\n\n\n\n","category":"method"}]
}
