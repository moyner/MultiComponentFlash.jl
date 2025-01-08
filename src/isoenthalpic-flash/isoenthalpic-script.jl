include("enthalpy_flash_util.jl")
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
cp = [1 1 1 1]

# TT = cp + cp
# ref
# Tref = cp

qq = enthalpy_phase(eos,2,1,cp)


# let index = 1
# while index < 300
#     println("Counter is: $index")
#     index += 1
# end
# end

# while index1 < 1000
#     println("Counter is: $index1")
#     index1 += 1
# end
