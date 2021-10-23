function tabulated_properties()
    # Various properties. Taken from CoolProp.
    data = Dict{String, MolecularProperty}()
    data["R507A"] = MolecularProperty(0.0988592, 3704900, 343.765, 0.0002014492399, 0.286)
    data["RC318"] = MolecularProperty(0.2000312, 2777500, 388.38, 0.0003226451742, 0.355345)
    data["SES36"] = MolecularProperty(0.18485, 2849000, 450.7, 0.0003571428571, 0.352)
    data["SulfurDioxide"] = MolecularProperty(0.0640638, 7886587.6, 430.64, 0.0001220256254, 0.2561299362)
    data["SulfurHexafluoride"] = MolecularProperty(0.1460554192, 3754983, 318.7232, 0.0001967606348, 0.21)
    data["Toluene"] = MolecularProperty(0.09213842, 4126000, 591.75, 0.000315556958, 0.2657)
    data["trans-2-Butene"] = MolecularProperty(0.05610632, 4027300, 428.61, 0.0002373605507, 0.2100768344)
    data["Water"] = MolecularProperty(0.018015268, 22064000, 647.096, 5.594803727e-05, 0.3442920843)
    data["Xenon"] = MolecularProperty(0.131293, 5842000, 289.733, 0.000119047619, 0.00363)
    data["1-Butene"] = MolecularProperty(0.05610632, 4005100, 419.29, 0.0002358490566, 0.1918606474)
    data["Acetone"] = MolecularProperty(0.05807914, 4700000, 508.1, 0.0002127659574, 0.3071)
    data["Air"] = MolecularProperty(0.02896546, 3786000, 132.5306, 8.452513778e-05, 0.0335)
    data["Ammonia"] = MolecularProperty(0.01703026, 11333000, 405.4, 7.569004444e-05, 0.25601)
    data["Argon"] = MolecularProperty(0.039948, 4863000, 150.687, 7.458551158e-05, -0.00219)
    data["Benzene"] = MolecularProperty(0.0781118, 4894000, 562.02, 0.0002562788314, 0.2108369733)
    data["CarbonDioxide"] = MolecularProperty(0.0440098, 7377300, 304.1282, 9.411847707e-05, 0.22394)
    data["CarbonMonoxide"] = MolecularProperty(0.0280101, 3494000, 132.86, 9.216589862e-05, 0.0497)
    data["CarbonylSulfide"] = MolecularProperty(0.0600751, 6370000, 378.77, 0.0001349527665, 0.0978)
    data["cis-2-Butene"] = MolecularProperty(0.05610632, 4225500, 435.75, 0.0002356267672, 0.2023595859)
    data["CycloHexane"] = MolecularProperty(0.08415948, 4082400, 553.6, 0.0003101736973, 0.20926)
    data["Cyclopentane"] = MolecularProperty(0.0701329, 4571200, 511.72, 0.0002617801047, 0.2011125734)
    data["CycloPropane"] = MolecularProperty(0.042081, 5579700, 398.3, 0.0001627891671, 0.1305495664)
    data["D4"] = MolecularProperty(0.29661576, 1347215.356, 586.5, 0.0009587727709, 0.5981248268)
    data["D5"] = MolecularProperty(0.3707697, 1160000, 619.15, 0.001216000031, 0.658)
    data["D6"] = MolecularProperty(0.444924, 961000, 645.78, 0.001594162692, 0.736)
    data["Deuterium"] = MolecularProperty(0.0040282, 1679600, 38.34, 5.803830528e-05, -0.1362902741)
    data["Dichloroethane"] = MolecularProperty(0.098959, 5227585.103, 561.6, 0.0002309468822, 0.2685785009)
    data["DiethylEther"] = MolecularProperty(0.0741216, 3649016.897, 466.7, 0.0002807636367, 0.2819293312)
    data["DimethylCarbonate"] = MolecularProperty(0.0900779, 4908800, 557, 0.00025, 0.346)
    data["DimethylEther"] = MolecularProperty(0.04606844, 5336800, 400.378, 0.0001683501684, 0.196)
    data["Ethane"] = MolecularProperty(0.03006904, 4872200, 305.322, 0.0001458387816, 0.099)
    data["Ethanol"] = MolecularProperty(0.04606844, 6268000, 514.71, 0.0001686340641, 0.644)
    data["EthylBenzene"] = MolecularProperty(0.106165, 3622400, 617.12, 0.0003648282243, 0.304)
    data["Ethylene"] = MolecularProperty(0.02805376, 5041800, 282.35, 0.0001309454817, 0.0866)
    data["EthyleneOxide"] = MolecularProperty(0.04405256, 7304686.16, 468.92, 0.0001394700139, 0.210195489)
    data["Fluorine"] = MolecularProperty(0.03799681, 5172400, 144.414, 6.409023906e-05, 0.0449)
    data["HeavyWater"] = MolecularProperty(0.020027508, 21672051.48, 643.847, 5.625704971e-05, 0.3643005688)
    data["Helium"] = MolecularProperty(0.004002602, 227600, 5.1953, 5.515719801e-05, -0.385)
    data["HFE143m"] = MolecularProperty(0.10004, 3635000, 377.921, 0.0002151397849, 0.2888713657)
    data["Hydrogen"] = MolecularProperty(0.00201588, 1296400, 33.145, 6.448284756e-05, -0.219)
    data["HydrogenChloride"] = MolecularProperty(0.0364609, 8288160.904, 324.55, 8.474576271e-05, 0.1281892174)
    data["HydrogenSulfide"] = MolecularProperty(0.03408088, 9000000, 373.1, 9.813542689e-05, 0.1005)
    data["IsoButane"] = MolecularProperty(0.0581222, 3629000, 407.817, 0.0002577481153, 0.1835317832)
    data["IsoButene"] = MolecularProperty(0.05610632, 4009800, 418.09, 0.0002398081535, 0.1925934522)
    data["Isohexane"] = MolecularProperty(0.08617536, 3040000, 497.7, 0.0003683241252, 0.2797)
    data["Isopentane"] = MolecularProperty(0.07214878, 3378000, 460.35, 0.0003057169061, 0.2274)
    data["Krypton"] = MolecularProperty(0.083798, 5525000, 209.48, 9.216589862e-05, -0.00089)
    data["m-Xylene"] = MolecularProperty(0.106165, 3534600, 616.89, 0.0003752345216, 0.326)
    data["MD2M"] = MolecularProperty(0.310685, 1227000, 599.4, 0.001093300515, 0.668)
    data["MD3M"] = MolecularProperty(0.384839, 945000, 628.36, 0.001458154971, 0.722)
    data["MD4M"] = MolecularProperty(0.45899328, 877000, 653.2, 0.001650000016, 0.8246471473)
    data["MDM"] = MolecularProperty(0.23653146, 1410044.756, 564.09, 0.000921288245, 0.5280658491)
    data["Methane"] = MolecularProperty(0.0160428, 4599200, 190.564, 9.862781099e-05, 0.01142)
    data["Methanol"] = MolecularProperty(0.03204216, 8215850, 512.5, 0.0001173705495, 0.5720322)
    data["MethylLinoleate"] = MolecularProperty(0.29447206, 1341000, 799, 0.001237011381, 0.8054063871)
    data["MethylLinolenate"] = MolecularProperty(0.29245618, 1369000, 772, 0.001180219521, 1.142605259)
    data["MethylOleate"] = MolecularProperty(0.29648794, 1246000, 782, 0.001230239282, 0.90584936)
    data["MethylPalmitate"] = MolecularProperty(0.27045066, 1350000, 755, 0.001114827202, 0.910316178)
    data["MethylStearate"] = MolecularProperty(0.29850382, 1239000, 775, 0.001258970162, 1.01756)
    data["MM"] = MolecularProperty(0.16237752, 1939000, 518.75, 0.0006290000472, 0.418)
    data["n-Butane"] = MolecularProperty(0.0581222, 3796000, 425.125, 0.0002549219298, 0.2008100946)
    data["n-Decane"] = MolecularProperty(0.14228168, 2103000, 617.7, 0.0006097560976, 0.4884)
    data["n-Dodecane"] = MolecularProperty(0.17033484, 1817000, 658.1, 0.0007518796992, 0.5741822212)
    data["n-Heptane"] = MolecularProperty(0.100202, 2736000, 540.13, 0.0004319051724, 0.349)
    data["n-Hexane"] = MolecularProperty(0.08617536, 3034000, 507.82, 0.000369565829, 0.299)
    data["n-Nonane"] = MolecularProperty(0.1282551, 2281000, 594.55, 0.0005524861878, 0.4433)
    data["n-Octane"] = MolecularProperty(0.1142285, 2497000, 569.32, 0.0004862867146, 0.395)
    data["n-Pentane"] = MolecularProperty(0.07214878, 3370000, 469.7, 0.0003109861207, 0.251)
    data["n-Propane"] = MolecularProperty(0.04409562, 4251200, 369.89, 0.0002, 0.1521)
    data["n-Undecane"] = MolecularProperty(0.15630826, 1990400, 638.8, 0.00066010223, 0.5390371014)
    data["Neon"] = MolecularProperty(0.020179, 2680000, 44.4918, 4.187253999e-05, -0.03844929927)
    data["Neopentane"] = MolecularProperty(0.07214878, 3196000, 433.74, 0.0003058103976, 0.1961)
    data["Nitrogen"] = MolecularProperty(0.02801348, 3395800, 126.192, 8.941423556e-05, 0.0372)
    data["NitrousOxide"] = MolecularProperty(0.0440128, 7245000, 309.52, 9.737098345e-05, 0.1613)
    data["Novec649"] = MolecularProperty(0.3160438, 1869026.583, 441.81, 0.0005208333333, 0.4710235936)
    data["o-Xylene"] = MolecularProperty(0.106165, 3737500, 630.259, 0.0003725088471, 0.312)
    data["OrthoDeuterium"] = MolecularProperty(0.0040282, 1679600, 38.34, 5.803830528e-05, -0.1362902741)
    data["OrthoHydrogen"] = MolecularProperty(0.00201594, 1310650, 33.22, 6.474779953e-05, -0.219)
    data["Oxygen"] = MolecularProperty(0.0319988, 5043000, 154.581, 7.336757153e-05, 0.0222)
    data["p-Xylene"] = MolecularProperty(0.106165, 3531500, 616.168, 0.0003712062719, 0.324)
    data["ParaDeuterium"] = MolecularProperty(0.0040282, 1679600, 38.34, 5.803830528e-05, -0.1362902741)
    data["ParaHydrogen"] = MolecularProperty(0.00201588, 1285800, 32.938, 6.435834728e-05, -0.219)
    data["Propylene"] = MolecularProperty(0.04207974, 4555000, 364.211, 0.0001832508704, 0.146)
    data["Propyne"] = MolecularProperty(0.04006, 5626000, 402.38, 0.0001635769703, 0.204)
    data["R11"] = MolecularProperty(0.137368, 4394000, 471.06, 0.0002431292035, 0.1887506483)
    data["R113"] = MolecularProperty(0.187375, 3392200, 487.21, 0.0003345982143, 0.252535)
    data["R114"] = MolecularProperty(0.170921, 3257000, 418.83, 0.0002947070612, 0.2523)
    data["R115"] = MolecularProperty(0.154466416, 3129036.968, 353.1, 0.0002512562814, 0.2484349931)
    data["R116"] = MolecularProperty(0.13801182, 3048000, 293.03, 0.0002250225023, 0.2566)
    data["R12"] = MolecularProperty(0.120913, 4136100, 385.12, 0.0002140053097, 0.1794783173)
    data["R123"] = MolecularProperty(0.152931, 3672000, 456.831, 0.0002780545193, 0.281922497)
    data["R1233zd(E)"] = MolecularProperty(0.1304944, 3623637.776, 439.6, 0.0002717391304, 0.3024522936)
    data["R1234yf"] = MolecularProperty(0.1140415928, 3382200, 367.85, 0.0002398081535, 0.276)
    data["R1234ze(E)"] = MolecularProperty(0.1140415928, 3636250, 382.52, 0.0002331002331, 0.313)
    data["R1234ze(Z)"] = MolecularProperty(0.1140415928, 3533000, 423.27, 0.0002426416868, 0.3274)
    data["R124"] = MolecularProperty(0.1364762, 3624295, 395.425, 0.0002437075, 0.2880950842)
    data["R125"] = MolecularProperty(0.1200214, 3617700, 339.173, 0.0002092487968, 0.3052)
    data["R13"] = MolecularProperty(0.104459, 3879000, 301.88, 0.0001792114695, 0.1745863278)
    data["R134a"] = MolecularProperty(0.102032, 4059280, 374.21, 0.0001993201985, 0.32684)
    data["R13I1"] = MolecularProperty(0.1959104, 3952566.072, 396.44, 0.000225703065, 0.1761811178)
    data["R14"] = MolecularProperty(0.0880046, 3750000, 227.51, 0.0001406584622, 0.1785)
    data["R141b"] = MolecularProperty(0.11694962, 4212000, 477.5, 0.0002550369804, 0.2195)
    data["R142b"] = MolecularProperty(0.10049503, 4055000, 410.26, 0.0002253267237, 0.2321)
    data["R143a"] = MolecularProperty(0.084041, 3761000, 345.857, 0.0001949906892, 0.2614896462)
    data["R152A"] = MolecularProperty(0.066051, 4520000, 386.411, 0.000179486413, 0.2752171145)
    data["R161"] = MolecularProperty(0.0480595, 5010000, 375.25, 0.0001592356688, 0.2162428411)
    data["R21"] = MolecularProperty(0.1029227, 5181200, 451.48, 0.0001956654009, 0.2061)
    data["R218"] = MolecularProperty(0.18801933, 2640000, 345.02, 0.0002994011976, 0.3172)
    data["R22"] = MolecularProperty(0.086468, 4990000, 369.295, 0.0001650649861, 0.22082)
    data["R227EA"] = MolecularProperty(0.17002886, 2925242.234, 374.9, 0.0002861230329, 0.357641242)
    data["R23"] = MolecularProperty(0.07001385, 4832000, 299.293, 0.0001329787234, 0.2629648925)
    data["R236EA"] = MolecularProperty(0.1520384, 3420000, 412.44, 0.0002691065662, 0.3687823804)
    data["R236FA"] = MolecularProperty(0.1520384, 3200000, 398.07, 0.0002757859901, 0.3769117976)
    data["R245ca"] = MolecularProperty(0.13404794, 3940740.489, 447.57, 0.0002551020408, 0.3545654191)
    data["R245fa"] = MolecularProperty(0.13404794, 3651000, 427.01, 0.0002580645161, 0.3776)
    data["R32"] = MolecularProperty(0.052024, 5782000, 351.255, 0.0001226981129, 0.2769)
    data["R365MFC"] = MolecularProperty(0.14807452, 3266207.068, 460, 0.0003125, 0.3774450972)
    data["R40"] = MolecularProperty(0.05048752, 6671686.977, 416.3, 0.0001390002242, 0.1500684964)
    data["R404A"] = MolecularProperty(0.0976038, 3734800, 345.27, 0.0002024291498, 0.293)
    data["R407C"] = MolecularProperty(0.0862036, 4631700, 359.345, 0.0001901140668, 0.363)
    data["R41"] = MolecularProperty(0.03403292, 5897000, 317.28, 0.0001075268817, 0.2004)
    data["R410A"] = MolecularProperty(0.0725854, 4901200, 344.494, 0.0001581277672, 0.296)
    return data
end

"""
    cubic_benchmark(name)

Get a benchmark equation-of-state instance together with some data for testing.
"""
function cubic_benchmark(name)
    psi = 6.894757293168360e+03
    Rankine = 0.555555555555556
    bar = 1e5
    data = Dict()
    eos_type = PengRobinson()

    bic = nothing
    if name == "spe5"
        T_c = [343, 665.7, 913.4, 1111.8, 1270.0, 1380.0]*Rankine;
        p_c = [667.8, 616.3, 436.9, 304.0, 200.0, 162.0]*psi;
        mw = [16.040, 44.100, 86.180, 142.290, 206.000, 282.000]/1000;
        acc = [0.0130, 0.1524, 0.3007, 0.4885, 0.6500, 0.8500];
        Z_c = [0.290, 0.277, 0.264, 0.257, 0.245, 0.235];
        # Convert critical volumes to critical densities
        V_c = Z_c.*8.314.*T_c./p_c;
        names = ["C1", "C3", "C6", "C10", "C15", "C20"];
        ncomp = length(names);
        bic = zeros(ncomp, ncomp);
        bic[1, 5] = 0.05;
        bic[1, 6] = 0.05;
        bic[2, 5] = 0.005;
        bic[2, 6] = 0.005;
        bic = Symmetric(bic)
        props = MolecularProperty.(mw, p_c, T_c, V_c, acc)

        z0 = [0.77, 0.20, 0.03, 0, 0, 0]
        z = [0.5, 0.03, 0.07, 0.20, 0.15, 0.05]
        p = 4000*psi
        T = 344.26
        p0 = p
        T0 = T

        initial_cond = (z = z0, p = p0, T = T0)
        displ_cond = (z = z, p = p, T = T)
        data["initial_condition"] = initial_cond
        data["injection_condition"] = displ_cond
    elseif name == "simple"
        names = ["Methane", "CarbonDioxide", "n-Decane"]
        props = MolecularProperty.(names)

        p = 75*bar
        T = 273.15 + 150
        z0 = [0.3, 0.1, 0.6]
        z = [0.1, 0.9, 0.0]
        initial_cond = (z = z0, p = p, T = T)
        displ_cond = (z = z, p = p, T = T)
        data["initial_condition"] = initial_cond
        data["injection_condition"] = displ_cond
    else
        error("Unknown benchmark $name")
    end
    mixture = MultiComponentMixture(props, names = names, A_ij = bic)
    eos = GenericCubicEOS(mixture, eos_type)

    return (eos, data)
end