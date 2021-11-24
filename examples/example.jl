using Catalyst
using Polyurethanes

#
# Reaction without evaporation
#
reaction_noevap = @reaction_network begin
    k1, U --> NCO + DIPA_l
    k2, NCO + OH --> POL
end k1 k2

#
# Reaction with evaporation
#
vap(OH,A,b) = A * exp(OH / b)
reaction_evap = @reaction_network begin
    (k1, km1), U <--> NCO + DIPA_l
    k2, NCO + OH --> POL
    vap(OH, A, b), DIPA_l --> DIPA_v
end k1 km1 k2 A b
    
# Directory of the data
dir = "./data"

#
# BD systems, without evaporation
#
bd_noevap = [
    System(
        title="BD IPDI 110C 1",
        label="NCO/DIPA = 1.00",
        experimental_data=readdlm("$dir/bd_ipdi_110C_1.0.dat"),
        c0=(U = 3.4675, NCO = 0., DIPA_l = 0., OH = 3.4675, POL = 0.),
        p0=(k1 = 6.8e-6, k2 = 2.2e-4),
        lower=(k1 = 0., k2 = 0.),
        upper=(k1 = +Inf, k2 = +Inf)
    )
    System(
        title="BD IPDI 110C 0.5",
        label="NCO/DIPA = 0.50",
        experimental_data=readdlm("$dir/bd_ipdi_110C_0.5.dat"),
        c0=(U = 2.246154, NCO = 2.246154, DIPA_l = 0.,	OH = 4.49231, POL = 0.),
        p0=(k1 = 6.8e-6, k2 = 7.7e-5),
        lower=(k1 = 0., k2 = 0.),
        upper=(k1 = +Inf, k2 = +Inf)
    )
    System(
        title="BD IPDI 110C 0.25",
        label="NCO/DIPA = 0.25",
        experimental_data=readdlm("$dir/bd_ipdi_110C_0.25.dat"),
        c0=(U = 1.603468, NCO = 4.810405, DIPA_l = 0.,	OH = 6.41387, POL = 0.),
        p0=(k1 = 6.8e-6, k2 = 7.7e-5),
        lower=(k1 = 0., k2 = 0.),
        upper=(k1 = +Inf, k2 = +Inf)
    ) 
]

#
# BD systems with evaporation
#
bd_evap = [ 
    System(
        title="BD IPDI 110C 1",
        label="NCO/DIPA = 1.00",
        experimental_data=readdlm("$dir/bd_ipdi_110C_1.0.dat"),
        c0=(U = 3.4675, NCO = 0., DIPA_l = 0., OH = 3.4675, POL = 0., DIPA_v = 0.),
        p0=(k1 = 6.8e-6, km1 = 4.5e-4, k2 = 2.2e-4, A = 1.13e-6, b = 0.5),
        lower=(k1 = 0., km1 = 0., k2 = 0., A = 0., b = -Inf),
        upper=(k1 = +Inf, km1 = +Inf, k2 = +Inf, A = +Inf, b = +Inf)
    )
    System(
        title="BD IPDI 110C 0.5",
        label="NCO/DIPA = 0.50",
        experimental_data=readdlm("$dir/bd_ipdi_110C_0.5.dat"),
        c0=(U = 2.246154, NCO = 2.246154, DIPA_l = 0.,	OH = 4.49231, POL = 0., DIPA_v = 0.),
        p0=(k1 = 6.8e-6, km1 = 4.5e-4, k2 = 7.7e-5, A = 1.13e-6, b = 0.5),
        lower=(k1 = 0., km1 = 0., k2 = 0., A = 0., b = -Inf),
        upper=(k1 = +Inf, km1 = +Inf, k2 = +Inf, A = +Inf, b = +Inf)
    )
    System(
        title="BD IPDI 110C 0.25",
        label="NCO/DIPA = 0.25",
        experimental_data=readdlm("$dir/bd_ipdi_110C_0.25.dat"),
        c0=(U = 1.603468, NCO = 4.810405, DIPA_l = 0.,	OH = 6.41387, POL = 0., DIPA_v = 0.),
        p0=(k1 = 6.8e-6, km1 = 4.5e-4, k2 = 7.7e-5, A = 1.13e-6, b = 0.5),
        lower=(k1 = 0., km1 = 0., k2 = 0., A = 0., b = -Inf),
        upper=(k1 = +Inf, km1 = +Inf, k2 = +Inf, A = +Inf, b = +Inf)
    ) 
]
      
#
# PEG at 100oC, with evaporation
#
peg_100_evap = [
    System(
        title="PEG IPDI 100C 1",
        label="NCO/DIPA = 1.00",
        experimental_data=readdlm("$dir/peg_ipdi_100C_1.0.dat"),
        c0=(U = 2.121362, NCO = 0., DIPA_l = 0., OH = 2.12136, POL = 0., DIPA_v = 0.),
        p0=(k1 = 6.8e-6, km1 = 4.5e-4, k2 = 2.2e-4, A = 1.13e-6, b = 0.5),
        lower=(k1 = 0., km1 = 0., k2 = 0., A = 0., b = -Inf),
        upper=(k1 = +Inf, km1 = +Inf, k2 = +Inf, A = +Inf, b = +Inf)
    )
    System(
        title="PEG IPDI 100C 0.25",
        label="NCO/DIPA = 0.25",
        experimental_data=readdlm("$dir/peg_ipdi_100C_0.25.dat"),
        c0=(U = 0.627473, NCO = 1.882418, DIPA_l = 0., OH = 2.50989, POL = 0., DIPA_v = 0.),
        p0=(k1 = 6.8e-6, km1 = 4.5e-4, k2 = 2.2e-4, A = 1.13e-6, b = 0.5),
        lower=(k1 = 0., km1 = 0., k2 = 0., A = 0., b = -Inf),
        upper=(k1 = +Inf, km1 = +Inf, k2 = +Inf, A = +Inf, b = +Inf)
    )
]

#
# PEG at 110oC, with evaporation
#
peg_110_evap = [
    System(
        title="PEG IPDI 110C 1",
        label="NCO/DIPA = 1.00",
        experimental_data=readdlm("$dir/peg_ipdi_110C_1.0.dat"),
        c0=(U = 2.121362, NCO = 0., DIPA_l = 0., OH = 2.12136, POL = 0., DIPA_v = 0.),
        p0=(k1 = 6.8e-6, km1 = 4.5e-4, k2 = 2.2e-4, A = 1.13e-6, b = 0.5),
        lower=(k1 = 0., km1 = 0., k2 = 0., A = 0., b = -Inf),
        upper=(k1 = +Inf, km1 = +Inf, k2 = +Inf, A = +Inf, b = +Inf)
    )
    System(
        title="PEG IPDI 110C 0.25",
        label="NCO/DIPA = 0.25",
        experimental_data=readdlm("$dir/peg_ipdi_110C_0.25.dat"),
        c0=(U = 0.627473, NCO = 1.882418, DIPA_l = 0., OH = 2.50989, POL = 0., DIPA_v = 0.),
        p0=(k1 = 6.8e-6, km1 = 4.5e-4, k2 = 2.2e-4, A = 1.13e-6, b = 0.5),
        lower=(k1 = 0., km1 = 0., k2 = 0., A = 0., b = -Inf),
        upper=(k1 = +Inf, km1 = +Inf, k2 = +Inf, A = +Inf, b = +Inf)
    )
]

#
# PEG 100oC, without evaporation
#
peg_100_noevap = [
    System(
        title="PEG IPDI 100C 1",
        label="NCO/DIPA = 1.00",
        experimental_data=readdlm("$dir/peg_ipdi_100C_1.0.dat"),
        c0=(U = 2.121362, NCO = 0., DIPA_l = 0., OH = 2.12136, POL = 0.),
        p0=(k1 = 6.8e-6, k2 = 2.2e-4),
        lower=(k1 = 0., k2 = 0.),
        upper=(k1 = +Inf, k2 = +Inf)
    )
    System(
        title="PEG IPDI 100C 0.25",
        label="NCO/DIPA = 0.25",
        experimental_data=readdlm("$dir/peg_ipdi_100C_0.25.dat"),
        c0=(U = 0.627473, NCO = 1.882418, DIPA_l = 0., OH = 2.50989, POL = 0.),
        p0=(k1 = 6.8e-6, k2 = 2.2e-4),
        lower=(k1 = 0., k2 = 0.),
        upper=(k1 = +Inf, k2 = +Inf)
    )
]
      
#
# PEG 110oC without evaporation
#
peg_110_noevap = [
    System(
        title="PEG IPDI 110C 1",
        label="NCO/DIPA = 1.00",
        experimental_data=readdlm("$dir/peg_ipdi_110C_1.0.dat"),
        c0=(U = 2.121362, NCO = 0., DIPA_l = 0., OH = 2.12136, POL = 0.),
        p0=(k1 = 6.8e-6, k2 = 2.2e-4),
        lower=(k1 = 0., k2 = 0.),
        upper=(k1 = +Inf, k2 = +Inf)
    )
    System(
        title="PEG IPDI 110C 0.25",
        label="NCO/DIPA = 0.25",
        experimental_data=readdlm("$dir/peg_ipdi_110C_0.25.dat"),
        c0=(U = 0.627473, NCO = 1.882418, DIPA_l = 0., OH = 2.50989, POL = 0.),
        p0=(k1 = 6.8e-6, k2 = 2.2e-4),
        lower=(k1 = 0., k2 = 0.),
        upper=(k1 = +Inf, k2 = +Inf)
    )
]

function run_system(basename,systems,reaction)
    # Fit each set independently
    for (isys,system) in pairs(systems)
        println("Running: ", system.title)
        result = fit(reaction,system)
        plt = plot(result)
        savefig(plt,"./plots/$(basename)_$isys.png")
        writecsv(result,"./simulations/$(basename)_$isys.csv")
    end
    # Obtain best overall fit
    println("Running fit for all sets...")
    result = fitall(reaction,systems)
    plt = plot(result)
    savefig(plt,"./plots/$basename.png")
    writecsv(result,"./simulations/$(basename)_all.csv")
end

for (basename, system, reaction) in [ 
                            ("bd_evap", bd_evap, reaction_evap), 
                            ("bd_noevap", bd_noevap, reaction_noevap), 
                            ("peg_100_evap", peg_100_evap, reaction_evap), 
                            ("peg_100_noevap", peg_100_noevap, reaction_noevap), 
                            ("peg_110_evap", peg_110_evap, reaction_evap), 
                            ("peg_110_noevap", peg_110_noevap, reaction_noevap),
                          ]
    run_system(basename, system, reaction)
end

