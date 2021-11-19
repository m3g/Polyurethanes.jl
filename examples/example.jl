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
bd_systems_noevap = [
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

# Fit each set independently
bd_noevap_1 = fit(reaction_noevap,bd_systems_noevap[1])
plot(bd_noevap_1)
savefig("./plots/bd_noevap_1.png")

bd_noevap_2 = fit(reaction_noevap,bd_systems_noevap[2])
plot(bd_noevap_2)
savefig("./plots/bd_noevap_2.png")

bd_noevap_3 = fit(reaction_noevap,bd_systems_noevap[3])
plot(bd_noevap_3)
savefig("./plots/bd_noevap_3.png")

# Obtain best overall fit
bd_noevap = fitall(reaction_noevap,bd_systems_noevap)
plot(bd_noevap)
savefig("./plots/bd_noevap.png")

#
# BD systems with evaporation
#
bd_systems_evap = [ 
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
      
# Fit each set independently
bd_evap_1 = fit(reaction_evap,bd_systems_evap[1])
plot(bd_evap_1)
savefig("./plots/bd_evap_1.png")

bd_evap_2 = fit(reaction_evap,bd_systems_evap[2])
plot(bd_evap_2)
savefig("./plots/bd_evap_2.png")

bd_evap_3 = fit(reaction_evap,bd_systems_evap[3])
plot(bd_evap_3)
savefig("./plots/bd_evap_3.png")

# Obtain best overall fit
bd_evap = fitall(reaction_evap,bd_systems_evap)
plot(bd_evap)
savefig("./plots/bd_evap.png")
      
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

# Fit each set independently
peg100_evap_1 = fit(reaction_evap,peg_100_evap[1])
plot(peg100_evap_1)
savefig("./plots/peg100_evap_1.png")

peg100_evap_2 = fit(reaction_evap,peg_100_evap[2])
plot(peg100_evap_2)
savefig("./plots/peg100_evap_2.png")

# Obtain best overall fit
peg100_evap = fitall(reaction_evap,peg_100_evap)
plot(peg100_evap)
savefig("./plots/peg100_evap.png")
      
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

# Fit each set independently
peg110_evap_1 = fit(reaction_evap,peg_110_evap[1])
plot(peg110_evap_1)
savefig("./plots/peg110_evap_1.png")

peg110_evap_2 = fit(reaction_evap,peg_110_evap[2])
plot(peg110_evap_2)
savefig("./plots/peg110_evap_2.png")

# Obtain best overall fit
peg110_evap = fitall(reaction_evap,peg_110_evap)
plot(peg110_evap)
savefig("./plots/peg110_evap.png")
      
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

# Fit each set independently
peg100_noevap_1 = fit(reaction_noevap,peg_100_noevap[1])
plot(peg100_noevap_1)
savefig("./plots/peg100_noevap_1.png")

peg100_noevap_2 = fit(reaction_noevap,peg_100_noevap[2])
plot(peg100_noevap_2)
savefig("./plots/peg100_noevap_2.png")

# Obtain best overall fit
peg100_noevap = fitall(reaction_noevap,peg_100_noevap)
plot(peg100_noevap)
savefig("./plots/peg100_noevap.png")
      
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

# Fit each set independently
peg110_noevap_1 = fit(reaction_noevap,peg_110_noevap[1])
plot(peg110_noevap_1)
savefig("./plots/peg110_noevap_1.png")

peg110_noevap_2 = fit(reaction_noevap,peg_110_noevap[2])
plot(peg110_noevap_2)
savefig("./plots/peg110_noevap_2.png")

# Obtain best overall fit
peg110_noevap = fitall(reaction_noevap,peg_110_noevap)
plot(peg110_noevap)
savefig("./plots/peg110_noevap.png")
