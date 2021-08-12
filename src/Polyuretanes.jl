module Polyuretanes

using Catalyst
using DifferentialEquations
using DelimitedFiles

export fit
export simulate
export plot_OH
export plot_all

function score(sim, experimental_data, sim_col)
    p = getindex.(sim.u, sim_col)
    score = 0.
    for (t, val) in eachrow(experimental_data)
        i = findfirst(ts -> ts > t, sim.t)
        dpdt = (p[i] - p[i - 1]) / (sim.t[i] - sim.t[i - 1])
        pred_at_t = p[i - 1] + dpdt * (t - sim.t[i - 1]) 
        score += (val - pred_at_t)^2
    end
    return score / (size(experimental_data, 1))
end

function optimizer(x, f;max_failed_trials=10000,showprogress=false)
    xbest = copy(x)
    xtrial = copy(x)
    failed_trials = 0
    step = 1.
    minstep = 1e-8
    fbest = f(x)
    itrial = 0
    irand = 1
    while failed_trials < max_failed_trials
        itrial += 1
        xtrial[irand] = xbest[irand] + -step * xbest[irand] + (2 * rand() * step) * xbest[irand]
        ftrial = f(xtrial)
        if ftrial < fbest
            failed_trials = 0
            xbest .= xtrial
            fbest = ftrial
            showprogress && println("Trial $itrial: fbest = $fbest step = $step")
            step = 1.
        else 
            irand = rand(1:length(x))
            step = max(minstep, step / 2)
            failed_trials += 1
        end
    end
    return fbest, xbest
end

function simulate(p0, c0, tspan, reaction_network)
    x = [values(p0)...]
    c0 = [values(c0)...]
    ode_problem = ODEProblem(reaction_network, c0, tspan, x)
    sim = solve(ode_problem, p=x)
    return sim
end

vap(OH,A,b) = A * exp(OH / b)

function reaction1()
    reaction = @reaction_network begin
        (k1, km1), U <--> NCO + DIPA_l
        k2, NCO + OH --> POL
        vap(OH, A, b), DIPA_l --> DIPA_v
    end k1 km1 k2 A b
    return reaction 
end
    
function reaction2()
    reaction = @reaction_network begin
        (k1, km1), U <--> NCO + DIPA_l
        2.2e-4, NCO + OH --> POL
        vap(OH, A, b), DIPA_l --> DIPA_v
    end k1 km1 A b
    return reaction 
end

struct Result
    sys
    sim
    score
    xbest
end

function Base.show(io::IO, r::Result)
    println("Score = ", r.score)
    println("Parameters: ")
    for val in r.xbest
        println(val[1], " = ", val[2])
    end
end

function plot_OH(result;title=nothing) 
    plot = Main.plot
    plot! = Main.plot!
    scatter! = Main.scatter!
    scalefontsizes = Main.scalefontsizes
    default = Main.default
  
    sim = result.sim
    experimental_data = result.sys.experimental_data
    default(
        fontfamily="Computer Modern",
        linewidth=2, framestyle=:box, label=nothing, grid=false
    )
    scalefontsizes()
    scalefontsizes(1.2)
  	plt = plot()
  	plot!(plt,
  	  	sim.t,getindex.(sim.u, 4),
  	  	xlabel="time / s",
  	  	ylabel="[OH] / mol/L",
  	  	label="Simulation"
  	)
  	scatter!(plt,
  		experimental_data[:,1],
  		experimental_data[:,2],
  		label="Experimental"
  	)
    plot!(plt, title=title)
  	return plt
end

function plot_all(
    results::AbstractVector{Result};
    title="",
    col=4,
    ylabel="[OH] / mol/L",
    xlabel="time / hours",
    tscale=3600
)
    plot = Main.plot
    plot! = Main.plot!
    scatter! = Main.scatter!
    scalefontsizes = Main.scalefontsizes
    default = Main.default
    default(
        fontfamily="Computer Modern",
        linewidth=2, framestyle=:box, label=nothing, grid=false
    )
    scalefontsizes()
    scalefontsizes(1.2)

    plt = plot()
    for (i, r) in pairs(results)
    sim = r.sim
        plot!(plt, sim.t / tscale, getindex.(sim.u, col), label=r.sys.label, color=i)
        scatter!(plt, r.sys.experimental_data[:,1] / tscale, r.sys.experimental_data[:,2], color=i)
    end
    plot!(plt, title=title)
    plot!(plt, xlabel=xlabel, ylabel=ylabel)

    return plt
end

Base.@kwdef struct ReactionSystem{C,P,T}
    title::String
    label::String
    experimental_data::Matrix{Float64}
    c0::C
    p0::P
    tspan::T
end

function datasets(dir="/home/leandro/.julia/dev/Polyuretanes/data")

    BD_IPDI_110C_100 = ReactionSystem(
        title="BD IPDI 110C 1",
        label="NCO/DIPA = 1",
        experimental_data=readdlm("$dir/bd_ipdi_110C_1.0.dat"),
        c0=(U = 6.68, NCO = 0., DIPA_l = 0., OH = 6.68, POL = 0., DIPA_v = 0.),
        p0=(k1 = 6.8e-6, km1 = 4.5e-4, k2 = 7.7e-5, A = 1.13e-6, b = 0.5),
        tspan=(0., 350e3)
    )
  
    BD_IPDI_110C_050 = ReactionSystem(
        title="BD IPDI 110C 0.5",
        label="NCP/DIPA = 0.5",
        experimental_data=readdlm("$dir/bd_ipdi_110C_0.5.dat"),
        c0=(U = 3.34, NCO = 3.34, DIPA_l = 0.,	OH = 6.68, POL = 0., DIPA_v = 0.),
        p0=(k1 = 6.8e-6, km1 = 4.5e-4, k2 = 7.7e-5, A = 1.13e-6, b = 0.5),
        tspan=(0., 350e3)
    )
      
    BD_IPDI_110C_025 = ReactionSystem(
        title="BD IPDI 110C 0.25",
        label="NCO/DIPA = 0.25",
        experimental_data=readdlm("$dir/bd_ipdi_110C_0.25.dat"),
        c0=(U = 1.67, NCO = 5.01, DIPA_l = 0.,	OH = 6.68, POL = 0., DIPA_v = 0.),
        p0=(k1 = 6.8e-6, km1 = 4.5e-4, k2 = 7.7e-5, A = 1.13e-6, b = 0.5),
        tspan=(0., 350e3)
    ) 
      
    return [BD_IPDI_110C_100, BD_IPDI_110C_050, BD_IPDI_110C_025] 

end

function objective_function(x, odes, systems, sim_col)
    f = 0.
    for i in 1:length(systems)
        ode_problem = odes[i]
        sim = solve(ode_problem, p=x, verbose=false)
        if sim.retcode == :Success
            f += score(sim, systems[i].experimental_data, sim_col)
        else
            return +Inf
        end
    end
    return f / length(systems)
end

fit(
    reaction=reaction1(), 
    system::ReactionSystem=datasets()[1],
    max_failed_trials=10000, showprogress=true
) = 
    fit(reaction, [system], max_failed_trials=max_failed_trials, showprogress=showprogress)

function fit(
    reaction=reaction1(),
    systems=datasets("/home/leandro/.julia/dev/Polyuretanes/data");
    max_failed_trials=10000,showprogress=true
)

    odes = []
    for system in systems
        c0 = [values(system.c0)...]
        x = Float64[]
        for par in params(reaction)
            push!(x, getfield(system.p0, Symbol(par)))
        end
        tspan = system.tspan
        ode_problem = ODEProblem(reaction, c0, tspan, x)
        push!(odes, ode_problem)
    end
    x0 = Float64[]
    for par in params(reaction)
        push!(x0, getfield(systems[1].p0, Symbol(par)))
    end
    @. x0 = x0 + -x0 + rand() * x0

    f(x) = objective_function(x, odes, systems, 4)
    fbest, xbest = optimizer(x0, f, max_failed_trials=max_failed_trials, showprogress=showprogress)

    r = Result[]
    for system in systems
        c0 = [values(system.c0)...]
        tspan = system.tspan
        sim = simulate(xbest, c0, tspan, reaction)
        sim_score = score(sim, system.experimental_data, 4)
        push!(
            r,Result(system, sim, sim_score,
            ntuple(i -> Symbol(params(reaction)[i]) => xbest[i], length(xbest)))
        )
    end

    return r, fbest
end


end # module