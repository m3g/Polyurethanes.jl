module Polyuretanes

using Plots, Plots.Measures
using Catalyst
using DifferentialEquations
using DelimitedFiles
using Printf

export fit
export fitall
export simulate
export plot

function score(sim, experimental_data, sim_col)
    p = getindex.(sim.u, sim_col)
    score = 0.
    for (t, val) in eachrow(experimental_data)
        i = findfirst(ts -> ts > t, sim.t)
        dpdt = (p[i] - p[i - 1]) / (sim.t[i] - sim.t[i - 1])
        pred_at_t = p[i - 1] + dpdt * (t - sim.t[i - 1]) 
        score += ((val - pred_at_t)/val)^2
    end
    return score / (size(experimental_data, 1))
end

function optimizer(x, f, lower, upper;opt=FitOptions())
    xbest = copy(x)
    xtrial = copy(x)
    failed_trials = 0
    step = 1.
    minstep = 1e-8
    fbest = f(x)
    itrial = 0
    irand = 1
    if opt.showprogress
        @printf("Trial: %10i fbest = %12.5e step = %12.5e", itrial, fbest, step)
    end
    while failed_trials < opt.max_failed_trials || itrial < opt.max_initial_trials
        itrial += 1
        xtrial[irand] = xbest[irand] + -step * xbest[irand] + (2 * rand() * step) * xbest[irand]
        xtrial[irand] = max(lower[irand],min(upper[irand],xtrial[irand]))
        ftrial = f(xtrial)
        if ftrial < fbest
            failed_trials = 0
            xbest .= xtrial
            fbest = ftrial
            if opt.showprogress
              @printf("\rTrial: %10i fbest = %12.5e step = %12.5e", itrial, fbest, step)
            end
            step = 1.
        else 
            irand = rand(1:length(x))
            step = max(minstep, step / 2)
            failed_trials += 1
        end
    end
    opt.showprogress && @printf("\n")
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

function plot(
    result::Result,
    species="OH"::String;
    labels=["Simulation","Experimental"],
    xlabel="time / hours",
    ylabel="[$species] / mol/L",
    title="",
    tscale=3600
) 

    icol = findfirst(x -> String(x) == species, keys(result.sys.c0))
  
    sim = result.sim
    experimental_data = result.sys.experimental_data
    default(
        fontfamily="Computer Modern",
        linewidth=2, framestyle=:box, label=nothing, grid=false,
        margin=5mm
    )
    scalefontsizes()
    scalefontsizes(1.2)
  	plt = Plots.plot()
  	plot!(plt,
  	  	sim.t/tscale,
        getindex.(sim.u, icol),
  	  	xlabel=xlabel,
  	  	ylabel=ylabel,
  	  	label=labels[1]
  	)
  	scatter!(plt,
  		experimental_data[:,1]/tscale,
  		experimental_data[:,2],
  		label=labels[2]
  	)
    plot!(plt, title=title)
  	return plt
end

function plot(
    results::AbstractVector{Result};
    species="OH",
    title="",
    ylabel="[$species] / mol/L",
    xlabel="time / hours",
    tscale=3600
)
    icol = findfirst(x -> String(x) == species, keys(results[1].sys.c0))

    default(
        fontfamily="Computer Modern",
        linewidth=2, framestyle=:box, label=nothing, grid=false,
        margin=5mm
    )
    scalefontsizes()
    scalefontsizes(1.2)

    plt = Plots.plot()
    for (i, r) in pairs(results)
    sim = r.sim
        plot!(plt, sim.t / tscale, getindex.(sim.u, icol), label=r.sys.label, color=i)
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
    lower::P
    upper::P
    tspan::T
end

function datasets(dir="/home/leandro/.julia/dev/Polyuretanes/data")

    BD_IPDI_110C_100 = ReactionSystem(
        title="BD IPDI 110C 1",
        label="NCO/DIPA = 1",
        experimental_data=readdlm("$dir/bd_ipdi_110C_1.0.dat"),
        c0=(U = 6.68, NCO = 0., DIPA_l = 0., OH = 6.68, POL = 0., DIPA_v = 0.),
        p0=(k1 = 6.8e-6, km1 = 4.5e-4, k2 = 2.2e-4, A = 1.13e-6, b = 0.5),
        lower=(k1 = -Inf, km1 = -Inf, k2 = 1.7e-4, A = -Inf, b = -Inf),
        upper=(k1 = +Inf, km1 = +Inf, k2 = 2.7e-4, A = +Inf, b = +Inf),
        tspan=(0., 350e3)
    )
  
    BD_IPDI_110C_050 = ReactionSystem(
        title="BD IPDI 110C 0.5",
        label="NCP/DIPA = 0.5",
        experimental_data=readdlm("$dir/bd_ipdi_110C_0.5.dat"),
        c0=(U = 3.34, NCO = 3.34, DIPA_l = 0.,	OH = 6.68, POL = 0., DIPA_v = 0.),
        p0=(k1 = 6.8e-6, km1 = 4.5e-4, k2 = 7.7e-5, A = 1.13e-6, b = 0.5),
        lower=(k1 = -Inf, km1 = -Inf, k2 = 1.7e-4, A = -Inf, b = -Inf),
        upper=(k1 = +Inf, km1 = +Inf, k2 = 2.7e-4, A = +Inf, b = +Inf),
        tspan=(0., 350e3)
    )
      
    BD_IPDI_110C_025 = ReactionSystem(
        title="BD IPDI 110C 0.25",
        label="NCO/DIPA = 0.25",
        experimental_data=readdlm("$dir/bd_ipdi_110C_0.25.dat"),
        c0=(U = 1.67, NCO = 5.01, DIPA_l = 0.,	OH = 6.68, POL = 0., DIPA_v = 0.),
        p0=(k1 = 6.8e-6, km1 = 4.5e-4, k2 = 7.7e-5, A = 1.13e-6, b = 0.5),
        lower=(k1 = -Inf, km1 = -Inf, k2 = 1.7e-4, A = -Inf, b = -Inf),
        upper=(k1 = +Inf, km1 = +Inf, k2 = 2.7e-4, A = +Inf, b = +Inf),
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

Base.@kwdef struct FitOptions
    max_failed_trials::Int = 10_000
    max_initial_trials::Int = 100_000
    showprogress::Bool = true
end

fit(reaction, 
    system::ReactionSystem;
    opt=FitOptions()
) = fit(reaction, [system], opt=opt)

function fit(
    reaction,
    systems::Vector{<:ReactionSystem};
    opt=FitOptions()
)

    println("Building differential equation system... ")
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
    lower = Float64[]
    upper = Float64[]
    for par in params(reaction)
        push!(x0, getfield(systems[1].p0, Symbol(par)))
        push!(lower, getfield(systems[1].lower, Symbol(par)))
        push!(upper, getfield(systems[1].upper, Symbol(par)))
    end
    @. x0 = x0 + -x0 + rand() * x0

    println("Setting up objective function... ")
    f(x) = objective_function(x, odes, systems, 4)
    fbest, xbest = optimizer(x0, f, lower, upper, opt=opt)

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

function fitall(
    reaction=reaction1(),
    systems::Vector{<:ReactionSystem}=datasets("/home/leandro/.julia/dev/Polyuretanes/data");
    opt=FitOptions()
)
    # Fit each system independently
    individual_results = []
    for system in systems
        println("----------------------------------------------")
        println("Fitting system $(system.title)...")
        println("----------------------------------------------")
        r = fit(reaction,system,opt=opt)
        push!(individual_results, r[1])
    end

    # Optimal overal fit
    println("----------------------------------------------")
    println("Fitting all systems...")
    println("----------------------------------------------")
    r = fit(reaction,systems,opt=opt)

    return individual_results, r
end

end # module