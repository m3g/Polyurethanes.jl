module Polyurethanes

using DocStringExtensions
using Plots, Plots.Measures
using Catalyst
using DifferentialEquations
using Printf
using LaTeXStrings

export fit
export fitall
export simulate
export plot, savefig, writecsv
export System

import DelimitedFiles: readdlm
export readdlm

struct Result
    sys
    sim
    score
    xbest
    fbest
end

function Base.show(io::IO, r::Result)
    println(io,"System = \"$(r.sys.label)\"")
    println(io,"Score = ", r.score)
    println(io,"Parameters: ")
    for val in r.xbest
        println(io,val[1], " = ", val[2])
    end
end

Base.@kwdef struct System{C,P}
    title::String
    label::String
    experimental_data::Matrix{Float64}
    c0::C
    p0::P
    lower::P
    upper::P
end

"""

$(TYPEDEF)

$(TYPEDFIELDS)

Optimizer options.

"""
Base.@kwdef struct FitOptions
    max_failed_trials::Int = 1_000
    max_initial_trials::Int = 10_000
    showprogress::Bool = true
end

"""

```
score(sim, experimental_data, sim_col)
```

The score of these fits is the average squared *relative* deviation of a 
each predicted value relative to a measured value. That is, for each 
observation `val` we compute `((val - pred_at_t)/val)^2` where `pred_at_t`
is thre predicted value from the simulation. The difference is divided 
by `val` such that data spawning different orders of magnitude does not
have a fit dominated by the larger values. The score is the average
relative deviation of all data points. 

"""
function score(sim, experimental_data, sim_col)
    p = getindex.(sim.u, sim_col)
    score = 0.
    for (t, val) in eachrow(experimental_data)
        i = findfirst(ts -> ts > t, sim.t)
        dpdt = (p[i] - p[i - 1]) / (sim.t[i] - sim.t[i - 1])
        pred_at_t = p[i - 1] + dpdt * (t - sim.t[i - 1]) 
        score += ((val - pred_at_t)/val)^2
    end
    score / (size(experimental_data, 1))
end

"""

```
optimizer(x, f, lower, upper; opt=FitOptions())
```

This is a very simple and slow, yet robust global optimizer. It just perturbs
the best point to obtain a new trial point, check if the function is properly
evaluated at the new point, and keeps the best point so far. 

"""
function optimizer(x, f, lower, upper; opt=FitOptions())
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
        ftrial = try 
            f(xtrial)
        catch
            +Inf
        end
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

function set_tspan(system)
    tmax = maximum(@view(system.experimental_data[:,1]))
    tspan = (0.,tmax + tmax/100)
    return tspan
end

fit(reaction, 
    system::System;
    opt=FitOptions()
) = fit(reaction, [system], opt=opt)

function fit(
    reaction,
    systems::Vector{<:System};
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
        tspan = set_tspan(system)
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
        tspan = set_tspan(system)
        sim = simulate(xbest, c0, tspan, reaction)
        sim_score = score(sim, system.experimental_data, 4)
        push!(r,
            Result(
                system, sim, sim_score,
                ntuple(i -> Symbol(params(reaction)[i]) => xbest[i], length(xbest)),
                fbest
            )
        )
    end
    if length(systems) == 1
        return r[1]
    else
        return r
    end
end

function fitall(
    reaction=reaction1(),
    systems::Vector{<:System}=datasets("/home/leandro/.julia/dev/Polyurethanes/data");
    opt=FitOptions()
)
    # Optimal overal fit
    println("----------------------------------------------")
    println("Fitting all systems...")
    println("----------------------------------------------")
    r = fit(reaction,systems,opt=opt)

    return r
end

function plot(
    result::Result,
    species="OH"::String;
    labels=["Simulation","Experimental"],
    xlabel="time / hours",
    ylabel=latexstring(raw"\textrm{["*"$species"*raw"] / mol~L^{-1}}"),
    title=result.sys.label,
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
    ylabel=latexstring(raw"\textrm{["*"$species"*raw"] / mol~L^{-1}}"),
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

function writecsv(
    result::Result, filename;
    species="OH"::String,
    ylabel="[$species] / mol/L",
    xlabel="time / hours",
    tscale=3600
) 
    icol = findfirst(x -> String(x) == species, keys(result.sys.c0))
    sim = result.sim
    title = replace(result.sys.title," " => "_")
    file = dirname(filename)*"/"*title*"-"*basename(filename)
    csv = open(file,"w")
    println(csv,"# $(result.sys.title) - $(result.sys.label)")
    println(csv,"# $xlabel, $ylabel")
    for i in eachindex(sim.t)
        t = sim.t[i]/tscale
        y = getindex(sim.u[i],icol)
        println(csv,"$t, $y")
    end
    close(csv)
end

function writecsv(
    results::AbstractVector{Result}, filename;
    kargs...
)
    for r in results
        writecsv(r,filename;kargs...)
    end
end

end # module