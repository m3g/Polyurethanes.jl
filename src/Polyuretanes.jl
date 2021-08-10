module Polyuretanes

using Catalyst
using DifferentialEquations
using DelimitedFiles

export fit
export simulate
export plot_OH
export plot_all

function score(sim,experimental_data,sim_col)
  score = 0.
  p = getindex.(sim.u,sim_col)
  for (t,val) in eachrow(experimental_data)
    i = findfirst(ts -> ts > t, sim.t)
    dpdt = (p[i]-p[i-1])/(sim.t[i]-sim.t[i-1])
    pred_at_t = p[i-1] + dpdt * (t-sim.t[i-1]) 
    score += (val - pred_at_t)^2
  end
  return score
end

function objective_function(x,ode_problem,experimental_data,sim_col)
  f = try 
    sim = solve(ode_problem,p=x,verbose=false)
	  score(sim,experimental_data,sim_col)
  catch
    +Inf
  end
  return f
end

function fit_parameters(p0,c0,tspan,reaction_network,experimental_data;
  max_failed_trials=10000,
  showprogress=true
)
	x = [values(p0)...]
	c0 = [values(c0)...]
	ode_problem = ODEProblem(reaction_network,c0,tspan,x)
	# the fourth column corresponds to the OH concentration
	f(x) = objective_function(x,ode_problem,experimental_data,4)
  fbest, x = optimizer(f,x,max_failed_trials=max_failed_trials,showprogress=showprogress)
  return fbest, x
end

function optimizer(x,f;max_failed_trials=10000,showprogress=false)
  xtrial = copy(x)
  failed_trials = 0
  step = 0.1
  fbest = f(x)
  while failed_trials < max_failed_trials
#    irand = rand(1:length(x))
#    xtrial[irand] = x[irand] + -step*x[irand] + 2*step*x[irand]*rand()
    @. xtrial = x + -step*x + (2*rand()*step)*x
    ftrial = f(xtrial)
    if ftrial < fbest
      failed_trials = 0
      x .= xtrial
      fbest = ftrial
      showprogress && @show fbest, step
    else 
      failed_trials += 1
    end
  end
  return fbest, x
end

function simulate(p0,c0,tspan,reaction_network)
	x = [values(p0)...]
	c0 = [values(c0)...]
	ode_problem = ODEProblem(reaction_network,c0,tspan,x)
  sim = solve(ode_problem,p=x)
  return sim
end

vap(OH,A,b) = A*exp(OH/b)

function reaction1()
  reaction = @reaction_network begin
    (k1,km1), U <--> NCO + DIPA_l
    k2, NCO + OH --> POL
    vap(OH,A,b), DIPA_l --> DIPA_v
  end k1 km1 k2 A b
  return reaction 
end

function reaction2()
  reaction = @reaction_network begin
    (k1,km1), U <--> NCO + DIPA_l
    2.2e-4, NCO + OH --> POL
    vap(OH,A,b), DIPA_l --> DIPA_v
  end k1 km1 A b
  return reaction 
end

struct Result
  sys
  sim
  score
  xbest
end

function Base.show(io::IO,r::Result)
  println("Score = ", r.score)
  println("Parameters: ")
  for val in r.xbest
    println(val[1]," = ", val[2])
  end
end

function fit(
  sys;
  reaction_network = reaction1(),
  max_failed_trials=10000,
  showprogress=true
)
  p0 = sys.p0
  c0 = sys.c0
  tspan = sys.tspan
  experimental_data = sys.experimental_data
  fbest, xbest = fit_parameters(
    p0,c0,tspan,reaction_network,experimental_data,
    max_failed_trials=max_failed_trials,
    showprogress=showprogress
  )
  sim = simulate(xbest,c0,tspan,reaction_network)
  return Result(sys, sim, fbest, ntuple(i -> Symbol(params(reaction_network)[i]) => xbest[i], length(xbest)))
end

function plot_OH(result;title=nothing) 
  plot = Main.plot
  plot! = Main.plot!
  scatter! = Main.scatter!
  scalefontsizes = Main.scalefontsizes
  default = Main.default

  sim = result.sim
  experimental_data = result.sys.experimental_data
  default(fontfamily="Computer Modern",
        linewidth=2, framestyle=:box, label=nothing, grid=false)
  scalefontsizes()
  scalefontsizes(1.2)
	plt = plot()
	plot!(plt,
		sim.t,getindex.(sim.u,4),
		xlabel="time / s",
		ylabel="[OH] / mol/L",
		label="Simulation"
	)
	scatter!(plt,
		experimental_data[:,1],
		experimental_data[:,2],
		label="Experimental"
	)
  plot!(plt,title=title)
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
  default(fontfamily="Computer Modern",
          linewidth=2, framestyle=:box, label=nothing, grid=false)
  scalefontsizes()
  scalefontsizes(1.2)

  plt = plot()
  for (i,r) in pairs(results)
    sim = r.sim
    plot!(plt,sim.t/tscale,getindex.(sim.u,col),label=r.sys.concentration,color=i)
    scatter!(plt,r.sys.experimental_data[:,1]/tscale,r.sys.experimental_data[:,2],color=i)
  end
  plot!(plt,title=title)
  plot!(plt,xlabel=xlabel,ylabel=ylabel)

  return plt
end

function datasets(dir)

  return (

    BD_IPDI_110C_100 = (
      title = "BD IPDI 110C 1.00 mol/L",
      concentration = "1.00 mol/L",
      experimental_data = readdlm("$dir/bd_ipdi_110C_1.0.dat"),
      c0 = ( U = 6.68, NCO = 0., DIPA_l = 0.,	OH = 6.68, POL = 0., DIPA_v = 0.),
      p0 = ( k1 = 6.8e-6, km1 = 4.5e-4, k2 = 7.7e-5, A = 1.13e-6, b = 0.5),
      tspan = (0.,350e3)
    ), 

    BD_IPDI_110C_050 = (
      title = "BD IPDI 110C 0.50 mol/L",
      concentration = "0.50 mol/L",
      experimental_data = readdlm("$dir/bd_ipdi_110C_0.5.dat"),
      c0 = ( U = 3.34, NCO = 3.34, DIPA_l = 0.,	OH = 6.68, POL = 0., DIPA_v = 0.),
      p0 = ( k1 = 6.8e-6, km1 = 4.5e-4, k2 = 7.7e-5, A = 1.13e-6, b = 0.5),
      tspan = (0.,350e3)
    ), 
    
    BD_IPDI_110C_025 = (
      title = "BD IPDI 110C 0.25 mol/L",
      concentration = "0.25 mol/L",
      experimental_data = readdlm("$dir/bd_ipdi_110C_0.25.dat"),
      c0 = ( U = 1.67, NCO = 5.01, DIPA_l = 0.,	OH = 6.68, POL = 0., DIPA_v = 0.),
      p0 = ( k1 = 6.8e-6, km1 = 4.5e-4, k2 = 7.7e-5, A = 1.13e-6, b = 0.5),
      tspan = (0.,350e3)
    ), 
    
  )

end

function objective_function_all(x,odes,systems,sim_col)
  f = 0.
  for i in 1:length(systems)
    ode_problem = odes[i]
    f += try 
      sim = solve(ode_problem,p=x,verbose=false)
      score(sim,systems[i].experimental_data,sim_col)
    catch
      return +Inf
    end
  end
  return f
end

function run_all(;max_failed_trials=10000,showprogress=true)

  systems = datasets("/home/leandro/.julia/dev/Polyuretanes/data")

#  reaction_network = reaction1()
  reaction_network = reaction2()

  odes = []
  for system in systems
    c0 = [values(system.c0)...]
    x = Float64[]
    for par in params(reaction_network)
      push!(x,getfield(system.p0,Symbol(par)))
    end
    tspan = system.tspan
    ode_problem = ODEProblem(reaction_network,c0,tspan,x)
    push!(odes,ode_problem)
  end
  x0 = Float64[]
  for par in params(reaction_network)
    push!(x0,getfield(systems[1].p0,Symbol(par)))
  end
  @. x0 = x0 + -x0 + rand()*x0

	f(x) = objective_function_all(x,odes,systems,4)
  fbest, xbest = optimizer(x0,f,max_failed_trials=max_failed_trials,showprogress=showprogress)

  r = Result[]
  for system in systems
    c0 = [values(system.c0)...]
    tspan = system.tspan
    sim = simulate(xbest,c0,tspan,reaction_network)
    sim_score = score(sim,system.experimental_data,4)
    push!(r,Result(system, sim, sim_score,
      ntuple(i -> Symbol(params(reaction_network)[i]) => xbest[i], length(xbest))))
  end

  return r, fbest
end

end # module