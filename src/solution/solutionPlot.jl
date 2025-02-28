"""
    plot(sol::Sol{T,O};idxs=Int[]::Vector{Int},note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool,marker=:circle::Symbol,title="") where{T,O}

generates a plot of the solution of the system (returned as a plot object).
    With the exception of the solution object, all arguments are optional.
  
# arguments:
- `sol::Sol{T,O}`: The solution struct.
- `xvars::Int...`: The indices of the variables to plot. If no indices are provided, all variables are plotted.
- `idxs::::Vector{Float64}`: The indices of the variables to plot. If no indices are provided, all variables are plotted.
- `note::String`: A note to add to the title of the plot.
- `xlims::Tuple{Float64, Float64}`: The x-axis limits of the plot.
- `ylims::Tuple{Float64, Float64}`: The y-axis limits of the plot.
- `legend::Bool`: A boolean indicating whether to display the legend.
- `marker::Symbol`: The marker to use for the plot.
- `title::String`: The title of the plot.
  
"""
function plot(sol::Sol{T,O};idxs=Int[]::Vector{Int},note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool,marker=:circle::Symbol,title="") where{T,O}
  if title==""
    title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.stats.totalSteps)_$(sol.stats.simulStepCount)_$(sol.stats.evCount) \n $note"
  end

  p1=plot()

  if idxs isa Tuple
      if length(idxs)==2
        if idxs[1]==0
          k=idxs[2]
          p1=plot!(sol.savedTimes[k], sol.savedVars[k],marker=(marker),label="x$k $(sol.stats.numSteps[k])")
        else
          sol1=solInterpolated(sol,idxs[1],sol.ft/100.0)
          sol2=solInterpolated(sol,idxs[2],sol.ft/100.0)
          p1=plot(sol1[2][1], sol2[2][1],marker=(marker),markersize=2,xlabel="u$(idxs[1])",ylabel="u$(idxs[2])")
          legend=false
        end
      elseif length(idxs)==3
        if idxs[1]==0
          
          sol2=solInterpolated(sol,idxs[2],sol.ft/100.0)
          sol3=solInterpolated(sol,idxs[3],sol.ft/100.0)
          p1=plot!(sol2[1][1], sol2[2][1],sol3[2][1],marker=(marker),markersize=2)
        else
          sol1=solInterpolated(sol,idxs[1],sol.ft/100.0)
          sol2=solInterpolated(sol,idxs[2],sol.ft/100.0)
          sol3=solInterpolated(sol,idxs[3],sol.ft/100.0)
          p1=plot!(sol1[2][1], sol2[2][1],sol3[2][1],marker=(marker),markersize=2)
        end
      else
        error("idxs=() takes 2 or 3 elements")
      end
  elseif idxs isa Vector
      if length(idxs)!=0
        for k in idxs 
          if k==1
            stle=:dash
            sze=2
          elseif k==2
            stle=:solid
            sze=4
          elseif k==3
            stle=:dot
            sze=3
          elseif k==4
            stle=:dashdot
            sze=2
          else
            sze=1
            stle=:solid
          end
          p1=plot!(p1,sol.savedTimes[k], sol.savedVars[k],line=(sze,stle),marker=(marker),label="x$k $(sol.stats.numSteps[k])")
        end
      else
        for k=1:T
          p1=plot!(p1,sol.savedTimes[k], sol.savedVars[k],marker=(marker),markersize=2,label="x$k $(sol.stats.numSteps[k])")
        end
      end
  else
      error("idxs must be a Tuple or a Vector{Int64}")
  end 
  if xlims!=(0.0,0.0) && ylims!=(0.0,0.0) 
    p1=plot!(p1,title=title, xlims=xlims ,ylims=ylims,legend=legend)
  elseif xlims!=(0.0,0.0) && ylims==(0.0,0.0) 
    p1=plot!(p1,title=title, xlims=xlims,legend=legend)
  elseif xlims==(0.0,0.0) && ylims!=(0.0,0.0) 
    p1=plot!(p1,title=title,ylims=ylims,legend=legend)
  else
    p1=plot!(p1, title=title,legend=legend)#user did not enter this option
  end
  p1
end




"""
    plot_Sol(sol::Sol{T,O},xvars::Int...;note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool,marker=:circle::Symbol,title="") where{T,O}

generates a plot of the solution of the system (returned as a plot object).
    With the exception of the solution object, all arguments are optional.
  
# arguments:
- `sol::Sol{T,O}`: The solution struct.
- `xvars::Int...`: The indices of the variables to plot. If no indices are provided, all variables are plotted.
- `note::String`: A note to add to the title of the plot.
- `xlims::Tuple{Float64, Float64}`: The x-axis limits of the plot.
- `ylims::Tuple{Float64, Float64}`: The y-axis limits of the plot.
- `legend::Bool`: A boolean indicating whether to display the legend.
- `marker::Symbol`: The marker to use for the plot.
- `title::String`: The title of the plot.

"""
function plot_Sol(sol::Sol{T,O},xvars::Int...;note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool,marker=:circle::Symbol,title="") where{T,O}
  idxs=Int64[]
  if xvars!=()
    for k in xvars 
      push!(idxs,k)
    end
  end
  plot(sol,idxs=idxs,note=note,xlims=xlims,ylims=ylims,legend=legend,marker=marker,title=title)
end



"""
    plot_SolSum(sol::Sol{T,O},xvars::Int...;interp=0.0001,note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool) where{T,O}

plots of the sum of the variables xvars.

# Arguments 
- `sol::Sol{T,O}`: The solution struct.
- `xvars::Int...`: The indices of the variables to sum. If no indices are provided, all variables are summed.
- `interp::Float64`: The interpolation step.
- `note::String`: A note to add to the title of the plot.
- `xlims::Tuple{Float64, Float64}`: The x-axis limits of the plot.
- `ylims::Tuple{Float64, Float64}`: The y-axis limits of the plot.
- `legend::Bool`: A boolean indicating whether to display the legend.


"""  
function plot_SolSum(sol::Sol{T,O},xvars::Int...;interp=0.0001,note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool) where{T,O}
  p1=plot()
  if xvars!=()
    solInterp=solInterpolated(sol,interp) # for now interpl all
    sumV=solInterp.savedVars[xvars[1]]
    sumT=solInterp.savedTimes[xvars[1]]# 
    numPoints=length(sumT)
    for i=1:numPoints
        for k=2: length(xvars) 
          sumV[i]+=solInterp.savedVars[xvars[k]][i]
        end
    end
      p1=plot!(p1,sumT, sumV,marker=(:circle),#= ,legend=:right =#)
  else
    println("pick idxs to plot their sum")
  end
  if xlims!=(0.0,0.0) && ylims!=(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.stats.totalSteps)_$(sol.stats.simulStepCount)_$(sol.stats.evCount) \n $note", xlims=xlims ,ylims=ylims,legend=legend)
  elseif xlims!=(0.0,0.0) && ylims==(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.stats.totalSteps)_$(sol.stats.simulStepCount)_$(sol.stats.evCount) \n $note", xlims=xlims #= ,ylims=ylims =#,legend=legend)
  elseif xlims==(0.0,0.0) && ylims!=(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.stats.totalSteps)_$(sol.stats.simulStepCount)_$(sol.stats.evCount) \n $note"#= , xlims=xlims  =#,ylims=ylims,legend=legend)
  else
    p1=plot!(p1, title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.stats.totalSteps)_$(sol.stats.simulStepCount)_$(sol.stats.evCount) \n $note",legend=legend)
  end
  p1
end



#= """
    save_Sol(sol::Sol{T,O},xvars::Int...;note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool) where{T,O}

Saves the plot of the system solution for the variables xvars.
  With the exception of the solution object, all arguments are optional.
  The default values are:\n
  - note = " "\n
  - xlims = (0.0,0.0)\n
  - ylims = (0.0,0.0)\n
  - legend = true

"""
function save_Sol(sol::Sol{T,O},xvars::Int...;note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool) where{T,O}
  p1= plot_Sol(sol,xvars...;note=note,xlims=xlims,ylims=ylims,legend=legend)
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(xvars)_$(sol.absQ)_$(note)_ft_$(sol.ft)_$(timestamp).png")
end
"""
    save_SolSum(sol::Sol{T,O},xvars::Int...;interp=0.0001,note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool) where{T,O}

Saves the plot of the sum of the variables xvars.

# Arguments 
- `sol::Sol{T,O}`: The solution struct.
- `xvars::Int...`: The indices of the variables to sum. If no indices are provided, all variables are summed.
- `interp::Float64`: The interpolation step.
- `note::String`: A note to add to the title of the plot.
- `xlims::Tuple{Float64, Float64}`: The x-axis limits of the plot.
- `ylims::Tuple{Float64, Float64}`: The y-axis limits of the plot.
- `legend::Bool`: A boolean indicating whether to display the legend.

"""  
function save_SolSum(sol::Sol{T,O},xvars::Int...;interp=0.0001,note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool) where{T,O}
  p1= plot_SolSum(sol,xvars...;interp=interp,note=note,xlims=xlims,ylims=ylims,legend=legend)
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(xvars)_$(sol.absQ)_$(note)_ft_$(sol.ft)_$(timestamp).png")
end
 =#


