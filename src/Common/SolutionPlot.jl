"""plot_Sol(sol::Sol{T,O},xvars::Int...;note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool) where{T,O}

  This function generates a plot of the solution of the system (returned as a plot object).
"""
function plot_Sol(sol::Sol{T,O},xvars::Int...;note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool) where{T,O}
  p1=plot()
  if xvars!=()
    for k in xvars
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
      p1=plot!(p1,sol.savedTimes[k], sol.savedVars[k],line=(sze,stle),marker=(:circle),label="x$k $(sol.numSteps[k])"#= ,legend=:right =#)
    end
  else
    for k=1:T
      if k==1
        mycolor=:red
      else
        mycolor=:purple
      end
      p1=plot!(p1,sol.savedTimes[k], sol.savedVars[k],marker=(:circle),markersize=2,label="x$k $(sol.numSteps[k])"#= ,legend=:false =#)
    end
  end
  if xlims!=(0.0,0.0) && ylims!=(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount)_$(sol.evCount) \n $note", xlims=xlims ,ylims=ylims)
  elseif xlims!=(0.0,0.0) && ylims==(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount)_$(sol.evCount)  \n $note", xlims=xlims #= ,ylims=ylims =#)
  elseif xlims==(0.0,0.0) && ylims!=(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount)_$(sol.evCount)  \n $note"#= , xlims=xlims  =#,ylims=ylims)
  else
    p1=plot!(p1, title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount)_$(sol.evCount)  \n $note",legend=legend)
  end
  p1
end
"""save_Sol(sol::Sol{T,O},xvars::Int...;note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool) where{T,O}

Save the plot of the system solution for the variables xvars.
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
function save_SolSum(sol::Sol{T,O},xvars::Int...;interp=0.0001,note=" "::String,xlims=(0.0,0.0)::Tuple{Float64, Float64},ylims=(0.0,0.0)::Tuple{Float64, Float64},legend=:true::Bool) where{T,O}
  p1= plot_SolSum(sol,xvars...;interp=interp,note=note,xlims=xlims,ylims=ylims,legend=legend)
  mydate=now()
  timestamp=(string(year(mydate),"_",month(mydate),"_",day(mydate),"_",hour(mydate),"_",minute(mydate),"_",second(mydate)))
  savefig(p1, "plot_$(sol.sysName)_$(sol.algName)_$(xvars)_$(sol.absQ)_$(note)_ft_$(sol.ft)_$(timestamp).png")
end


#plot the sum of variables
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
    println("pick vars to plot their sum")
  end
  if xlims!=(0.0,0.0) && ylims!=(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount)_$(sol.evCount) \n $note", xlims=xlims ,ylims=ylims)
  elseif xlims!=(0.0,0.0) && ylims==(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount)_$(sol.evCount)  \n $note", xlims=xlims #= ,ylims=ylims =#)
  elseif xlims==(0.0,0.0) && ylims!=(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount)_$(sol.evCount)  \n $note"#= , xlims=xlims  =#,ylims=ylims)
  else
    p1=plot!(p1, title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.totalSteps)_$(sol.simulStepCount)_$(sol.evCount)  \n $note",legend=legend)
  end
  p1
end

