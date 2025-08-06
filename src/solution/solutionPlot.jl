
@recipe function f(sol::LightSol; idxs=nothing)
     xlabel := "Time"
     ylabel := "Value"
     legend := :topright
     seriestype := :line

     if idxs === nothing
         # Plot all variables (each with its own time)
         for i in eachindex(sol.savedVars)
             x = sol.savedTimes[i]
             y = sol.savedVars[i]
             @series begin
                 label := "Var $i"
                 x, y
             end
         end
 
     elseif isa(idxs, AbstractVector)
         # Plot selected variables
         for i in idxs
             x = sol.savedTimes[i]
             y = sol.savedVars[i]
             @series begin
                 label := "Var $i"
                 x, y
             end
         end
 
     elseif isa(idxs, Tuple) && length(idxs) == 2
         i, j = idxs
 
         # Interpolate both to common time grid
         sol1 = solInterpolated(sol, i, sol.ft / 1000.0)
         sol2 = solInterpolated(sol, j, sol.ft / 1000.0)
 
         x = sol1.savedVars[1]
         y = sol2.savedVars[1]
 
         xlabel := "Var $i"
         ylabel := "Var $j"
 
         @series begin
             label := "Var $j vs Var $i"
             x, y
         end
 
        elseif isa(idxs, Tuple) && length(idxs) == 3
          i, j, k = idxs
          sol1 = solInterpolated(sol, i, sol.ft / 1000.0)
          sol2 = solInterpolated(sol, j, sol.ft / 1000.0)
          sol3 = solInterpolated(sol, k, sol.ft / 1000.0)
          x = sol1.savedVars[1]
          y = sol2.savedVars[1]
          z = sol3.savedVars[1]
          xlabel := "Var $i"
          ylabel := "Var $j"
          zlabel := "Var $k"
          @series begin
              label := "3D Trajectory"
              x, y, z
          end
  
      else
          error("Invalid idxs argument. Must be nothing, a vector, a tuple of length 2 (for 2D), or length 3 (for 3D).")
      end
 end
 
 




"""
    plot_SolSum(sol::Sol{T,O},xvars::Int...;interp=0.0001,note=" "::String,xlims::Tuple{Float64, Float64}=(0.0,0.0),ylims::Tuple{Float64, Float64}=(0.0,0.0),legend::Bool=true) where{T,O}

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
function plot_SolSum(sol::Sol{T,O},xvars::Int...;interp=0.0001,note=" "::String,xlims::Tuple{Float64, Float64}=(0.0,0.0),ylims::Tuple{Float64, Float64}=(0.0,0.0),legend::Bool=true) where{T,O}
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
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.stats.totalSteps)_$(sol.stats.simulSteps)_$(sol.stats.eventSteps) \n $note", xlims=xlims ,ylims=ylims,legend=legend)
  elseif xlims!=(0.0,0.0) && ylims==(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.stats.totalSteps)_$(sol.stats.simulSteps)_$(sol.stats.eventSteps) \n $note", xlims=xlims #= ,ylims=ylims =#,legend=legend)
  elseif xlims==(0.0,0.0) && ylims!=(0.0,0.0) 
    p1=plot!(p1,title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.stats.totalSteps)_$(sol.stats.simulSteps)_$(sol.stats.eventSteps) \n $note"#= , xlims=xlims  =#,ylims=ylims,legend=legend)
  else
    p1=plot!(p1, title="$(sol.sysName)_$(sol.algName)_$(sol.absQ)_$(sol.stats.totalSteps)_$(sol.stats.simulSteps)_$(sol.stats.eventSteps) \n $note",legend=legend)
  end
  p1
end




