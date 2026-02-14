
"""
    updateQ(::Val{2}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}, a::Float64, dxaux::Vector{MVector{2,Float64}}, qaux::Vector{MVector{2,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64})

Update the quantized state for the second-order quantizer.

# Arguments
- `::Val{2}`: Type parameter indicating the second-order quantizer.
- `i::Int`: Index of the state variable to update.
- `xv::Vector{Taylor0}`: Vector of state variables.
- `qv::Vector{Taylor0}`: Vector of quantized state variables.
- `quantum::Vector{Float64}`: Vector of quantum values of the state variables.
- `exactA::Function`: Function to compute the exact value of a jacobian entry.
- `d::Vector{Float64}`: Vector of discrete variables.
- `cacheA::MVector{1,Float64}`: Cache for jacobian entry computation.
- `dxaux::Vector{MVector{2,Float64}}`: Auxiliary vector for saving old x values.
- `qaux::Vector{MVector{2,Float64}}`: Auxiliary vector for saving old quantized values.
- `tx::Vector{Float64}`: Vector of state update times.
- `tq::Vector{Float64}`: Vector of quantized state update times.
- `simt::Float64`: Current simulation time.
- `ft::Float64`: Final time of the simulation.
- `nextStateTime::Vector{Float64}`: Vector of times for the next state updates.


"""
function updateQ(::Val{2}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}, a::Float64, dxaux::Vector{MVector{2,Float64}}, qaux::Vector{MVector{2,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64})
    q = qv[i][0]
    q1 = qv[i][1]
    x = xv[i][0]
    x1 = xv[i][1]
    x2 = xv[i][2] * 2 #u1=uv[i][i][1]; u2=uv[i][i][2]
    qaux[i][1] = q + (simt - tq[i]) * q1#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2] = q1                     #appears only here...updated here and used in updateQevent
    u1 = x1 - a * qaux[i][1]
    u2 = x2 - a * q1
    dxaux[i][1] = x1
    dxaux[i][2] = x2
    ddx = x2
    quan = quantum[i]
    h = 0.0

 

    if isnan(a)
        @warn("a is NaN: The Jacobian is not defined at this instant $(simt). This may be due to an undefined operation. Consider computing the Jacobian coefficient manually.")
        a = 0.0
    end
    if a != 0.0
        if a * a * x + a * u1 + u2 <= 0.0
            if -(a * a * x + a * u1 + u2) / (a * a) < quan # asymptote<delta...no sol ...no need to check
                h = Inf
                q = -(a * u1 + u2) / (a * a)#q=x+asymp
            else
                q = x + quan
                h = minPosRoot(quan,-a * quan,(a * a * (x + quan) + a * u1 + u2) / 2,Val(2)) #needed for dq
                #math was used to show h cannot be 1/a
            end
        elseif a * a * x + a * u1 + u2 >= 0.0
            if -(a * a * x + a * u1 + u2) / (a * a) > -quan # asymptote>-delta...no sol ...no need to check
                h = Inf
                q = -(a * u1 + u2) / (a * a)
            else
                h = minPosRoot( -quan,a * quan,(a * a * (x - quan) + a * u1 + u2) / 2,Val(2))
                q = x - quan
            end
        else#a*a*x+a*u1+u2==0 -->f=0....q=x+f(h)=x+h*g(h) -->g(h)==0...dx==0 -->aq+u==0
            q = -u1 / a
            h = Inf
        end
        if h != Inf
            if h != 1 / a
                q1 = (a * q + u1 + h * u2) / (1 - h * a)
            else#h=1/a
                error("report bug: updateQ: h cannot=1/a; h=$h , 1/a=$(1/h)")
            end
        else #h==inf make ddx==0 dq=-u2/a
            q1 = -u2 / a
        end
    else
        #println("a==0")
        if x2 != 0.0
            h = sqrt(abs(2 * quan / x2))   #sqrt necessary with u2
            q = x - h * h * x2 / 2
            q1 = x1 + h * x2
        else
            # println("x2==0")
            if x1 != 0.0
                #quantum[i]=1quan
                h = abs(1 * quan / x1)   # *10 just to widen the step otherwise it would behave like 1st order
                q = x + h * x1
                q1 = 0.0
            else
                h = Inf
                q = x
                q1 = 0.0
            end
        end
    end

 
    qv[i][0] = q
    qv[i][1] = q1
    nextStateTime[i] = simt + h
    #return h
end

"""
    updateQInit(::Val{2}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}, a::Float64,dxaux::Vector{MVector{2,Float64}}, qaux::Vector{MVector{2,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64}) 

Initialize the quantized state variables for the LIQSS2 method. It is similar to the [`updateQ`](@ref) function but does not accept q to be set to x when all derivatives are zero, which is the case when an equilibrium ocurrs during the simulation.


"""
function updateQInit(::Val{2}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}, a::Float64,dxaux::Vector{MVector{2,Float64}}, qaux::Vector{MVector{2,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64}) 
 
    q = qv[i][0]
    q1 = qv[i][1]
    x = xv[i][0]
    x1 = xv[i][1]
    x2 = xv[i][2] * 2 #u1=uv[i][i][1]; u2=uv[i][i][2]
    qaux[i][1] = q + (simt - tq[i]) * q1
    qaux[i][2] = q1                    
    u1 = x1 - a * qaux[i][1]
    u2 = x2 - a * q1
    dxaux[i][1] = x1
    dxaux[i][2] = x2
    ddx = x2
    quan = quantum[i]
    h = 0.0
    if isnan(a)
        @warn("a is NaN: The Jacobian is not defined at this instant $(simt). This may be due to  an undefined operation. Consider computing the Jacobian coefficient manually.")
        a = 0.0
    end
    if a != 0.0
        if a * a * x + a * u1 + u2 <= 0.0
            if -(a * a * x + a * u1 + u2) / (a * a) < quan # asymptote<delta...no sol ...no need to check
                h = Inf
                q = -(a * u1 + u2) / (a * a)#q=x+asymp
            else
                q = x + quan
                h = minPosRoot(quan,-a * quan,(a * a * (x + quan) + a * u1 + u2) / 2,Val(2))
                #math was used to show h cannot be 1/a
            end
        elseif a * a * x + a * u1 + u2 >= 0.0
            if -(a * a * x + a * u1 + u2) / (a * a) > -quan # asymptote>-delta...no sol ...no need to check
                h = Inf
                q = -(a * u1 + u2) / (a * a)
            else
                h = minPosRoot( -quan,a * quan,(a * a * (x - quan) + a * u1 + u2) / 2,Val(2))
                q = x - quan
            end
        else#a*a*x+a*u1+u2==0 -->f=0....q=x+f(h)=x+h*g(h) -->g(h)==0...dx==0 -->aq+u==0
            q = -u1 / a
            h = Inf
        end
        if h != Inf
            if h != 1 / a
                q1 = (a * q + u1 + h * u2) / (1 - h * a)
            else#h=1/a
                error("report bug: updateQ: h cannot=1/a; h=$h , 1/a=$(1/h)")
            end
        else #h==inf make ddx==0 dq=-u2/a
            q1 = -u2 / a
        end
    else
        if x2 != 0.0
            h = sqrt(abs(2 * quan / x2))   #sqrt necessary with u2
            q = x - h * h * x2 / 2
            q1 = x1 + h * x2
        else
            if x1 != 0.0
                h = abs(1 * quan / x1)   # *10 just to widen the step otherwise it would behave like 1st order
                q = x + h * x1
                q1 = 0.0
            else
                h = Inf
                q = x+quan
                q1 = 0.0
            end
        end
    end
    qv[i][0] = q
    qv[i][1] = q1
    nextStateTime[i] = simt + h
    return h
end

"""
    Liqss_reComputeNextTime(::Val{2}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64})

Recomputes the next time for a given state in a second-order quantized state system.

# Arguments
- `::Val{2}`: A type parameter indicating the order of the quantized state system (second-order in this case).
- `i::Int`: The index of the state for which the next time is being recomputed.
- `simt::Float64`: The current simulation time.
- `nextStateTime::Vector{Float64}`: A vector containing the next state times for all states.
- `xv::Vector{Taylor0}`: A vector containing the current state values represented as Taylor series.
- `qv::Vector{Taylor0}`: A vector containing the quantized state values represented as Taylor series.
- `quantum::Vector{Float64}`: A vector containing the quantum values the states.

# Returns
- This function does not return a value. It updates the `nextStateTime` vector in place.
"""
function Liqss_reComputeNextTime(::Val{2}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}) 
    q = qv[i][0]
    x = xv[i][0]
    q1 = qv[i][1]
    x1 = xv[i][1]
    x2 = xv[i][2]
    quani = quantum[i]
    β = 0
    if abs(q - x) >= 2 * quani # this happened when var i and j s turns are now...var i depends on j, j is asked here for next time...or if you want to increase quant*10 later it can be put back to normal and q & x are spread out by 10quan
        nextStateTime[i] = simt + 1e-12
    else
        nextStateTime[i] = simt + minPosRoot(q - x, q1 - x1, -x2, Val(2))
        if q - x > 0.0#1e-9
            timetemp = simt + minPosRoot(q - x - 2 * quantum[i], q1 - x1, -x2, Val(2))
            if timetemp < nextStateTime[i]
                nextStateTime[i] = timetemp
            end
        elseif q - x < 0.0#-1e-9
            timetemp = simt + minPosRoot(q - x + 2 * quantum[i], q1 - x1, -x2, Val(2))
            if timetemp < nextStateTime[i]
                nextStateTime[i] = timetemp
            end
        end
        if nextStateTime[i] <= simt # this is coming from the fact that a variable can reach 2quan distance when it is not its turn, then computation above gives next=simt+(p-p)/dx...p-p should be zero but it can be very small negative
            nextStateTime[i] = simt + Inf#1e-14
        end
    end
end
