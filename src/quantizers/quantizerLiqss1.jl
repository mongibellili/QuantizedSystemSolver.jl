"""
    updateQ(::Val{1}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}, exactA::Function, d::Vector{Float64}, cacheA::MVector{1,Float64}, dxaux::Vector{MVector{1,Float64}}, qaux::Vector{MVector{1,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64},f::F) where{F} 

Update the quantized state for the LIQSS1 (Linearly Implicit Quantized State System 1) method.

# Arguments
- `::Val{1}`: Type parameter indicating the LIQSS1 method.
- `i::Int`: Index of the state variable to update.
- `xv::Vector{Taylor0}`: Vector of current state values.
- `qv::Vector{Taylor0}`: Vector of quantized state values.
- `quantum::Vector{Float64}`: Vector of quantum values for the state variables.
- `exactA::Function`: Function to compute the exact value of a jacobian entry.
- `d::Vector{Float64}`: Vector of discrete variables.
- `cacheA::MVector{1,Float64}`: Cache for jacobian entry computation.
- `dxaux::Vector{MVector{1,Float64}}`: Auxiliary vector for saving old x values.
- `qaux::Vector{MVector{1,Float64}}`: Auxiliary vector for saving old quantized values.
- `tx::Vector{Float64}`: Vector of times at which the state variables were updated.
- `tq::Vector{Float64}`: Vector of times at which the quantized state variables were updated.
- `simt::Float64`: Current simulation time.
- `ft::Float64`: Final time of the simulation.
- `nextStateTime::Vector{Float64}`: Vector of times at which the state variables will be updated next.

# Returns
- None. The function updates the quantized state and other info in place.
"""
function updateQ(::Val{1}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}, exactA::Function, d::Vector{Float64}, cacheA::MVector{1,Float64}, dxaux::Vector{MVector{1,Float64}}, qaux::Vector{MVector{1,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64},f::F) where{F} 
    cacheA[1] = 0.0
    exactA(qv, d, cacheA, i, i, simt,f)
    a = cacheA[1]
    q = qv[i][0]
    x = xv[i][0]
    x1 = xv[i][1]
    qaux[i][1] = q
    u = x1 - a * q
    dxaux[i][1] = x1
    h = 0.0
    Δ = quantum[i]
    if a != 0.0
        α = -(a * x + u) / a
        h1denom = a * (x + Δ) + u
        h2denom = a * (x - Δ) + u
        h1 = Δ / h1denom
        h2 = -Δ / h2denom
        if a < 0
            if α > Δ
                h = h1
                q = x + Δ
            elseif α < -Δ
                h = h2
                q = x - Δ
            else
                h = Inf
                q = -u / a
            end
        else
            if α > Δ
                h = h2
                q = x - Δ
            elseif α < -Δ
                h = h1
                q = x + Δ
            else
                if a * x + u > 0
                    q = x - Δ
                    h = h2
                else
                    q = x + Δ
                    h = h1
                end
            end
        end
    else
        #= if x1 > 0.0
            q = x + Δ
        else
            q = x - Δ
        end
        if x1 != 0
            h = (abs(Δ / x1))
        else
            h = Inf
        end =#
        if x1 != 0.0
            #quantum[i]=1quan
            h = abs(1 * Δ / x1)   # *10 just to widen the step otherwise it would behave like 1st order
            q = x + h * x1
        else
            h = Inf
            q = x
        end
    end
    qv[i][0] = q
    nextStateTime[i] = simt + h
    return h
end
"""
    updateQInit(::Val{1}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}, exactA::Function, d::Vector{Float64}, cacheA::MVector{1,Float64}, dxaux::Vector{MVector{1,Float64}}, qaux::Vector{MVector{1,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64},f::F) where{F} 

Initialize the quantized state for the LIQSS1 method.

# Arguments
- `::Val{1}`: Type parameter indicating the LIQSS1 method.
- `i::Int`: Index of the state variable to update.
- `xv::Vector{Taylor0}`: Vector of state variables.
- `qv::Vector{Taylor0}`: Vector of quantized state variables.
- `quantum::Vector{Float64}`: Vector of quantum values for the state variables.
- `exactA::Function`: Function to compute the exact value of the state variable.
- `d::Vector{Float64}`: Vector of derivatives of the state variables.
- `cacheA::MVector{1,Float64}`: Cache for intermediate computations.
- `dxaux::Vector{MVector{1,Float64}}`: Auxiliary vector for derivatives.
- `qaux::Vector{MVector{1,Float64}}`: Auxiliary vector for quantized states.
- `tx::Vector{Float64}`: Vector of times at which state variables were last updated.
- `tq::Vector{Float64}`: Vector of times at which quantized state variables were last updated.
- `simt::Float64`: Current simulation time.
- `ft::Float64`: Final simulation time.
- `nextStateTime::Vector{Float64}`: Vector of times at which the next state update is scheduled.

# Description
This function initializes the quantized state for the LIQSS1 method by updating the quantized state variables and their associated times based on the provided state variables, derivatives, and quantum values.
"""
function updateQInit(::Val{1}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}, exactA::Function, d::Vector{Float64}, cacheA::MVector{1,Float64}, dxaux::Vector{MVector{1,Float64}}, qaux::Vector{MVector{1,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64},f::F) where{F} 
    cacheA[1] = 0.0
    exactA(qv, d, cacheA, i, i, simt,f)
    a = cacheA[1]
    q = qv[i][0]
    x = xv[i][0]
    x1 = xv[i][1]
    qaux[i][1] = q
    u = x1 - a * q
    dxaux[i][1] = x1
    h = 0.0
    Δ = quantum[i]
    if a != 0.0
        α = -(a * x + u) / a
        h1denom = a * (x + Δ) + u
        h2denom = a * (x - Δ) + u
        h1 = Δ / h1denom
        h2 = -Δ / h2denom
        if a < 0
            if α > Δ
                h = h1
                q = x + Δ
            elseif α < -Δ
                h = h2
                q = x - Δ
            else
                h = Inf
                q = -u / a
            end
        else
            if α > Δ
                h = h2
                q = x - Δ
            elseif α < -Δ
                h = h1
                q = x + Δ
            else
                if a * x + u > 0
                    q = x - Δ
                    h = h2
                else
                    q = x + Δ
                    h = h1
                end
            end
        end
    else
      #=   if x1 > 0.0
            q = x + Δ
        else
            q = x - Δ
        end
        if x1 != 0
            h = (abs(Δ / x1))
        else
            h = Inf 
        end =#
        if x1 != 0.0
            #quantum[i]=1quan
            h = abs(1 * Δ / x1)   # *10 just to widen the step otherwise it would behave like 1st order
            q = x + h * x1
        else
            h = Inf
            q = x+Δ
        end
    end
    qv[i][0] = q
    nextStateTime[i] = simt + h
    return h
end


"""
    Liqss_reComputeNextTime(::Val{1}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64})

Recomputes the next time for the LIQSS1 quantizer.

# Arguments
- `::Val{1}`: Type parameter indicating the LIQSS1 method.
- `i::Int`: Index of the state variable.
- `simt::Float64`: Current simulation time.
- `nextStateTime::Vector{Float64}`: Vector containing the next state times for each state variable.
- `xv::Vector{Taylor0}`: Vector of current state values represented as Taylor series.
- `qv::Vector{Taylor0}`: Vector of quantized state values represented as Taylor series.
- `quantum::Vector{Float64}`: Vector of quantum values for the state variables.

# Returns
- Updates the `nextStateTime` vector with the recomputed next time for the specified state variable.
"""
function Liqss_reComputeNextTime(::Val{1}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64})
    dt = 0.0
    q = qv[i][0]
    x = xv[i][0]
    x1 = xv[i][1]
    if abs(q - x) >= 2 * quantum[i]
        nextStateTime[i] = simt + 1e-12
    else
        if x1 != 0.0
            dt = (q - x) / x1
            if dt > 0.0
                nextStateTime[i] = simt + dt
            elseif dt < 0.0
                if x1 > 0.0
                    nextStateTime[i] = simt + (q - x + 2 * quantum[i]) / x1
                else
                    nextStateTime[i] = simt + (q - x - 2 * quantum[i]) / x1
                end
            end
        else
            nextStateTime[i] = Inf
        end
    end
    if nextStateTime[i] <= simt
        nextStateTime[i] = simt + Inf
    end
end
