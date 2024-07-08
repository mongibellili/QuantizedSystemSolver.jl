function updateQ(::Val{1}, i::Int, xv::Vector{Taylor0}, qv::Vector{Taylor0}, quantum::Vector{Float64}, exactA::Function, d::Vector{Float64}, cacheA::MVector{1,Float64}, dxaux::Vector{MVector{1,Float64}}, qaux::Vector{MVector{1,Float64}}, tx::Vector{Float64}, tq::Vector{Float64}, simt::Float64, ft::Float64, nextStateTime::Vector{Float64})
    cacheA[1] = 0.0
    exactA(qv, d, cacheA, i, i, simt)
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
        if x1 > 0.0
            q = x + Δ
        else
            q = x - Δ
        end
        if x1 != 0
            h = (abs(Δ / x1))
        else
            h = Inf
        end
    end
    qv[i][0] = q
    nextStateTime[i] = simt + h
    return h
end
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
