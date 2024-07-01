
function updateQ(::Val{3},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64},exacteA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{3,Float64}},qaux::Vector{MVector{3,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})

    #a=av[i][i]
    cacheA[1]=0.0;exacteA(qv,d,cacheA,i,i);a=cacheA[1]
    q=qv[i][0];q1=qv[i][1];q2=2*qv[i][2];x=xv[i][0];x1=xv[i][1];x2=2*xv[i][2];x3=6*xv[i][3];u1=uv[i][i][1];u2=uv[i][i][2];u3=uv[i][i][3]
    elapsed=simt-tq[i]
    qaux[i][1]=q+elapsed*q1+elapsed*elapsed*q2/2#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1+elapsed*q2   ;qaux[i][3]=q2     #never used
    olddx[i][1]=x1  
   # tq[i]=simt
   # elapsed=simt-tu[i]
   # u1=u1+elapsed*u2+elapsed*elapsed*u3/2  
    u1=x1-a*qaux[i][1]
    uv[i][i][1]=u1
   #=  u2=u2+elapsed*u3 
    uv[i][i][2]=u2 =#
   uv[i][i][2]=x2-a*qaux[i][2]  #---------------------------------------------------------------
    u2=uv[i][i][2]
   uv[i][i][3]=x3-a*q2
   u3=uv[i][i][3]
   # tu[i]=simt
    dddx=x3
    dxaux[i][1]=x1
    dxaux[i][2]=x2
    dxaux[i][3]=x3

    quan=quantum[i]
    h=0.0
   # println("before q update",abs(q - x) > 2 * quan)
     if a!=0.0
        if dddx ==0.0
            dddx=a*a*a*(q)+a*a*u1+a*u2+u3 #*2
            if dddx==0.0
                dddx=1e-40# changing -40 to -6 nothing changed
                println("dddx=0")  #this appeared once with sys1 liqss3
            end
        end
     

        h = ft-simt
        α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
        λ=x+h*u1+h*h*u2/2+h*h*h*u3/6
        β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+λ
        γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)
        q = β/γ
        if (abs(q - x) >  2*quan) # removing this did nothing...check @btime later
          h = cbrt(abs((6*quan) / dddx));
          #h= cbrt(abs((q-x) / x3));#h=cbrt(abs(6*(q-x) / x3))# shifts up a little
          α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
          β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+x+h*u1+h*h*u2/2+h*h*h*u3/6
        γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)
        q = β/γ
        end

        
        maxIter=515
        while (abs(q - x) >  2*quan) && (maxIter>0)
            maxIter-=1
         # h = h *(0.98*quan / abs(q - x));
          h = h *cbrt(quan / abs(q - x));
          α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
          β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+x+h*u1+h*h*u2/2+h*h*h*u3/6
        γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)

        q = β/γ
        end
        q1=(a*(1-h*a)*q+u1*(1-h*a)-h*h*(a*u2+u3)/2)/(1-h*a+h*h*a*a/2)
        q2=(a*q1+u2+h*u3)/(1-h*a)

        if maxIter <200
            @show maxIter
        end
        if h==0.0
            @show h
        end
 
    else
       
        if x3!=0.0
            h=cbrt(abs(6*quan/x3))
            q=x+h*h*h*x3/6
            q1=x1-x3*h*h/2   #*2
            q2=x2+h*x3
       else
            if x2!=0
                h=sqrt(abs(2*quan/x2))
                q=x-h*h*x2/2
                q1=x1+h*x2  
            else
                if x1!=0
                    h=abs(quan/x1)
                    q=x+h*x1
                else
                    q=x
                    h=Inf
                    
                end
                q1=x1
                
            end
            q2=x2
        end 
    end
    qv[i][0]=q
    qv[i][1]=q1 
    qv[i][2]=q2/2  
    nextStateTime[i]=simt+h
    return nothing
end

function nupdateQ(::Val{3},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64},av::Vector{Vector{Float64}},uv::Vector{Vector{MVector{O,Float64}}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})where{O}
    a=av[i][i]
    q=qv[i][0];q1=qv[i][1];q2=2*qv[i][2];x=xv[i][0];x1=xv[i][1];x2=2*xv[i][2];x3=6*xv[i][3];u1=uv[i][i][1];u2=uv[i][i][2];u3=uv[i][i][3]
    elapsed=simt-tq[i]
    qaux[i][1]=q+elapsed*q1+elapsed*elapsed*q2/2#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1+elapsed*q2   ;qaux[i][3]=q2     #never used
    olddx[i][1]=x1  
   # tq[i]=simt
    #elapsed=simt-tu[i]
    #u1=u1+elapsed*u2+elapsed*elapsed*u3/2  
    u1=x1-a*qaux[i][1]
    uv[i][i][1]=u1
   # u2=u2+elapsed*u3 
   # uv[i][i][2]=u2
    uv[i][i][2]=x2-a*qaux[i][2]
    u2=uv[i][i][2]
   uv[i][i][3]=x3-a*q2
   u3=uv[i][i][3]
    #tu[i]=simt
    dddx=x3
   
    quan=quantum[i]
    h=0.0

     if a!=0.0
        if dddx ==0.0
            dddx=a*a*a*(q)+a*a*u1+a*u2+u3 #*2
            if dddx==0.0
                dddx=1e-40# changing -40 to -6 nothing changed
               println("nupdate dddx=0")  #this appeared once with sys1 liqss3
            end
        end
 
        h = ft-simt
        α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
        λ=x+h*u1+h*h*u2/2+h*h*h*u3/6
        β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+λ
        γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)
        q = β/γ
        if (abs(q - x) > 2* quan) # removing this did nothing...check @btime later
          h = cbrt(abs((6*quan) / dddx));
         
          #h= cbrt(abs((q-x) / x3));#h=cbrt(abs(6*(q-x) / x3))# shifts up a little
          α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
          β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+x+h*u1+h*h*u2/2+h*h*h*u3/6
        γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)
        q = β/γ
        end

        
        maxIter=515
        while (abs(q - x) > 2* quan) && (maxIter>0)
            maxIter-=1
         # h = h *(0.98*quan / abs(q - x));
          h = h *cbrt(quan / abs(q - x))
          α=h*(1-a*h+h*h*a*a/3)/(1-h*a)
          β=-α*(u1-u1*h*a-h*h*(a*u2+u3)/2)/(1-a*h+h*h*a*a/2)-h*h*(0.5-h*a/6)*(u2+h*u3)/(1-a*h)+x+h*u1+h*h*u2/2+h*h*h*u3/6
        γ=1-a*h+α*a*(1-a*h)/(1-a*h+h*h*a*a/2)

        q = β/γ
        end
        q1=(a*(1-h*a)*q+u1*(1-h*a)-h*h*(a*u2+u3)/2)/(1-h*a+h*h*a*a/2)
        q2=(a*q1+u2+h*u3)/(1-h*a)
        if maxIter <200
            @show maxIter
        end
        if h==0.0
            @show h
        end
 
    else
       
        if x3!=0.0
            h=cbrt(abs(6*quan/x3))
            q=x+h*h*h*x3/6
            q1=x1-x3*h*h/2   #*2
            q2=x2+h*x3
       else
            if x2!=0
                h=sqrt(abs(2*quan/x2))
                q=x-h*h*x2/2
                q1=x1+h*x2  
            else
                if x1!=0
                    h=abs(quan/x1)
                    q=x+h*x1
                else
                    q=x
                    
                end
                q1=x1
                
            end
            q2=x2
        end 
    end
    qv[i][0]=q
    qv[i][1]=q1 
    qv[i][2]=q2/2  
    nextStateTime[i]=simt+h
    return nothing
end


 



function Liqss_reComputeNextTime(::Val{3}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,a::Vector{Vector{Float64}} =#)
    q=qv[i][0];x=xv[i][0];q1=qv[i][1];x1=xv[i][1];x2=xv[i][2];q2=qv[i][2];x3=xv[i][3]
    coef=@SVector [q - x , q1-x1,q2-x2,-x3]# x and q might get away even further(change of sign) and we want that to be no more than another quan
    nextStateTime[i] = simt + minPosRoot(coef, Val(3))
    if abs(q-x) >= 2*quani # this happened when var i and j s turns are now...var i depends on j, j is asked here for next time...or if you want to increase quant*10 later it can be put back to normal and q & x are spread out by 10quan
        nextStateTime[i] = simt+1e-15
    else
        if q-x >0
            coef=setindex(coef, q-x-2*quantum[i],1)
            timetemp = simt + minPosRoot(coef, Val(3))
            if timetemp < nextStateTime[i] 
                nextStateTime[i]=timetemp
            end
        elseif  q-x <0
            coef=setindex(coef, q-x+2*quantum[i],1)
            timetemp = simt + minPosRoot(coef, Val(3))
            if timetemp < nextStateTime[i] 
                nextStateTime[i]=timetemp
            end
        else
            #nextStateTime[i]=simt+Inf#1e-19
        end
        if nextStateTime[i]<simt # this is coming from the fact that a variable can reach 2quan distance when it is not its turn, then computation above gives next=simt+(p-p)/dx...p-p should be zero but it can be very small negative
            nextStateTime[i]=simt+1e-15
        # @show simt,nextStateTime[i],i,x,q,quantum[i],xv[i][1]
        end
    end
    
end

