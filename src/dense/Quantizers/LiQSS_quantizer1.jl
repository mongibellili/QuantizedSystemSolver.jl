 #iters
#= function updateQ(::Val{1},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,av::Vector{Vector{Float64}} =#,exactA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
   # a=av[i][i]
     cacheA[1]=0.0;exactA(qv,d,cacheA,i,i,simt)
    
     a=cacheA[1]
 

    #=    if i==1
        if abs(a+1.1)>1e-3
            @show i,a
        end
    else
        if abs(a+20.0)>1e-3
            @show i,a
        end
    end =#

     q=qv[i][0];x=xv[i][0];x1=xv[i][1];
     qaux[i][1]=q
    #=  olddx[i][1]=x1 =#
     u=x1-a*q
     #uv[i][i][1]=u
     dx=x1
     dxaux[i][1]=x1
     h=0.0
     Δ=quantum[i]
     debugH=4.3;debugL=14.1
    if a !=0.0
        if dx==0.0
            dx=u+(q)*a
            if dx==0.0
                dx=1e-26
            end
        end
    #for order1 finding h is easy but for higher orders iterations are cheaper than finding exact h using a quadratic,cubic...
    #exacte for order1: h=-2Δ/(u+xa-2aΔ) or h=2Δ/(u+xa+2aΔ)
        h = ft-simt
        q = (x + h * u) /(1 - h * a)
        if (abs(q - x) >  1*quantum[i]) # removing this did nothing...check @btime later
          h = (abs( quantum[i] / dx));
          q= (x + h * u) /(1 - h * a)
        end
        while (abs(q - x) >  1*quantum[i]) 
          h = h * 0.99*(quantum[i] / abs(q - x));
          q= (x + h * u) /(1 - h * a)
        end
       
       #=  α=-(a*x+u)/a
        h1denom=a*(x+Δ)+u;h2denom=a*(x-Δ)+u
        if h1denom==0.0 h1denom=1e-26 end
        if h2denom==0.0 h2denom=1e-26 end
        h1=Δ/h1denom
        h2=-Δ/h2denom
        if a<0
           if α>Δ
               h=h1;q=x+Δ
           elseif α<-Δ
               h=h2;q=x-Δ
           else
               h=max(ft-simt,-1/a)
               q=(x+h*u)/(1-h*a)#1-h*a non neg because h > -1/a > 1/a
           end
        else #a>0
            if α>Δ
                h=h2;q=x-Δ
            elseif α<-Δ
                h=h1;q=x+Δ
            else
                if a*x+u>0
                    h=max(ft-simt,-(a*x+u+a*Δ)/(a*(a*x+u-a*Δ)))# midpoint between asymptote and delta
                else
                    h=max(ft-simt,-(a*x+u-a*Δ)/(a*(a*x+u+a*Δ)))
                end
                q=(x+h*u)/(1-h*a)#1-h*a non neg because values of h from graph show that h >1/a
            end
        end =#
    else
        dx=u
        if dx>0.0
            q=x+quantum[i]# 
        else
            q=x-quantum[i]
        end
        if dx!=0
        h=(abs(quantum[i]/dx))
        else
            h=Inf
        end
    end
    qv[i][0]=q
   # println("inside single updateQ: q & qaux[$i][1]= ",q," ; ",qaux[i][1])
   nextStateTime[i]=simt+h
    return h
end   =# 
  
 
 #analytic favor q-x
function updateQ(::Val{1},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,av::Vector{Vector{Float64}} =#,exactA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
    cacheA[1]=0.0;exactA(qv,d,cacheA,i,i,simt);a=cacheA[1]
     q=qv[i][0];x=xv[i][0];x1=xv[i][1];
     qaux[i][1]=q
     u=x1-a*q
     dxaux[i][1]=x1
     h=0.0
     Δ=quantum[i]
     if a !=0.0
         α=-(a*x+u)/a
        #=  h1denom=a*(x+Δ)+u;h2denom=a*(x-Δ)+u
         if h1denom==0.0 h1denom=1e-26 end
         if h2denom==0.0 h2denom=1e-26 end
         h1=Δ/h1denom
         h2=-Δ/h2denom =#
         h1denom=a*(x+Δ)+u;h2denom=a*(x-Δ)+u
        #=  if h1denom==0.0 h1denom=1e-26 end
         if h2denom==0.0 h2denom=1e-26 end =#
         h1=Δ/h1denom # h1denom ==0 only when case α==Δ h1 and h2 not used h=Inf. ..
         h2=-Δ/h2denom #idem
         if a<0
            if α>Δ
                h=h1;q=x+Δ
            elseif α<-Δ
                h=h2;q=x-Δ
            else
                h=Inf
                q=-u/a
            end
        else #a>0
            if α>Δ
                h=h2;q=x-Δ
            elseif α<-Δ
                h=h1;q=x+Δ
            else
                if a*x+u>0
                    q=x-Δ   #
                    h=h2
                else
                    q=x+Δ 
                    h=h1
                end
            end
        end
     else #a=0
         if x1>0.0
             q=x+Δ # 
         else
             q=x-Δ 
         end
         if x1!=0
            h=(abs(Δ /x1))
         else
             h=Inf
         end
     end
    qv[i][0]=q
    nextStateTime[i]=simt+h
    return h
end     
  
 #analytic favor q-x but care about h large
#= function updateQ(::Val{1},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,av::Vector{Vector{Float64}} =#,exactA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
    
    # a=av[i][i]
     cacheA[1]=0.0;exactA(qv,cacheA,i,i)
    
     a=cacheA[1]
 

  #=    if i==1
        if abs(a+1.1)>1e-3
            @show i,a
        end
    else
        if abs(a+20.0)>1e-3
            @show i,a
        end
    end =#

     q=qv[i][0];x=xv[i][0];x1=xv[i][1];
     qaux[i][1]=q
    #=  olddx[i][1]=x1 =#
     u=x1-a*q
     #uv[i][i][1]=u
     dx=x1
     dxaux[i][1]=x1
     h=0.0
     Δ=quantum[i]
     if a !=0.0
         if dx==0.0
             dx=u+(q)*a
             if dx==0.0
                 dx=1e-26
             end
         end
         α=-(a*x+u)/a
         h1denom=a*(x+Δ)+u;h2denom=a*(x-Δ)+u
         if h1denom==0.0 h1denom=1e-26 end
         if h2denom==0.0 h2denom=1e-26 end
         h1=Δ/h1denom
         h2=-Δ/h2denom
         if a<0
            if α>Δ
                h=h1;q=x+Δ
            elseif α<-Δ
                h=h2;q=x-Δ
            else
               #=  h=Inf
                if a*x+u>0
                    q=x+α
                else
                    q=x+α
                end =#
                
                
                   #=  h=-1/a
                    q=(x+h*u)/(1-h*a) =#
                
                   
                h=max(ft-simt,-1/a)
                q=(x+h*u)/(1-h*a)#1-h*a non neg because a<0
                
            end
        else #a>0
            if α>Δ
                h=h2;q=x-Δ
            elseif α<-Δ
                h=h1;q=x+Δ
            else
               #=  h=Inf
                if a*x+u>0
                    q=x+α
                else
                    q=x+α
                end =#
                
                #= if a*x+u>0
                    h=h2;q=x-Δ
                else
                    h=h1;q=x+Δ
                end =#
               #=  h=1/a+ft-simt # i have if simt>=ft break in intgrator
                q=(x+h*u)/(1-h*a) =#
               #=  if abs(α)<Δ/2
                    h=2/a
                    q=x+2α
                else
                    h=3/a
                    q=(x+h*u)/(1-h*a)
                end =#
                if a*x+u>0
                    h=max(ft-simt,-(a*x+u+a*Δ)/(a*(a*x+u-a*Δ)))# midpoint between asymptote and delta
                  #=   q=x-Δ   #
                    h=h2 =#
                else
                    h=max(ft-simt,-(a*x+u-a*Δ)/(a*(a*x+u+a*Δ)))
                  #=   q=x+Δ 
                    h=h1 =#
                end
                q=(x+h*u)/(1-h*a)#1-h*a non neg because values of h from graph show that h >1/a

            end
        end
           
           
      
        
    

     else
         dx=u
         if dx>0.0
             q=x+quantum[i]# 
         else
             q=x-quantum[i]
         end
         if dx!=0
         h=(abs(quantum[i]/dx))
         else
             h=Inf
         end
     end
     qv[i][0]=q
    #=  if simt>=0.22816661756287676
     dxithr=a*q+u
     @show dxithr
     end =#

    # println("inside single updateQ: q & qaux[$i][1]= ",q," ; ",qaux[i][1])
    nextStateTime[i]=simt+h
     return h
end   =#

#= 
function updateQ(::Val{1},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,av::Vector{Vector{Float64}} =#,exactA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
    
   # a=av[i][i]
    cacheA[1]=0.0;exactA(qv,cacheA,i,i)
   
    a=cacheA[1]
   quan=quantum[i]
    q=qv[i][0];x=xv[i][0];x1=xv[i][1];
    qaux[i][1]=q
   # olddx[i][1]=x1
    u=x1-a*q
    #uv[i][i][1]=u
    dx=x1
    dxaux[i][1]=x1
    h=0.0
    if a !=0.0
        if dx==0.0
            dx=u+(q)*a
            if dx==0.0
                dx=1e-26
            end
        end
    #for order1 finding h is easy but for higher orders iterations are cheaper than finding exact h using a quadratic,cubic...
    #exacte for order1: h=-2Δ/(u+xa-2aΔ) or h=2Δ/(u+xa+2aΔ)
        h = ft-simt
        q = (x + h * u) /(1 - h * a)
        if (abs(q - x) >  1*quan) # removing this did nothing...check @btime later
          h = (abs( quan / dx));
          q= (x + h * u) /(1 - h * a)
        end
        while (abs(q - x) >  1*quan) 
          h = h * 0.99*(quan / abs(q - x));
          q= (x + h * u) /(1 - h * a)
        end
 
    else
        dx=u
        if dx>0.0
            q=x+quan# 
        else
            q=x-quan
        end
        if dx!=0
        h=(abs(quan/dx))
        else
            h=Inf
        end
    end
    qv[i][0]=q
   # println("inside single updateQ: q & qaux[$i][1]= ",q," ; ",qaux[i][1])
   nextStateTime[i]=simt+h
    return h
end   =#






function Liqss_reComputeNextTime(::Val{1}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,a::Vector{Vector{Float64}} =#)
    dt=0.0; q=qv[i][0];x=xv[i][0];x1=xv[i][1]
    if abs(q-x) >= 2*quantum[i] # this happened when var i and j s turns are now...var i depends on j, j is asked here for next time...or if you want to increase quant*10 later it can be put back to normal and q & x are spread out by 10quan
        nextStateTime[i] = simt+1e-15
    else
                if x1 !=0.0 #&& abs(q-x)>quantum[i]/10
                    dt=(q-x)/x1
                    if dt>0.0
                        nextStateTime[i]=simt+dt# later guard against very small dt
                    elseif dt<0.0
                        if x1>0.0  
                            nextStateTime[i]=simt+(q-x+2*quantum[i])/x1
                        else
                            nextStateTime[i]=simt+(q-x-2*quantum[i])/x1
                        end
                    end
                else
                    nextStateTime[i]=Inf
                end
    end
    if nextStateTime[i]<=simt 
        nextStateTime[i]=simt+Inf#1e-12
    end
  # this is coming from the fact that a variable can reach 2quan distance when it is not its turn, then computation above gives next=simt+(p-p)/dx...p-p should be zero but it can be very small negative
end





  #= function updateQ(::Val{1},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,av::Vector{Vector{Float64}} =#,exactA::Function,cacheA::MVector{1,Float64},dxaux::Vector{MVector{1,Float64}},qaux::Vector{MVector{1,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
    
    # a=av[i][i]
     cacheA[1]=0.0;exactA(qv,cacheA,i,i)
    
     a=cacheA[1]
 

  #=    if i==1
        if abs(a+1.1)>1e-3
            @show i,a
        end
    else
        if abs(a+20.0)>1e-3
            @show i,a
        end
    end =#

     q=qv[i][0];x=xv[i][0];x1=xv[i][1];
     qaux[i][1]=q
    #=  olddx[i][1]=x1 =#
     u=x1-a*q
     #uv[i][i][1]=u
     dx=x1
     dxaux[i][1]=x1
     h=0.0
     if a !=0.0
         if dx==0.0
             dx=u+(q)*a
             if dx==0.0
                 dx=1e-26
             end
         end
     #for order1 finding h is easy but for higher orders iterations are cheaper than finding exact h using a quadratic,cubic...
     #exacte for order1: h=-2Δ/(u+xa-2aΔ) or h=2Δ/(u+xa+2aΔ)
        #=  h = (ft+100.0-simt)
         q = (x + h * u) /(1 - h * a)
         if (abs(q - x) >  quantum[i]) =#
            h1=-1.0
            h0 = (ft-simt)
            q0 = (x + h0 * u) /(1 - h0 * a)
            if (abs(q - x) <  quantum[i])
                h1=h0
            end
        # end
        #=  if (abs(q - x) >  2*quantum[i]) # removing this did nothing...check @btime later
           h = (abs( quantum[i] / dx));
           q= (x + h * u) /(1 - h * a)
         end =#
        #=  while (abs(q - x) >  2*quantum[i]) 
           h = h * 1.98*(quantum[i] / abs(q - x));
           q= (x + h * u) /(1 - h * a)
         end =#
         coefQ=1
         #if (abs(q - x) >  coefQ*quantum[i])
            h2=coefQ*quantum[i]/(a*(x+coefQ*quantum[i])+u)
            
            if h2<0
                h2=-coefQ*quantum[i]/(a*(x-coefQ*quantum[i])+u)
                q=x-coefQ*quantum[i]
            else
                q=x+coefQ*quantum[i]
            end
           
         #end
        if h2<0
            h=Inf
        end
        h=max(h1,h2)
       #=   if h==ft+100.0-simt
            @show i,simt,x,x1,q,a*q+u
            
         end =#

     else
         dx=u
         if dx>0.0
             q=x+quantum[i]# 
         else
             q=x-quantum[i]
         end
         if dx!=0
         h=(abs(quantum[i]/dx))
         else
             h=Inf
         end
     end
     qv[i][0]=q
    # println("inside single updateQ: q & qaux[$i][1]= ",q," ; ",qaux[i][1])
    nextStateTime[i]=simt+h
     return nothing
end  
=#
 
    