

#iterations
#= function updateQ(::Val{2},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64},exactA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{2,Float64}},qaux::Vector{MVector{2,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
    cacheA[1]=0.0;exactA(qv,d,cacheA,i,i);a=cacheA[1]
   # exactA(xv,cacheA,i,i);a=cacheA[1]
    q=qv[i][0] ;q1=qv[i][1]; x=xv[i][0];  x1=xv[i][1]; x2=xv[i][2]*2; #u1=uv[i][i][1]; u2=uv[i][i][2]
    qaux[i][1]=q+(simt-tq[i])*q1#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1                     #appears only here...updated here and used in updateQevent

    u1=x1-a*qaux[i][1]
    u2=x2-a*q1
    dxaux[i][1]=x1
    dxaux[i][2]=x2
   
    ddx=x2
    quan=quantum[i]
    h=0.0
    if a!=0.0
        if ddx ==0.0
             ddx=a*a*q+a*u1 +u2
            if ddx==0.0 
                ddx=1e-40# changing -40 to -6 nothing changed
               # println("ddx=0")
            end
        end
        h = ft-simt; q=(x-h*a*x-h*h*(a*u1+u2)/2)/(1 - h * a + h * h * a * a / 2)
               
        
        if (abs(q - x) > 2.0* quan) # removing this did nothing...check @btime later
          h = sqrt(abs(2*quan / ddx)) # sqrt highly recommended...removing it leads to many sim steps..//2* is necessary in 2*quan when using ddx
          #= q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                   (1 - h * a + h * h * a * a / 2) =#
                #   qtemp=(x-h*a*x-h*h*(a*u1+u2)/2)/(1 - h * a + h * h * a * a / 2)
                   q=x+(-h*h*(a*a*x+a*u1+u2)/2)/((1 - h * a + h * h * a * a / 2))
                  #=  if abs(q-qtemp)>1e-16
                    @show qtemp,q
                   end =#
                 
        end
        maxIter=10000
       # tempH=h
        while (abs(q - x) >2.0*  quan) && (maxIter>0) && (h>0)
            
          h = h *sqrt(quan / abs(q - x))
        #  h = h *0.99*(1.8*quan / abs(q - x))
          #= q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                   (1 - h * a + h * h * a * a / 2) =#
                 #  qtemp=(x-h*a*x-h*h*(a*u1+u2)/2)/(1 - h * a + h * h * a * a / 2)
                   q=x+(-h*h*(a*a*x+a*u1+u2)/2)/((1 - h * a + h * h * a * a / 2))
                   #= if abs(q-qtemp)>1e-16
                    @show q,qtemp
                   end =#
                  
          maxIter-=1
        end
        if maxIter==0
            println("updateQ maxIterReached")
        end
      #=   if  (abs(q - x) > 2* quan)
            coef=@SVector [quan, -a*quan,(a*a*(x+quan)+a*u1+u2)/2]#
                h1= minPosRoot(coef, Val(2))
              coef=@SVector [-quan, a*quan,(a*a*(x-quan)+a*u1+u2)/2]#
                h2= minPosRoot(coef, Val(2))
  
              if h1<h2
                  h=h1;q=x+quan
              else
                  h=h2;q=x-quan
              end
              qtemp=(x-h*a*x-h*h*(a*u1+u2)/2)/(1 - h * a + h * h * a * a / 2)
              if abs(qtemp-q)>1e-12
                println("error quad vs qexpression")
              end
              if h1!=Inf && h2!=Inf
                println("quadratic eq double mpr")
              end
  
              if h1==Inf && h2==Inf
                println("quadratic eq NO mpr")
              end


          end  =#
         #= 
                                                                                        if  (abs(q - x) > 2* quan)

                                                                                        # if x2<0.0  #x2 might changed direction...I should use ddx but again ddx=aq+u and q in unknown the sign is unknown
                                                                                              #=   pertQuan=quan-1e-13  #
                                                                                                q=x+pertQuan      # guess q to be thrown up         
                                                                                                coef=@SVector [pertQuan, -a*pertQuan,(a*a*(x+pertQuan)+a*u1+u2)/2]#
                                                                                                h= minPosRoot(coef, Val(2))
                                                                                               pertQuan=quan+1e-13         
                                                                                                coef=@SVector [pertQuan, -a*pertQuan,(a*a*(x+pertQuan)+a*u1+u2)/2]#
                                                                                                h2= minPosRoot(coef, Val(2))
                                                                                                if h2<h
                                                                                                    q=x+pertQuan 
                                                                                                    h=h2
                                                                                                end
                                                                                                ########### q to be thrown down
                                                                                             pertQuan=quan-1e-13
                                                                                                coef=@SVector [-pertQuan, a*pertQuan,(a*a*(x-pertQuan)+a*u1+u2)/2]#
                                                                                                h2= minPosRoot(coef, Val(2))
                                                                                                if h2<h
                                                                                                    q=x-pertQuan 
                                                                                                    h=h2
                                                                                                end
                                                                                                pertQuan=quan+1e-13
                                                                                                coef=@SVector [-pertQuan, a*pertQuan,(a*a*(x-pertQuan)+a*u1+u2)/2]#
                                                                                                h2= minPosRoot(coef, Val(2))
                                                                                                if h2<h
                                                                                                    q=x-pertQuan 
                                                                                                    h=h2
                                                                                                end =#
                                                                                             pertQuan=quan
                                                                                                coef=@SVector [-pertQuan, a*pertQuan,(a*a*(x-pertQuan)+a*u1+u2)/2]#
                                                                                                h= minPosRoot(coef, Val(2))
                                                                                                q=x-pertQuan
                                                                                               #=  if h2<h
                                                                                                    q=x-pertQuan 
                                                                                                    h=h2
                                                                                                end =#
                                                                                                coef=@SVector [pertQuan, -a*pertQuan,(a*a*(x+pertQuan)+a*u1+u2)/2]#
                                                                                                h2= minPosRoot(coef, Val(2))
                                                                                                if h2<h
                                                                                                    q=x+pertQuan 
                                                                                                    h=h2
                                                                                                end
                                                                                        # end 
                                                                                        end =#
                                                                                
             

         

         α1=1-h*a
        if abs(α1)==0.0
            α1=1e-30*sign(α1)
        end
        q1=(a*q+u1+h*u2)/α1  #later investigate 1=h*a


    else
        #println("a==0")
        if x2!=0.0  
           h=sqrt(abs(2*quan/x2))   #sqrt necessary with u2
           q=x-h*h*x2/2
           q1=x1+h*x2
        else
           # println("x2==0")
            if x1!=0.0
                h=abs(quan/x1)
                q=x+h*x1
                q1=0#x#/2
            else
                h=Inf
                q=x
                q1=x1
            end
        end 

    end
 
    qv[i][0]=q
    qv[i][1]=q1  
   nextStateTime[i]=simt+h

    return h
   
end   =#


#analytic
function updateQ(::Val{2},i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64},exactA::Function,d::Vector{Float64},cacheA::MVector{1,Float64},dxaux::Vector{MVector{2,Float64}},qaux::Vector{MVector{2,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})
    cacheA[1]=0.0;exactA(qv,d,cacheA,i,i,simt);a=cacheA[1]
   # exactA(xv,cacheA,i,i);a=cacheA[1]
    q=qv[i][0] ;q1=qv[i][1]; x=xv[i][0];  x1=xv[i][1]; x2=xv[i][2]*2; #u1=uv[i][i][1]; u2=uv[i][i][2]
    qaux[i][1]=q+(simt-tq[i])*q1#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1                     #appears only here...updated here and used in updateQevent




    u1=x1-a*qaux[i][1]
    u2=x2-a*q1
    dxaux[i][1]=x1
    dxaux[i][2]=x2
   
    ddx=x2
    quan=quantum[i]
    h=0.0
    if a!=0.0                                             
            # quan=quan*1.5
        if a*a*x+a*u1+u2<=0.0
            if -(a*a*x+a*u1+u2)/(a*a)<quan # asymptote<delta...no sol ...no need to check
                h=Inf
                q=-(a*u1+u2)/(a*a)#q=x+asymp
            else
                q=x+quan
                coefi=NTuple{3,Float64}(((a*a*(x+quan)+a*u1+u2)/2,-a*quan,quan))
                h=minPosRootv1(coefi) #needed for dq
                #math was used to show h cannot be 1/a
            end
        elseif a*a*x+a*u1+u2>=0.0
           if -(a*a*x+a*u1+u2)/(a*a)>-quan # asymptote>-delta...no sol ...no need to check
                h=Inf
                q=-(a*u1+u2)/(a*a)
           else
                coefi=NTuple{3,Float64}(((a*a*(x-quan)+a*u1+u2)/2,a*quan,-quan))
                h=minPosRootv1(coefi)
              #  h=mprv2(coefi)
                
                q=x-quan
            end
        else#a*a*x+a*u1+u2==0 -->f=0....q=x+f(h)=x+h*g(h) -->g(h)==0...dx==0 -->aq+u==0
            q=-u1/a
            h=Inf
        end
         

        if h!=Inf
            if h!=1/a
                 q1=(a*q+u1+h*u2)/(1-h*a) 
            else#h=1/a
                error("report bug: updateQ: h cannot=1/a; h=$h , 1/a=$(1/h)")
            end
        else #h==inf make ddx==0 dq=-u2/a
            q1=-u2/a
        end


       #=  if h!=Inf
         α1=1-h*a
        if abs(α1)==0.0
            α1=1e-30*sign(α1)
        end
        q1=(a*q+u1+h*u2)/α1  
        else #h==inf make ddx==0 dq=-u2/a
            q1=-u2/a
        end =#


    else
        #println("a==0")
        if x2!=0.0  
           h=sqrt(abs(2*quan/x2))   #sqrt necessary with u2
           q=x-h*h*x2/2
           q1=x1+h*x2
        else
           # println("x2==0")
            if x1!=0.0
                #quantum[i]=1quan
                h=abs(1*quan/x1)   # *10 just to widen the step otherwise it would behave like 1st order
                q=x+h*x1
                q1=0.0
            else
                h=Inf
                q=x
                q1=0.0
            end
        end 

    end
  #=   if 2.500745083181994e-5<=simt<3.00745083181994e-5
    @show simt,i,q,x,h,a
    @show a*a*x+a*u1+u2
    end =#
    qv[i][0]=q
    qv[i][1]=q1  
   nextStateTime[i]=simt+h
 
    return h
   
end    
#= 
function nupdateQ(::Val{2}#= ,cacheA::MVector{1,Float64},map::Function =#,i::Int, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64},av::Vector{Vector{Float64}},uv::Vector{Vector{MVector{O,Float64}}},qaux::Vector{MVector{O,Float64}},olddx::Vector{MVector{O,Float64}},tx::Vector{Float64},tq::Vector{Float64},simt::Float64,ft::Float64, nextStateTime::Vector{Float64})where{O}
    #a=getA(Val(Sparsity),cacheA,av,i,i,map)
    a=av[i][i]
    q=qv[i][0] ;q1=qv[i][1]; x=xv[i][0];  x1=xv[i][1]; x2=xv[i][2]*2; u1=uv[i][i][1]; u2=uv[i][i][2]
    qaux[i][1]=q+(simt-tq[i])*q1#appears only here...updated here and used in updateApprox and in updateQevent later
    qaux[i][2]=q1                     #appears only here...updated here and used in updateQevent
    olddx[i][1]=x1#appears only here...updated here and used in updateApprox   
    #u1=u1+(simt-tu[i])*u2 # for order 2: u=u+tu*deru  this is necessary deleting causes scheduler error
    u1=x1-a*qaux[i][1]
    uv[i][i][1]=u1
    uv[i][i][2]=x2-a*q1
    u2=uv[i][i][2]
    #tu[i]=simt  
    # olddx[i][2]=2*x2# 
    ddx=x2  
    quan=quantum[i]
    h=0.0
   #=  if simt == 0.004395600232045285
        @show a
        @show x1
        @show u1
        @show u2

    end =#
    if a!=0.0
        if ddx ==0.0
             ddx=a*a*q+a*u1 +u2
            if ddx==0.0 
                ddx=1e-40# changing -40 to -6 nothing changed
                #println("ddx=0")
            end
        end
        h = ft-simt
        #tempH1=h
        #q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /(1 - h * a + h * h * a * a / 2)
                 q=(x-h*a*x-h*h*(a*u1+u2)/2)/(1 - h * a + h * h * a * a / 2)
        
        if (abs(q - x) > 2* quan) # removing this did nothing...check @btime later
          h = sqrt(abs(2*quan / ddx)) # sqrt highly recommended...removing it leads to many sim steps..//2* is necessary in 2*quan when using ddx
          q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                   (1 - h * a + h * h * a * a / 2)
          
        end
        maxIter=1000
        tempH=h
        while (abs(q - x) > 2* quan) && (maxIter>0) && (h>0)
            
          h = h *sqrt(quan / abs(q - x))
          q = ((x + h * u1 + h * h / 2 * u2) * (1 - h * a) + (h * h / 2 * a - h) * (u1 + h * u2)) /
                   (1 - h * a + h * h * a * a / 2)
          maxIter-=1
        end

        q1=(a*q+u1+h*u2)/(1-h*a)  #later investigate 1=h*a
    else
        if x2!=0.0  
           h=sqrt(abs(2*quan/x2))   #sqrt necessary with u2
           q=x-h*h*x2/2
           q1=x1+h*x2

        else
           # println("x2==0")
            if x1!=0.0
                h=abs(quan/x1)
                q=x+h*x1
                q1=x1
            else
                h=Inf
                q=x
                q1=x1
            end
        end 
    end
 
    qv[i][0]=q
    qv[i][1]=q1  
   nextStateTime[i]=simt+h

    return h
end =#




function Liqss_reComputeNextTime(::Val{2}, i::Int, simt::Float64, nextStateTime::Vector{Float64}, xv::Vector{Taylor0},qv::Vector{Taylor0}, quantum::Vector{Float64}#= ,a::Vector{Vector{Float64}} =#)
    q=qv[i][0];x=xv[i][0];q1=qv[i][1];x1=xv[i][1];x2=xv[i][2];quani=quantum[i]
    β=0
    if abs(q-x) >= 2*quani # this happened when var i and j s turns are now...var i depends on j, j is asked here for next time...or if you want to increase quant*10 later it can be put back to normal and q & x are spread out by 10quan
        nextStateTime[i] = simt+1e-12
       #=  if simt==2.500745083181994e-5
        println("quantizer-recomputeNext: abs(q-x) >= 2*quani")
        @show simt,i,x,q,abs(q-x),quani
        end  =#
       #=  if simt>3.486047550372409
            @show djsksfs
        end =#
    #elseif q!=x#abs(q-x)>quani/10000.0# #q=x ..let it be caught later by 2*Quantum
    else
                    coef=@SVector [q-x#= -1e-13 =#, q1-x1,-x2]#

                    
                        nextStateTime[i] = simt + minPosRoot(coef, Val(2))
                    
                    if q-x >0.0#1e-9
                        coef=setindex(coef, q-x-2*quantum[i],1)
                        timetemp = simt + minPosRoot(coef, Val(2))
                        if timetemp < nextStateTime[i] 
                            nextStateTime[i]=timetemp
                        end
                    #=  if 0.000691<simt<0.00071
                            @show qv,xv,nextStateTime
                        end =#
                    elseif  q-x <0.0#-1e-9
                        coef=setindex(coef, q-x+2*quantum[i],1)
                        timetemp = simt + minPosRoot(coef, Val(2))
                        if timetemp < nextStateTime[i] 
                            nextStateTime[i]=timetemp
                        end
                        #= if 0.000691<simt<0.00071
                            @show coef
                            @show nextStateTime
                            @show simt,i,quantum[i]
                        end =#
                    #else# q=x ..let it be caught later by 2*Quantum
                    #  nextStateTime[i]=simt+Inf#1e-19
                        #nextStateTime[i]=simt+1e-19
                    #=  if q-x==0.0
                        nextStateTime[i]=simt+Inf#1e-19 #
                    else
                        nextStateTime[i]=simt+1e-12#Inf#1e-19 #
                    end =#
                    end


                    if nextStateTime[i]<=simt # this is coming from the fact that a variable can reach 2quan distance when it is not its turn, then computation above gives next=simt+(p-p)/dx...p-p should be zero but it can be very small negative
                        nextStateTime[i]=simt+Inf#1e-14
                    # @show simt,nextStateTime[i],i,x,q,quantum[i],xv[i][1]
                    end
                #=  if 2.500745083181994e-5<=simt<3.00745083181994e-5
                        @show simt,i,nextStateTime,x,q
                    end  =#
    #else
      #=   adavnce=nextStateTime[i]-simt
        @show i,abs(q-x),quani/100.0,nextStateTime[i],simt,adavnce =#

    end
end
 

