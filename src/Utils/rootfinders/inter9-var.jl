#using BenchmarkTools

@inline function smallest(args::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds m = args[1]
	for i in 2:N
		@inbounds arg = args[i]
		m = arg < m ? arg : m
	end
	m
end

@inline function horner(x::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	for i in 2:N
		@inbounds coeff = coeffs[i]
		v = muladd(x, v, coeff)
	end
	v
end

@inline function horner2(x::F, y::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	w = v
	for i in 2:N
		@inbounds coeff = coeffs[i]
		v = muladd(x, v, coeff)
		w = muladd(y, w, coeff)
	end
	v, w
end

@inline function hornerd(x::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds v = coeffs[1]
	d = zero(F)
	for i in 2:N
		@inbounds coeff = coeffs[i]
		d = muladd(x, d, v)
		v = muladd(x, v, coeff)
	end
	v, d
end

@inline function intervalhorner(low::F, high::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds colow = coeffs[1]
	cohigh = colow
	if low < zero(F)
		for i in 2:N
			@inbounds coeff = coeffs[i]
			colow, cohigh = if colow > zero(F)
				muladd(low, cohigh, coeff), muladd(high, colow, coeff)
			elseif cohigh < zero(F)
				muladd(high, cohigh, coeff), muladd(low, colow, coeff)
			else
				muladd(low, cohigh, coeff), muladd(low, colow, coeff)
			end
		end
	else
		for i in 2:N
			@inbounds coeff = coeffs[i]
			colow, cohigh = if colow > zero(F)
				muladd(low, colow, coeff), muladd(high, cohigh, coeff)
			elseif cohigh < zero(F)
				muladd(high, colow, coeff), muladd(low, cohigh, coeff)
			else
				muladd(high, colow, coeff), muladd(high, cohigh, coeff)
			end
		end	
	end
    colow, cohigh
end

@inline function posintervalhorner(low::F, high::F, coeffs::NTuple{N,F}) where {N, F <: AbstractFloat}
	@inbounds colow = coeffs[1]
	cohigh = colow
	for i in 2:N
		@inbounds coeff = coeffs[i]
		colow, cohigh = if colow > zero(F)
			muladd(low, colow, coeff), muladd(high, cohigh, coeff)
		elseif cohigh < zero(F)
			muladd(high, colow, coeff), muladd(low, cohigh, coeff)
		else
			muladd(high, colow, coeff), muladd(high, cohigh, coeff)
		end
	end
	colow, cohigh
end

function smallestpositiverootintervalnewtonregulafalsi(coeffs::NTuple{N,F}, doms::Ptr{NTuple{2,F}}) where {N, F <: AbstractFloat}
	if N == 1
		return typemax(F)
	elseif N == 2
		@inbounds ret = -coeffs[2] / coeffs[1]
		if ret < zero(F)
			return typemax(F)
		else
			return ret
		end
	end
    @inbounds _coeff1 = inv(coeffs[1])
    poly = ntuple(N) do i
        @inbounds _coeff1 * coeffs[i]
    end
    poly′ = ntuple(N - 1) do i
        @inbounds (N-i) * poly[i]
    end
	MM = smallest(poly)
	if MM > zero(F) return typemax(F) end
	domlow, domhigh = zero(F), one(F) - MM
	index = 0
    while true
		#@show domlow domhigh
		mid = 0.5(domlow + domhigh)
		comid = horner(mid, poly)
		codom′low, codom′high = posintervalhorner(domlow, domhigh, poly′)
		#@show comid codom′low codom′high
		if codom′low < zero(F) < codom′high
			leftlow, lefthigh, rightlow, righthigh = if comid < zero(F)
				domlow, mid - comid / codom′low, mid - comid / codom′high, domhigh
			else
				domlow, mid - comid / codom′high, mid - comid / codom′low, domhigh
			end
			#@show leftlow lefthigh rightlow righthigh
			if leftlow < lefthigh
				if rightlow < righthigh
					index += 1
					unsafe_store!(doms, (rightlow, righthigh), index)
				end
				domlow, domhigh = leftlow, lefthigh
				continue
			elseif rightlow < righthigh
				domlow, domhigh = rightlow, righthigh
				continue
			end
		else
			codomlow, codomhigh = horner2(domlow, domhigh, poly)
			#@show domlow domhigh codomlow codomhigh
			if codomlow * codomhigh < zero(F)
				while true
					comid, comid′ = hornerd(mid, poly)
					delta = comid / comid′
					newmid = mid - delta
					if abs(delta) < 1.0e-8mid 
						return newmid
					elseif domlow < newmid < domhigh
						mid = newmid
					else
						if comid * codomlow < zero(F)
							domhigh, codomhigh = mid, comid
						else
							domlow, codomlow = mid, comid
						end
						mid = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regula falsi
					end
				end
			end
		end
		if index == 0 break end
		domlow, domhigh = unsafe_load(doms, index)
		index -= 1
	end
	return typemax(F)
end

function allrealrootintervalnewtonregulafalsi(coeffs::NTuple{N,F}, res::Ptr{F}, doms::Ptr{NTuple{2,F}}) where {N, F <: AbstractFloat}
	if N == 1
		unsafe_store!(res, typemax(F), 1)
		return 1#typemax(F)
	elseif N == 2
		if coeffs[1]!=zero(F)
		  @inbounds ret = -coeffs[2] / coeffs[1]
		  unsafe_store!(res, ret, 1)
		  return 1#ret
		else
			unsafe_store!(res, typemax(F), 1)
			return 1#typemax(F)
		end
	end
	if coeffs[1]==0.0 #abs(coeffs[1])<1e-25
		nCoef=coeffs[2:end]
		#Base.GC.enable(false)
		ret=allrealrootintervalnewtonregulafalsi(nCoef, res, doms)
		#Base.GC.enable(true)
		return ret
    end
    @inbounds _coeff1 = inv(coeffs[1])
    poly = ntuple(N) do i
        @inbounds _coeff1 * coeffs[i]
    end
    poly′ = ntuple(N - 1) do i
        @inbounds (N-i) * poly[i]
    end
	mm = zero(F)
	MM = zero(F)
	s = -one(F)
	for i in 2:N
		@inbounds coeff = poly[i]
		if coeff < MM; MM = coeff end
		_coeff = s * coeff
		if _coeff < mm; mm = _coeff end
		s = -s
	end
	if mm == zero(F) && MM == zero(F) return 0 end
	index = 0
	domlow, domhigh = if mm < zero(F)
		if MM < zero(F)
			index = 1
			unsafe_store!(doms, (zero(F), one(F) - MM), 1)
		end
		mm - one(F), zero(F)
	else
		zero(F), one(F) - MM
	end
    counter = 0
    while true
		#@show domlow domhigh
		mid = 0.5(domlow + domhigh)
		comid = horner(mid, poly)
		codom′low, codom′high = intervalhorner(domlow, domhigh, poly′)
		#@show mid comid codom′low codom′high
		if codom′low < zero(F) < codom′high
			leftlow, lefthigh, rightlow, righthigh = if comid < zero(F)
				domlow, mid - comid / codom′low, mid - comid / codom′high, domhigh
			else
				domlow, mid - comid / codom′high, mid - comid / codom′low, domhigh
			end
			#@show leftlow lefthigh rightlow righthigh
			if leftlow < lefthigh
				if rightlow < righthigh
					index += 1
					unsafe_store!(doms, (rightlow, righthigh), index)
				end
				domlow, domhigh = leftlow, lefthigh
				continue
			elseif rightlow < righthigh
				domlow, domhigh = rightlow, righthigh
				continue
			end
		else
			codomlow, codomhigh = horner2(domlow, domhigh, poly)
            #@show domlow domhigh codomlow codomhigh
			if codomlow * codomhigh < zero(F)
				while true
					comid, comid′ = hornerd(mid, poly)
					delta = comid / comid′
                    #@show mid comid comid′
					newmid = mid - delta
					if abs(delta) < 1.0e-14abs(mid) 
						counter += 1
                        unsafe_store!(res, newmid, counter)
                        break
					elseif domlow < newmid < domhigh
						mid = newmid
					else
						if comid * codomlow < zero(F)
							domhigh, codomhigh = mid, comid
						else
							domlow, codomlow = mid, comid
						end
						mid = (domlow*codomhigh - domhigh*codomlow) / (codomhigh - codomlow) # regula falsi
					end
				end
			end
		end
		if index == 0 || counter == N-1 break end
		domlow, domhigh = unsafe_load(doms, index)
		index -= 1
	end
	return counter
end








function classicRoot(coeff::NTuple{3,Float64}) # credit goes to github.com/CIFASIS/qss-solver
    mpr=(-1.0,-1.0) #
	a=coeff[1];b=coeff[2];c=coeff[3]
    if a == 0.0 #|| (10000 * abs(a)) < abs(b)# coef3 is the coef of t^2
        if b != 0.0
			if -c / b>0.0
			 mpr = (-c / b,-1.0)
			end
		  end
       
    else 
       #double disc;
        disc = b * b - 4.0 * a * c#b^2-4ac
        if disc > 0.0 # no real roots
         
        
          #double sd, r1;
          sd = sqrt(disc);
		  r1 = (-b + sd) / (2.0 * a);
       
          r2 = (-b - sd) / (2.0 * a);

			mpr = (r1,r2)
		
		elseif disc == 0.0
			r1 = (-b ) / (2.0 * a);
			mpr = (r1,-1.0)
		end
        
    end
    return mpr
end


function quadRootv2(coeff::NTuple{3,Float64}) # 
	mpr=(-1.0,-1.0) #size 2 to use mpr[2] in quantizer
	a=coeff[1];b=coeff[2];c=coeff[3]
	if a == 0.0 || (10000 * abs(a)) < abs(b)# coef3 is the coef of t^2
		if b != 0.0
		  if -c / b>0.0
		   mpr = (-c / b,-1.0)
		  end
		end
	elseif b==0.0
		if -c/a>0
		mpr = (sqrt(-c / a),-1.0)
		end
	elseif c==0.0
		mpr=(-1.0,-b/a)
	else 
	   #double disc;
	   Δ = 1.0 - 4.0*c*a / (b*b)
		if Δ >0.0
			#= q = -0.5*(1.0+sign(b)*sqrt(Δ))*b
			r1 = q / a
		   
			r2=c / q =#
			sq=sqrt(Δ)
			r1=-0.5*(1.0+sq)*b/a
			r2=-0.5*(1.0-sq)*b/a
		 
			mpr = (r1,r2)
		elseif Δ ==0.0
			r1=-0.5*b/a
			mpr = (r1,r1-1e-12)
		end
	end
	return mpr
  end



  function mprv2(coeff::NTuple{3,Float64}) # 
	mpr=Inf
	a=coeff[1];b=coeff[2];c=coeff[3]
	if a == 0.0 #|| (10000 * abs(a)) < abs(b)# coef3 is the coef of t^2
		if b != 0.0
		  if -c / b>0.0
		   mpr = -c / b
		  end
		end
	elseif b==0.0
		if -c/a>0
		mpr = sqrt(-c / a)
		end
	elseif c==0.0
		if -b/a>0.0
			mpr = -b/a
		end
		
	else 
	   #double disc;
	   Δ = 1.0 - 4.0*c*a / (b*b)
		if Δ >0.0
			#= q = -0.5*(1.0+sign(b)*sqrt(Δ))*b
			r1 = q / a
		   
			r2=c / q =#
			sq=sqrt(Δ)
			r1=-0.5*(1.0+sq)*b/a
			if r1 > 0 
				mpr = r1;
			end
			r1=-0.5*(1.0-sq)*b/a

			if ((r1 > 0) && (r1 < mpr)) 
				mpr = r1;
			  end
		elseif Δ ==0.0
			r1=-0.5*b/a
			if r1 > 0 
				mpr = r1;
			end
		end
	end
	

	return mpr
  end


 #=    coef=NTuple{3,Float64}((-2.5768276401549883e-12, -239999.2196676244,1735305783508325))
  x=mprv2(coef)
  @show x
  sl=coef[3]+x*coef[2]+coef[1]*x*x
@show sl =#
 # 1735305783508325, -239999.2196676244, -2.5768276401549883
             #=    coeffs2=NTuple{3,Float64}((0.7,23999.2196676244,-2.5768276401549883e-12))
				#coeffs2=NTuple{3,Float64}((-2.5768276401549883e-12,-239999.2196676244,17))
				
				#coeffs2=NTuple{3,Float64}((14.691504647354595,-747452.6968034876,1.0e-6))
						function iter(res1::Ptr{Float64}, pp::Ptr{NTuple{2,Float64}},coeffs2::NTuple{3,Float64})
											#coeffs2=NTuple{3,Float64}((1.0, -2.0, -1.06))
										
										pp=pointer(Vector{NTuple{2,Float64}}(undef, 7))
										res1 = pointer(Vector{Float64}(undef, 2))
										unsafe_store!(res1, -1.0, 1);unsafe_store!(res1, -1.0, 2)
										allrealrootintervalnewtonregulafalsi(coeffs2,res1,pp)
										resTup=(unsafe_load(res1,1),unsafe_load(res1,2))
										#@show resTup
										resfilterd=filter((x) -> x >0.0 , resTup)
										
									#	display(resfilterd)
						end
						function  anal1(coeffs2::NTuple{3,Float64})
							#coeffs2=NTuple{3,Float64}((1.0, -2.0, -1.06))
						    classicRoot(coeffs2)
						end
						function  anal2(coeffs2::NTuple{3,Float64})
							#coeffs2=NTuple{3,Float64}((1.0, -2.0, -1.06))
						    quadRootv2(coeffs2)
						end

						pp=pointer(Vector{NTuple{2,Float64}}(undef, 7))
						res1 = pointer(Vector{Float64}(undef, 2))
						@show iter(res1,pp,coeffs2)
						@show anal1(coeffs2)
						@show anal2(coeffs2) 
						x= smallestpositiverootintervalnewtonregulafalsi(coeffs2,pp)
						@show x
						sl=coeffs2[3]+x*coeffs2[2]+coeffs2[1]*x*x
					  @show sl  =#
#= @btime iter(res1,pp)
@btime anal1()
@btime anal2() =#


#= coeffs=NTuple{3,Float64}((29160.956861496, 67.56376290717117, 0.014452560693316113))
#coeffs2=NTuple{3,Float64}((201011.0777843986, 106.4755863000737, 0.014452560693316113))
coeffs2=NTuple{4,Float64}((1.021243315770821e8, -3914.116824448214, -0.040323840000000166,1e-6))
pp=pointer(Vector{NTuple{2,Float64}}(undef, 11))
res1 = pointer(Vector{Float64}(undef, 3))
res2 = pointer(Vector{Float64}(undef, 3)) =#
#= count=allrealrootintervalnewtonregulafalsi(coeffs2,res2,pp)
@show count
resTup2=(unsafe_load(res2,1),unsafe_load(res2,2),unsafe_load(res2,3))
@show resTup2
resfilterd2=filter((x) -> x >0.0 , resTup2)
display(resfilterd2) =#


#coefi=NTuple{7,Float64}((2.4855534831418154e12, 1.2375936757103965e15, -1.1620693333574219e9, -1.21270658290325e12, -4.129153524641228e7, -80.08000016083042, -4.0000000000000105e-5))

#coefi=NTuple{7,Float64}((2.4794107580075693e12, 1.228418015847665e15, -6.091596024738817e15, -1.5165486073419932e18, -4.545138347948389e12, -6.054155392961131e6, -3.024053636764291))
#= coefi=NTuple{7,Float64}((-2.5668078295588949092296481e+16, 1.0474193043858423245087585e+18,1.7114669426922051413515583e+13, -2.8508122708869649897226691e+14, -8.5524376124873223589469262e+08, -1140.3250641304001419575054, -0.00057016254603761695740615778))

pp=pointer(Vector{NTuple{2,Float64}}(undef, 7))
respp = pointer(Vector{Float64}(undef, 6))

unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2);unsafe_store!(respp, -1.0, 3);unsafe_store!(respp, -1.0, 4);unsafe_store!(respp, -1.0, 5);unsafe_store!(respp, -1.0, 6)
allrealrootintervalnewtonregulafalsi(coefi,respp,pp)
resi1,resi2,resi3,resi4,resi5,resi6=unsafe_load(respp,1),unsafe_load(respp,2),unsafe_load(respp,3),unsafe_load(respp,4),unsafe_load(respp,5),unsafe_load(respp,6) 

for h in(resi1,resi2,resi3,resi4,resi5,resi6)
	#@show h
	if h>0
			f(x)=coefi[1]*x^6+coefi[2]*x^5+coefi[3]*x^4+coefi[4]*x^3+coefi[5]*x^2+coefi[6]*x+coefi[7]
			g(x)=coefi[7]+coefi[6]*x+coefi[5]*x^2+coefi[4]*x^3+coefi[3]*x^4+coefi[2]*x^5+coefi[1]*x^6
			@show h
			@show f(h),g(h)
	end
end =#

#= 
coefi=NTuple{7,BigFloat}((-2.5668078295588949092296481e+16, 1.0474193043858423245087585e+18,1.7114669426922051413515583e+13, -2.8508122708869649897226691e+14, -8.5524376124873223589469262e+08, -1140.3250641304001419575054, -0.00057016254603761695740615778))

pp=pointer(Vector{NTuple{2,BigFloat}}(undef, 7))
respp = pointer(Vector{BigFloat}(undef, 6))
Base.GC.enable(false)
unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2);unsafe_store!(respp, -1.0, 3);unsafe_store!(respp, -1.0, 4);unsafe_store!(respp, -1.0, 5);unsafe_store!(respp, -1.0, 6)
allrealrootintervalnewtonregulafalsi(coefi,respp,pp)
resi1,resi2,resi3,resi4,resi5,resi6=unsafe_load(respp,1),unsafe_load(respp,2),unsafe_load(respp,3),unsafe_load(respp,4),unsafe_load(respp,5),unsafe_load(respp,6) 
Base.GC.enable(true)
for h in(resi1,resi2,resi3,resi4,resi5,resi6)
	#@show h
	if h>0
			f(x)=coefi[1]*x^6+coefi[2]*x^5+coefi[3]*x^4+coefi[4]*x^3+coefi[5]*x^2+coefi[6]*x+coefi[7]
			g(x)=coefi[7]+coefi[6]*x+coefi[5]*x^2+coefi[4]*x^3+coefi[3]*x^4+coefi[2]*x^5+coefi[1]*x^6
			@show h
			@show f(h),g(h)
	end
end =#