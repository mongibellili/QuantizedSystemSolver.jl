

@inline function smallest(args::NTuple{N,F}) where {N, F <: Float64}
	@inbounds m = args[1]
	for i in 2:N
		@inbounds arg = args[i]
		m = arg < m ? arg : m
	end
	m
end

@inline function horner(x::F, coeffs::NTuple{N,F}) where {N, F <: Float64}
	@inbounds v = coeffs[1]
	for i in 2:N
		@inbounds coeff = coeffs[i]
		v = muladd(x, v, coeff)
	end
	v
end

@inline function horner2(x::F, y::F, coeffs::NTuple{N,F}) where {N, F <: Float64}
	@inbounds v = coeffs[1]
	w = v
	for i in 2:N
		@inbounds coeff = coeffs[i]
		v = muladd(x, v, coeff)
		w = muladd(y, w, coeff)
	end
	v, w
end

@inline function hornerd(x::F, coeffs::NTuple{N,F}) where {N, F <: Float64}
	@inbounds v = coeffs[1]
	d = zero(F)
	for i in 2:N
		@inbounds coeff = coeffs[i]
		d = muladd(x, d, v)
		v = muladd(x, v, coeff)
	end
	v, d
end

@inline function intervalhorner(low::F, high::F, coeffs::NTuple{N,F}) where {N, F <: Float64}
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

@inline function posintervalhorner(low::F, high::F, coeffs::NTuple{N,F}) where {N, F <: Float64}
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


#function allrealrootintervalnewtonregulafalsi(coeffs::NTuple{N,F}, res::Ptr{F}, doms::Ptr{NTuple{2,F}}) where {N, F <: Float64}
function allrealrootintervalnewtonregulafalsi(coeffs::NTuple{N,F}, res::Ptr{F}, doms::Ptr{NTuple{2,F}}) where {N, F <: Float64}
	if N == 1
		return 0
	elseif N == 2
		#println("fix setIndex Ptr error")
		@inbounds res[1] = -coeffs[2] / coeffs[1]
		return 1
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
					if abs(delta) < 1.0e-15abs(mid) 
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


a,b,c,d,e,f,g=2.4794107580075693e12, 1.228418015847665e15, -6.091596024738817e15, -1.5165486073419932e18, -4.545138347948389e12, -6.054155392961131e6, -3.024053636764291

coefi=NTuple{7,Float64}((2.4794107580075693e12, 1.228418015847665e15, -6.091596024738817e15, -1.5165486073419932e18, -4.545138347948389e12, -6.054155392961131e6, -3.024053636764291))

#coefi=NTuple{7,Float64}((Float64(a), Float64(b), Float64(c), Float64(d), Float64(e), Float64(f), Float64(g)))
#coefi=NTuple{7,Float64}((Float64(2.4855534831418154e12), Float64(1.2375936757103965e15), Float64(-1.1620693333574219e9), Float64(-1.21270658290325e12), Float64(-4.129153524641228e7), Float64(-80.08000016083042), Float64(-4.0000000000000105e-5)))
#pp=pointer(Vector{NTuple{2,Float64}}(undef, 7))
#respp = pointer(Vector{Float64}(undef, 6))

pp=pointer(Vector{NTuple{2,Float64}}(undef, 7))
respp = pointer(Vector{Float64}(undef, 6))

unsafe_store!(respp, -1.0, 1);unsafe_store!(respp, -1.0, 2);unsafe_store!(respp, -1.0, 3);unsafe_store!(respp, -1.0, 4);unsafe_store!(respp, -1.0, 5);unsafe_store!(respp, -1.0, 6)
#= allrealrootintervalnewtonregulafalsi(coefi,respp,pp)
resi1,resi2,resi3,resi4,resi5,resi6=unsafe_load(respp,1),unsafe_load(respp,2),unsafe_load(respp,3),unsafe_load(respp,4),unsafe_load(respp,5),unsafe_load(respp,6) 

for h in(resi1,resi2,resi3,resi4,resi5,resi6)
ff(x)=coefi[1]*x^6+coefi[2]*x^5+coefi[3]*x^4+coefi[4]*x^3+coefi[5]*x^2+coefi[6]*x+coefi[7]
@show h
@show ff(h)
end =#

h=smallestpositiverootintervalnewtonregulafalsi(coefi,pp)

ff(x)=coefi[1]*x^6+coefi[2]*x^5+coefi[3]*x^4+coefi[4]*x^3+coefi[5]*x^2+coefi[6]*x+coefi[7]
@show h
@show ff(h)
