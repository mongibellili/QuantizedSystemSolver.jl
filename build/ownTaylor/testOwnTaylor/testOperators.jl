
#using qss

function testmul(t::Taylor0{Float64})
    t1=Taylor0([1.0,2.0,0.0],2)
    res=t*t1
end

function testmul()
    t1=Taylor0([1.0,2.0,0.0],2)
    res=t1*t1
end