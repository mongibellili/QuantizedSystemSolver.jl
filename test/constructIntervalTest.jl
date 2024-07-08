acceptedi=Vector{Vector{Float64}}(undef,3)
for i =1:3
   acceptedi[i]=[0.0,0.0]#zeros(2)
end
constructIntrval(acceptedi,1.0,2.0,3.0,4.0)
@test acceptedi[1]  ==[0.0, 1.0]
@test acceptedi[2]==[2.0, 3.0]
@test acceptedi[3] ==[4.0, Inf]

constructIntrval(acceptedi,-1.0,2.0,3.0,4.0)
constructIntrval(acceptedi,1.0,-2.0,3.0,4.0)
constructIntrval(acceptedi,1.0,2.0,-3.0,4.0)
constructIntrval(acceptedi,1.0,2.0,3.0,-4.0)

constructIntrval(acceptedi,-1.0,-2.0,3.0,4.0)
constructIntrval(acceptedi,-1.0,2.0,-3.0,4.0)
constructIntrval(acceptedi,-1.0,2.0,3.0,-4.0)
constructIntrval(acceptedi,1.0,-2.0,-3.0,4.0)
constructIntrval(acceptedi,1.0,2.0,-3.0,-4.0)

constructIntrval(acceptedi,-1.0,-2.0,-3.0,4.0)
constructIntrval(acceptedi,-1.0,-2.0,3.0,-4.0)
constructIntrval(acceptedi,-1.0,2.0,-3.0,-4.0)
constructIntrval(acceptedi,1.0,-2.0,-3.0,-4.0)