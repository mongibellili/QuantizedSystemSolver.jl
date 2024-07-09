acceptedi=Vector{Vector{Float64}}(undef,3)
for i =1:3
   acceptedi[i]=[0.0,0.0]#zeros(2)
end
constructIntrval(acceptedi,1.0,2.0,3.0,4.0)
@test acceptedi[1]  ==[0.0, 1.0]
@test acceptedi[2]==[2.0, 3.0]
@test acceptedi[3] ==[4.0, Inf] 

constructIntrval(acceptedi,1.0,2.0,4.0,3.0)
constructIntrval(acceptedi,1.0,4.0,3.0,2.0)

constructIntrval(acceptedi,2.0,1.0,3.0,4.0)
constructIntrval(acceptedi,2.0,1.0,4.0,3.0)
constructIntrval(acceptedi,4.0,1.0,3.0,1.0)

constructIntrval(acceptedi,-1.0,2.0,3.0,4.0)
constructIntrval(acceptedi,-1.0,4.0,3.0,2.0)
constructIntrval(acceptedi,-1.0,2.0,4.0,3.0)
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

simt, x, β, c, α, b=0.05,1.0,3.5,-0.5,1.0,-6.0
getQfromAsymptote(simt, x, β, c, α, b) 
simt, x, β, c, α, b=0.05,1.0,0.0,-0.5,1.0,-6.0
getQfromAsymptote(simt, x, β, c, α, b) 
simt, x, β, c, α, b=0.05,1.0,0.0,0.0,1.0,-6.0
getQfromAsymptote(simt, x, β, c, α, b) 
simt, x, β, c, α, b=0.05,1.0,0.0,0.0,0.0,-6.0
getQfromAsymptote(simt, x, β, c, α, b) 
h, xi, quani, xj, quanj, aii, ajj, aij, aji, uij, uji=1.0,0.1,0.0001,-3.0,0.003,2.0,0.5,10.0,1.5,1.0,1.0
iterationH(h, xi, quani, xj, quanj, aii, ajj, aij, aji, uij, uji)