using Symbolics, LinearAlgebra

@show 2
@variables x,il1,il2,rd1,rd2,rr,rpr,v,is1,is2,I
 A=[5.1e-6+4.2e-9+0.453e-6*x 4.2e-9+0.453e-6*x    ;4.2e-9+0.453e-6*x  5.1e-6+4.2e-9+0.453e-6*x  ]
B=(inv(A))
display(B)

C=[-(rd1+3.88e-3)*il1-(rr+rpr+0.453e-6*v)*I+rd1*is1;-(rd2+3.88e-3)*il2-(rr+rpr+0.453e-6*v)*I+rd2*is2;]

#display(C)
D=B*C
display(D)
#= D1=simplify(D[1], expand=true)
display(D1);println()
D2=simplify(D[2], expand=true)
display(D2);println()
D3=simplify(D[3], expand=true)
display(D3);println() =#

#= D=[5.1e-6 0 0;0 5.1e-6 0;0 0 5.1e-6]
E=inv(D)
display(E) =#
#= C=substitute.(B, (Dict(x=> c),))
display(C) =#


