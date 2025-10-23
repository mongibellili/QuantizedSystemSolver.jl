function test_ir_assignment(du,u,p,t)
      a,b=2.3,4.5
      @_inline x=foo(u[1])
      y=bar(2.0,a)
      @_inline z=[1.1,2.2,3.3]
      c=y[1]
      d,y[1]=1.0,2.0+b
      y[2]=4.0
      e,y[3]=p
      vec2=[1.0,2.0]
      g,h=[2.22,3.33],vec2    
      g[1]=5.0
      du[1]=a*u[2] + y[1]*u[1]+y[2]*h[1];
      du[2]= u[1] + x*u[2]+e+g[1]+z[2];
end  
tspan=(0.0,2.0)
u = [-0.5, 0.0,1.0] 
ir_group=IR(test_ir_assignment, u, tspan,inline_mode=AUTO) 
all_statements=ir_group.ir.statements

@test all_statements[1].lhs == :a
@test isapprox(all_statements[1].rhs, 2.3)
@test all_statements[1].keep_assignment== false

@test all_statements[2].lhs == :b
@test isapprox(all_statements[2].rhs, 4.5)
@test all_statements[2].keep_assignment== false

@test all_statements[3].lhs == :x
@test all_statements[3].rhs == :(foo((q[1])[0])) 
@test all_statements[3].keep_assignment == false

@test all_statements[4].lhs == :y
@test all_statements[4].rhs == :(bar(2.0, 2.3))
@test all_statements[4].keep_assignment == false

@test all_statements[5].lhs == :z
@test all_statements[5].rhs == :([1.1, 2.2, 3.3])
@test all_statements[5].keep_assignment == false

@test all_statements[6].lhs == :c
@test all_statements[6].rhs == :((bar(2.0, 2.3))[1])
@test all_statements[6].keep_assignment == true

@test all_statements[7].lhs == :d
@test all_statements[7].rhs == 1.0
@test all_statements[7].keep_assignment == false

@test all_statements[8].lhs == :(y[1])
@test all_statements[8].rhs == :(2.0 + 4.5)
@test all_statements[8].keep_assignment == true

@test all_statements[9].lhs == :(y[2])
@test all_statements[9].rhs == 4.0
@test all_statements[9].keep_assignment == true

@test all_statements[10].lhs == :((e, y[3])) 
@test all_statements[10].rhs == :p
@test all_statements[10].keep_assignment == true


@test all_statements[11].lhs == :vec2
@test all_statements[11].rhs == :([1.0, 2.0])
@test all_statements[11].keep_assignment == true

@test all_statements[12].lhs == :g
@test all_statements[12].rhs == :([2.22, 3.33])
@test all_statements[12].keep_assignment == true

@test all_statements[13].lhs == :h
@test all_statements[13].rhs == :vec2
@test all_statements[13].keep_assignment == false

@test all_statements[14].lhs == :(g[1])
@test all_statements[14].rhs == 5.0
@test all_statements[14].keep_assignment == true


#= function foo(i)
  i+1.0
end
function bar(a,b)
  [a+b, a-b, a*b]
end =#

#= :a = 2.3 (optimized out)
:b = 4.5 (optimized out)
:x = :(foo((q[1])[0])) (optimized out)
:y = :(bar(2.0, 2.3))
:z = :([1.1, 2.2, 3.3]) (optimized out)
:c = :(y[1]) (optimized out)
:d = 1.0 (optimized out)
:(y[1]) = :(2.0 + 4.5)
:(y[2]) = 4.0
:((e, y[3])) = :p
:vec2 = :([1.0, 2.0])
:g = :([2.22, 3.33])
:h = :vec2 (optimized out)
:(g[1]) = 5.0 =#
