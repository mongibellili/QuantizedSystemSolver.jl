
function test_ir_assignment(du,u,p,t)
    a=u[2] 
    M = [a-1.0, a+1.0] 
    for k in 1:2
      a=3.2
      du[k] = M[1]+a
    end
    du[3] = foo(a)
    if t>2.0
      u[2] = a+ M[2]
    end
end  
tspan=(0.0,2.0)
u = [-0.5, 0.0,1.0] 
ir_group=IR(test_ir_assignment, u, tspan,inline_mode=AUTO)
all_statements=ir_group.ir.statements

for_body = all_statements[3].body
@test for_body[1].lhs == :a
@test isapprox(for_body[1].rhs, 3.2)
@test for_body[1].keep_assignment== false
@test for_body[2].lhs == :(du[k])
@test for_body[2].rhs == :(M[1] + 3.2)
@test for_body[2].keep_assignment== true



if_condition=all_statements[5].condition
@test if_condition == :(t - 2.0)
if_expression=all_statements[5].body  # body contains the whole expr: condition + body
if_body = if_expression.args[2] #  just the body of the if
@test if_body.args[1] == :(q[2] = q[2] + M[2])
