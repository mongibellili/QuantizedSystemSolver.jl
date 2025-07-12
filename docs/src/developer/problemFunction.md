# Problem Function

*Constructing the Function Code:* After processing all ODEs and events,
the 'odeProblemFunc' function dynamically generates a Julia function (diffEqfunctionF) code needed to store
the system of ODEs and events. This code is built into a function that
handles different cases (i.e., which equation to evaluate based on an
index of a state change or an event).

## scoping and world age issues
There are several ways to construct this function:

 **case_p1:**
diffEqfunctionF=@RuntimeGeneratedFunction(diffEqfunctionExpression) # during the package's compilation time
 **case_p2:**
diffEqfunctionF=RuntimeGeneratedFunction(module,module,diffEqfunctionExpression) #generated at runtime in another module
 **case_p3:**
f1= Base.eval(module, diffEqfunctionExpression)
diffEqfunctionF= (args...) -> Base.invokelatest(f1, args...)
 **case_p4:**
diffEqfunctionF= Base.eval(module, diffEqfunctionExpression)
 **case_p5:**
diffEqfunctionF=mk_function(diffEqfunctionExpression)


*inside user space*

 **Case_u1:**  all on **top level**
```julia
function mainModel()
end
f=ODEProblem(mainModel)
solve(f)
```
**Case_u2:**  all on **local scope**
```julia
function test()
      function mainModel()
      end
      f=ODEProblem(mainModel)
      solve(f)
end
```
**Case_u3:**  mainModel function on **top level** and solve called from **local scope**  
```julia
function mainModel()
end
function test()
       f=ODEProblem(mainModel)
      solve(f)
end
```

### Discussion:
First of all, if there are no helperF, f can see anything it needs. For the structures of the 3 cases above, there is no issue of visibility. Below only world age issues are discussed:

 **Case_p1:** @RuntimeGeneratedFunction
✅ Case_u1, Case_u2 and Case_u3 are fine, because @RuntimeGeneratedFunction (macro version) **compiles early** at compile-time and thus not subject to world age.

 **Case_p2:** RuntimeGeneratedFunction
✅ Case_u1: works, because at top level, solve(f) is interpreted — not precompiled — so it sees everything defined up to the current world age.
❌ Case_u2 and Case_u3: an error is thrown (The applicable method may be too new: running in world age) because inside a function, solve(f) is compiled earlier — in the world it was defined — so it may not see newer methods added afterward.
```
t0: 📦 ODEProblem defined
t1: 🧑‍💻 function test() compiled (world age 100)
t2: 🧑‍💻 test() called
t3: 🧑‍💻 local mainModel defined
t4: ⚙️ f = ODEProblem(mainModel)
         └─ RuntimeGeneratedFunction(Main, Main, expr) adds method at world age 105
t5: ⚙️ solve(f)
         └─ fails: trying to call a method compiled at world age 105 from code compiled at 100
         └─ 🔥 WorldAgeError
```


**Case_p3:** Base.invokelatest
✅ Case_u1, Case_u2 and Case_u3 are fine, because it uses Julia's dynamic dispatch, bypassing the world age restrictions. However, it produces performance penalties.

 **case_p4:** Base.eval
This approach faces world age issues.

 **case_p5:** mk_function
Although mk_function elegantly avoids runtime eval, it comes with significant caveats:
  - Only supports unnamed function expressions — named or typed functions fail
  - Cannot include type annotations (e.g., u::Vector{Float64})
  - Does not support references to external functions
  - High latency & memory usage
  - Risk of segfaults when functions grow enormous

## Visibility of helper functions

### A helper function **OUTSIDE** the mainModel function:

**Case-by-case analysis:**

 
*inside user space*

 **Case_u1: **  Both functions live on **top level**
```julia
function helperF()
end
function mainModel()
    du=helperF()
end
f=ODEProblem(mainModel)
solve(f)
```
**Case_u2: **  Both function on **local scope**
```julia
function test()
      function helperF()
      end
      function mainModel()
           du=helperF()
      end
      f=ODEProblem(mainModel)
      solve(f)
end
```
**Case_u3: **  mainModel function on **local scope** and helper function on **top level**
```julia
function helperF()
end
function test()
      function mainModel()
           du=helperF()
      end
       f=ODEProblem(mainModel)
      solve(f)
end
```

### Discussion:

 **Case_p1:** @RuntimeGeneratedFunction
the f function is compiled  inside the package module, any helper function in the user scope is not seen during the evaluation of the package. 
Case_u1, Case_u2 and Case_u3 an error is thrown (helperF is undefined).

 **Case_p2:** RuntimeGeneratedFunction
The f function is compiled  inside the calling module (Main for example).
Case_u1: works.
Case_u2 and Case_u3: an error is thrown (The applicable method may be too new: running in world age)

 **Case_p3:** Base.invokelatest
Case_u1 and Case_u3: work.
Case_u2: works an error is thrown (UndefVarError: `helperF` not defined)

### A helper function **INSIDE** the mainModel function:

For example:
```julia
function foo(du,u,p,t)
  function bar(k, θ)
  end
  #differential equations
end
```
This can be accomplished in 2 ways:

1. function bar should be **placed** inside the main function as is, but we have the following issues:
  - @RuntimeGeneratedFunction does not allow closure: it gives :
   `$(Expr(:opaque_closure, :((k, θ)->begin`  ....instead of `function fname(k, θ)`
  - GeneralizedGenerated tricky to integrate
  - invokelatest is slow
  

2.  function bar should be **passed** to the main function as is.

### dependency not extracted from helper function

It is difficult to extract Jac and dependencies from the bar(k, θ) functions. For example:
```julia
function bar(k, θ)
         if 1 < k < 4
            return Y[k-1] * sin(θ[k] - θ[k-1])  # Contribution from (k, k-1)
         else
            return Y[k] * sin(θ[k] )  # Contribution from (k)
        end
       
    end
```
A solution that has not been implemented is : during parsing, everytime a function call g(i,j) is detected, @code_string g(i,j) must be used to get the body of the helper function.


## Conclusion


```julia
  if helper_functions_outside_model_definition
        warn_if_symbolic_jac
        if is_top_level                     
            RuntimeGeneratedFunction
        else
            warn_scope_performance
            invokelatest
        end
    else
      @RuntimeGeneratedFunction
    end
```