# This file is part of the TaylorSeries.jl Julia package, MIT license
#
# Parameters for HomogeneousPolynomial and TaylorN

#= 
"""
    ParamsTaylor0

DataType holding the current variable name for `Taylor0`.

**Field:**

- `var_name   :: String`  Names of the variables

These parameters can be changed using [`set_Taylor0_varname`](@ref)
""" =#
mutable struct ParamsTaylor0
    var_name   :: String
end

const _params_Taylor0_ = ParamsTaylor0("t")

#= """
    set_Taylor0_varname(var::String)

Change the displayed variable for `Taylor0` objects.
""" =#
set_Taylor0_varname(var::String) = _params_Taylor0_.var_name = strip(var)




# Control the display of the big 𝒪 notation
const bigOnotation = Bool[true]
const _show_default = [false]

#= """
    displayBigO(d::Bool) --> nothing

Set/unset displaying of the big 𝒪 notation in  the output
of `Taylor0` and `TaylorN` polynomials. The initial value is
`true`.
""" =#
displayBigO(d::Bool) = (bigOnotation[end] = d; d)
#= 
"""
    use_Base_show(d::Bool) --> nothing

Use `Base.show_default` method (default `show` method
in Base), or a custom display. The initial value is
`false`, so customized display is used.
""" =#
use_show_default(d::Bool) = (_show_default[end] = d; d)
