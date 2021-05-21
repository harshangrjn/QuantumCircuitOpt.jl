"""
    variable_domain(var::JuMP.VariableRef)

Computes the valid domain of a given JuMP variable taking into account bounds
and the varaible's implicit bounds (e.g. binary).
"""
function variable_domain(var::JuMP.VariableRef)

  lb = -Inf
  (JuMP.has_lower_bound(var)) && (lb = JuMP.lower_bound(var))
  (JuMP.is_binary(var)) && (lb = max(lb, 0.0))

  ub = Inf
  (JuMP.has_upper_bound(var)) && (ub = JuMP.upper_bound(var))
  (JuMP.is_binary(var)) && (ub = min(ub, 1.0))

  return (lower_bound=lb, upper_bound=ub)
  
end

"""
    relaxation_bilinear(m::JuMP.Model, xy::JuMP.VariableRef, x::JuMP.VariableRef, y::JuMP.VariableRef)

general relaxation of binlinear term (McCormick), which can be used to obtain specific variants in partiuclar cases of variables (like binary)
```
z >= JuMP.lower_bound(x)*y + JuMP.lower_bound(y)*x - JuMP.lower_bound(x)*JuMP.lower_bound(y)
z >= JuMP.upper_bound(x)*y + JuMP.upper_bound(y)*x - JuMP.upper_bound(x)*JuMP.upper_bound(y)
z <= JuMP.lower_bound(x)*y + JuMP.upper_bound(y)*x - JuMP.lower_bound(x)*JuMP.upper_bound(y)
z <= JuMP.upper_bound(x)*y + JuMP.lower_bound(y)*x - JuMP.upper_bound(x)*JuMP.lower_bound(y)
```
"""
function relaxation_bilinear(m::JuMP.Model, xy::JuMP.VariableRef, x::JuMP.VariableRef, y::JuMP.VariableRef)
    lb_x, ub_x = variable_domain(x)
    lb_y, ub_y = variable_domain(y)
    
    if (lb_x == 0) && (ub_x == 1) && ((lb_y != 0) || (ub_y != 1))
      
      JuMP.@constraint(m, xy >= lb_y*x)
      JuMP.@constraint(m, xy >= y + ub_y*x - ub_y)
      JuMP.@constraint(m, xy <= ub_y*x)
      JuMP.@constraint(m, xy <= y + lb_y*x - lb_y)

    elseif (lb_y == 0) && (ub_y == 1) && ((lb_x != 0) || (ub_x != 1))
      
      JuMP.@constraint(m, xy >= lb_x*y)
      JuMP.@constraint(m, xy >= ub_x*y + x - ub_x)
      JuMP.@constraint(m, xy <= lb_x*y + x - lb_x)
      JuMP.@constraint(m, xy <= ub_x*y)

    elseif (lb_x == 0) && (ub_x == 1) && (lb_y == 0) && (ub_y == 1)
      
      JuMP.@constraint(m, xy <= x)
      JuMP.@constraint(m, xy <= y)
      JuMP.@constraint(m, xy >= x + y - 1)

    else
      
      JuMP.@constraint(m, xy >= lb_x*y + lb_y*x - lb_x*lb_y)
      JuMP.@constraint(m, xy >= ub_x*y + ub_y*x - ub_x*ub_y)
      JuMP.@constraint(m, xy <= lb_x*y + ub_y*x - lb_x*ub_y)
      JuMP.@constraint(m, xy <= ub_x*y + lb_y*x - ub_x*lb_y)

    end

    return m
  end