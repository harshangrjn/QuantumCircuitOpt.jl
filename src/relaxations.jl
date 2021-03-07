function variable_domain(var::JuMP.VariableRef)
  lb = -Inf
  (JuMP.has_lower_bound(var)) && (lb = JuMP.lower_bound(var))
  (JuMP.is_binary(var)) && (lb = max(lb, 0.0))

  ub = Inf
  (JuMP.has_upper_bound(var)) && (ub = JuMP.upper_bound(var))
  (JuMP.is_binary(var)) && (ub = min(ub, 1.0))

  return (lower_bound=lb, upper_bound=ub)
end

function get_mccormick_relaxation(m::JuMP.Model, xy::JuMP.VariableRef, x::JuMP.VariableRef, y::JuMP.VariableRef)
    lb_x, ub_x = variable_domain(x)
    lb_y, ub_y = variable_domain(y)
    
    if (lb_x == 0) && (ub_x == 1)
      @constraint(m, xy >= lb_y*x)
      @constraint(m, xy >= y + ub_y*x - ub_y)
      @constraint(m, xy <= ub_y*x)
      @constraint(m, xy <= y + lb_y*x - lb_y)
    elseif (lb_y == 0) && (ub_y == 1)
      @constraint(m, xy >= lb_x*y)
      @constraint(m, xy >= ub_x*y + x - ub_x)
      @constraint(m, xy <= lb_x*y + x - lb_x)
      @constraint(m, xy <= ub_x*y)
    elseif (lb_x == 0) && (ub_x == 1) && (lb_y == 0) && (ub_y == 1)
      @constraint(m, xy <= x)
      @constraint(m, xy <= y)
      @constraint(m, xy >= x + y - 1)
    else
      @constraint(m, xy >= lb_x*y + lb_y*x - lb_x*lb_y)
      @constraint(m, xy >= ub_x*y + ub_y*x - ub_x*ub_y)
      @constraint(m, xy <= lb_x*y + ub_y*x - lb_x*ub_y)
      @constraint(m, xy <= ub_x*y + lb_y*x - ub_x*lb_y)
    end

    return m
  end