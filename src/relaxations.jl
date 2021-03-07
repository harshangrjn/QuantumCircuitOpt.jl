function mcc(m,xy,x,y)
    lb = [lower_bound(x), lower_bound(y)]
    ub = [upper_bound(x), upper_bound(y)]
    @constraint(m, xy >= lb[1]*y + lb[2]*x - lb[1]*lb[2])
    @constraint(m, xy >= ub[1]*y + ub[2]*x - ub[1]*ub[2])
    @constraint(m, xy <= lb[1]*y + ub[2]*x - lb[1]*ub[2])
    @constraint(m, xy <= ub[1]*y + lb[2]*x - ub[1]*lb[2])
    return m
  end
  
  function mcc_cont_bin(m,xy,x,y)
    lb = [lower_bound(x), lower_bound(y)]
    ub = [upper_bound(x), upper_bound(y)]
    @assert lb[2] == 0
    @assert ub[2] == 1
    @constraint(m, xy >= lb[1]*y)
    @constraint(m, xy >= ub[1]*y + x - ub[1])
    @constraint(m, xy <= lb[1]*y + x - lb[1])
    @constraint(m, xy <= ub[1]*y)
    return m
  end