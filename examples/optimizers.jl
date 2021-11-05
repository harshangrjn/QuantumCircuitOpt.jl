#=====================================#
# MIP solvers (commercial, but fast)  #
#=====================================#

function get_gurobi()
    return JuMP.optimizer_with_attributes(Gurobi.Optimizer, 
                                          MOI.Silent() => false, 
                                          # "MIPFocus" => 3, # Focus on optimality over feasibility 
                                          "Presolve" => 1) 
end

function get_cplex()
     return JuMP.optimizer_with_attributes(CPLEX.Optimizer, 
                                           MOI.Silent() => false, 
                                           # "CPX_PARAM_EPGAP" => 1E-4,
                                           # "CPX_PARAM_MIPEMPHASIS" => 2 # Focus on optimality over feasibility 
                                           "CPX_PARAM_PREIND" => 1)
end

#======================================#
# MIP solvers (open-source, but slow)  #
#======================================#

function get_cbc()
    return JuMP.optimizer_with_attributes(Cbc.Optimizer, 
                                          MOI.Silent() => false) 
end

function get_glpk()
    return JuMP.optimizer_with_attributes(GLPK.Optimizer, 
                                          MOI.Silent() => false)
end

#========================================================#
# Continuous nonlinear programming solver (open-source)  #
#========================================================#

function get_ipopt()
     return JuMP.optimizer_with_attributes(Ipopt.Optimizer, 
                                       MOI.Silent() => true, 
                                       "sb" => "yes", 
                                       "max_iter" => Int(1E4))
end


#=================================================================#
# Local mixed-integer nonlinear programming solver (open-source)  #
#=================================================================#
function get_juniper()
     return JuMP.optimizer_with_attributes(Juniper.Optimizer, 
                                         MOI.Silent() => false, 
                                         "mip_solver" => get_gurobi(), 
                                         "nl_solver" => get_ipopt())
end

#=================================================================#
# Global mixed-integer nonlinear programming solver (open-source) #
#=================================================================#
function get_alpine()
        return JuMP.optimizer_with_attributes(Alpine.Optimizer, 
                                            "nlp_solver" => get_ipopt(),
                                            "minlp_solver" => get_juniper(),  
                                            "mip_solver" => get_gurobi(),
                                            "optimizer_presolve_bt" => false,
                                            "optimizer_presolve_max_iter" => 10,
                                            "optimizer_presolve_bp" => false,
                                            "disc_ratio" => 10)
end