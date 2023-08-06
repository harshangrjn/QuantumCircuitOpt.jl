#===================================#
# MIP solvers (commercial, faster)  #
#===================================#

# https://github.com/jump-dev/Gurobi.jl
function get_gurobi()
    return JuMP.optimizer_with_attributes(Gurobi.Optimizer, 
                                          MOI.Silent() => false, 
                                        #   "MIPFocus" => 3, # Focus on optimality over feasibility 
                                          "Presolve" => 1) 
end

# https://github.com/jump-dev/CPLEX.jl
function get_cplex()
     return JuMP.optimizer_with_attributes(CPLEX.Optimizer, 
                                           MOI.Silent() => false, 
                                           # "CPX_PARAM_EPGAP" => 1E-4,
                                           # "CPX_PARAM_MIPEMPHASIS" => 2 # Focus on optimality over feasibility 
                                           "CPX_PARAM_PREIND" => 1)
end

#====================================#
# MIP solvers (open-source, slower)  #
#====================================#
# https://github.com/jump-dev/HiGHS.jl
function get_highs()
    return JuMP.optimizer_with_attributes(
        HiGHS.Optimizer,
        "presolve" => "on",
        "log_to_console" => true,
    )
end

# https://github.com/jump-dev/Cbc.jl
function get_cbc()
    return JuMP.optimizer_with_attributes(Cbc.Optimizer, 
                                          MOI.Silent() => false) 
end

# https://github.com/jump-dev/Glpk.jl
function get_glpk()
    return JuMP.optimizer_with_attributes(GLPK.Optimizer, 
                                          MOI.Silent() => false)
end

#========================================================#
# Continuous nonlinear programming solver (open-source)  #
#========================================================#
# https://github.com/jump-dev/Ipopt.jl
function get_ipopt()
     return JuMP.optimizer_with_attributes(Ipopt.Optimizer, 
                                       MOI.Silent() => true, 
                                       "sb" => "yes", 
                                       "max_iter" => Int(1E4))
end

#=================================================================#
# Local mixed-integer nonlinear programming solver (open-source)  #
#=================================================================#
# https://github.com/lanl-ansi/Juniper.jl
function get_juniper()
     return JuMP.optimizer_with_attributes(Juniper.Optimizer, 
                                         MOI.Silent() => false, 
                                         "mip_solver" => get_gurobi(), 
                                         "nl_solver" => get_ipopt())
end

#=================================================================#
# Global mixed-integer nonlinear programming solver (open-source) #
#=================================================================#
# https://github.com/lanl-ansi/Alpine.jl
function get_alpine()
        return JuMP.optimizer_with_attributes(Alpine.Optimizer, 
                                            "nlp_solver" => get_ipopt(),
                                            "minlp_solver" => get_juniper(),  
                                            "mip_solver" => get_gurobi()
        )
end