"""
    get_solver(params)
    This function returns the JuMP optimizer with its appropriate attributes, based on 
        user-defined inputs in `params` dictionary.  
"""
function get_solver(params)
    optimizers_list = ["cplex", "cbc", "ipopt", "juniper", "alpine"]

    if !(params["optimizer"] in optimizers_list)
        Memento.error(_LOGGER, "Specified optimizer does not belong to the predefined list. Add your optimizer separately with it's attributes")
    end

    # Mixed-integer programming optimizers
    if params["optimizer"] == "cplex"    # commercial solver
       return get_cplex() 
    elseif params["optimizer"] == "cbc"  # open-source solver
       return get_cbc()
    
    # Local mixed-integer nonlinear programming optimizers
    elseif params["optimizer"] == "ipopt"    # open-source solver
       return get_ipopt()
    elseif params["optimizer"] == "juniper"  # open-source solver
       return get_juniper()
    
    # Global NLP/MINLP solver
    elseif params["optimizer"] == "alpine"   # open-source solver
        const alpine = optimizer_with_attributes(Alpine.Optimizer, 
                                                 "nlp_solver" => get_ipopt(),
                                                 "minlp_solver" => get_juniper(),  
                                                 "mip_solver" => get_cplex(),
                                                 "presolve_bt" => false,
                                                 "presolve_max_iter" => 10,
                                                 "presolve_bp" => false,
                                                 "disc_ratio" => 10)
        return alpine
    end
end

function get_cplex()
    const cplex = optimizer_with_attributes(CPLEX.Optimizer, 
                                            MOI.Silent() => !params["optimizer_log"], 
                                            "CPX_PARAM_PREIND" => params["presolve"]) 
    return cplex 
end

function get_ipopt()
    const ipopt = optimizer_with_attributes(Ipopt.Optimizer, 
                                            MOI.Silent() => !params["optimizer_log"], 
                                            "sb" => "yes", 
                                            "max_iter" => 1E4)
    return ipopt 
end

function get_cbc()
    const cbc = optimizer_with_attributes(Cbc.Optimizer, 
                                          MOI.Silent() => !params["optimizer_log"])
    return cbc 
end

function get_juniper()
    const juniper = optimizer_with_attributes(Juniper.Optimizer, 
                                              MOI.Silent() => !params["optimizer_log"], 
                                              "mip_solver" => get_cplex(), 
                                              "nl_solver" => get_ipopt())
    return juniper
end