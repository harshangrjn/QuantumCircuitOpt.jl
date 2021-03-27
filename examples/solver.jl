"""
    get_solver(params)
    This function returns the JuMP optimizer with its appropriate attributes, based on 
        user-defined inputs in `params` dictionary.  
"""
function get_solver(params::Dict{String,Any})
    optimizers_list = ["cplex", "cbc", "ipopt", "juniper", "alpine"]

    if !(params["optimizer"] in optimizers_list)
        Memento.error(_LOGGER, "Specified optimizer does not belong to the predefined list. Add your optimizer separately with it's attributes")
    end

    # Mixed-integer programming optimizers
    if params["optimizer"] == "cplex"    # commercial 
       return get_cplex(params) 
    elseif params["optimizer"] == "cbc"  # open-source
       return get_cbc(params)
    
    # Local mixed-integer nonlinear programming optimizers
    elseif params["optimizer"] == "ipopt"    # open-source 
       return get_ipopt(params)
    elseif params["optimizer"] == "juniper"  # open-source 
       return get_juniper(params)
    
    # Global NLP/MINLP optimizer
    elseif params["optimizer"] == "alpine"   # open-source 
         alpine = optimizer_with_attributes(Alpine.Optimizer, 
                                            "nlp_solver" => get_ipopt(params),
                                            "minlp_solver" => get_juniper(params),  
                                            "mip_solver" => get_cplex(params),
                                            "presolve_bt" => false,
                                            "presolve_max_iter" => 10,
                                            "presolve_bp" => false,
                                            "disc_ratio" => 10)
        return alpine
    end
end

function get_cplex(params::Dict{String,Any})
     cplex = optimizer_with_attributes(CPLEX.Optimizer, 
                                       MOI.Silent() => !params["optimizer_log"], 
                                       "CPX_PARAM_PREIND" => params["presolve"]) 
    return cplex 
end

function get_ipopt(params::Dict{String,Any})
     ipopt = optimizer_with_attributes(Ipopt.Optimizer, 
                                       MOI.Silent() => !params["optimizer_log"], 
                                       "sb" => "yes", 
                                       "max_iter" => 1E4)
    return ipopt 
end

function get_cbc(params::Dict{String,Any})
     cbc = optimizer_with_attributes(Cbc.Optimizer, 
                                     MOI.Silent() => !params["optimizer_log"])
    return cbc 
end

function get_juniper(params::Dict{String,Any})
     juniper = optimizer_with_attributes(Juniper.Optimizer, 
                                         MOI.Silent() => !params["optimizer_log"], 
                                         "mip_solver" => get_cplex(params), 
                                         "nl_solver" => get_ipopt(params))
    return juniper
end
