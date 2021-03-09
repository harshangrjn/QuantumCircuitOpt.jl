function create_QCO_mip_formulation(m::QCOoptimizer)
    m.Qmodel = Model(get_option(m, :mip_solver)) # Construct JuMP Model

    build_QCO_mip_variables(m)
    build_QCO_mip_constraints(m)
    build_QCO_mip_objective(m)

    return
end

function solve_QCO_mip_formulation(m::QCOoptimizer)
    create_QCO_mip_formulation(m)
    
    start_time = time()
    optimize!(m.Qmodel)
    status = termination_status(m.Qmodel)
    m.sol_time = time() - start_time
    m.objective_value = objective_value(m.Qmodel)
    m.best_bound = objective_bound(m.Qmodel)
    m.solution = JuMP.value.(x)
end
