""
function build_QCModel_result(qcm::QuantumCircuitModel, solve_time::Number)
    # try-catch is needed until solvers reliably support ResultCount()
    result_count = 1
    try
        result_count = JuMP.result_count(qcm.model)
    catch
        Memento.warn(_LOGGER, "the given optimizer does not provide the ResultCount() attribute, assuming the solver returned a solution which may be incorrect.");
    end

    solution = Dict{String,Any}()

    if result_count > 0
        solution = QCO.build_QCModel_solution(qcm)
    else
        Memento.warn(_LOGGER, "Quantum circuit model has no results, solution cannot be built")
    end

    result = Dict{String,Any}(
        "optimizer" => JuMP.solver_name(qcm.model),
        "termination_status" => JuMP.termination_status(qcm.model),
        "primal_status" => JuMP.primal_status(qcm.model),
        "objective" => QCO.get_objective_value(qcm.model),
        "objective_lb" => QCO.get_objective_bound(qcm.model),
        "solve_time" => solve_time,
        "solution" => solution,
    )

    # Objective slack Penalty
    if qcm.data["decomposition_type"] == "approximate"
        result["objective_slack_penalty"] = qcm.options.objective_slack_penalty
    end

    return result
end

""
function get_objective_value(model::JuMP.Model)
    obj_val = NaN

    try
        obj_val = JuMP.objective_value(model)
    catch
        Memento.warn(_LOGGER, "Objective value is unbounded. Problem may be infeasible or not constrained properly");
    end

    return obj_val
end


""
function get_objective_bound(model::JuMP.Model)
    obj_lb = -Inf

    try
        obj_lb = JuMP.objective_bound(model)
    catch
    end

    return obj_lb
end

function build_QCModel_solution(qcm::QuantumCircuitModel)
    solution = Dict{String,Any}()
    for i in keys(qcm.variables)
        solution[String(i)] = JuMP.value.(qcm.variables[i])
    end
    return solution
end