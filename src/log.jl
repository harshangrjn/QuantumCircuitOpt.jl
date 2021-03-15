function visualize_QCModel_solution(results::Dict, data::Dict; gate_sequence = false)
    if results["primal_status"] != MOI.FEASIBLE_POINT
        Memento.warn(_LOGGER, "Non-feasible primal status. Gate decomposition may not be exact!")
    end

    gates_sol = Array{String,1}()
    for d = 1:data["depth"]
        id = findall(isone.(round.(results["solution"]["z_onoff_var"][:,d], digits=3)))[1]
        if data["elementary_gates"][id] != "Identity"
            push!(gates_sol, data["elementary_gates"][id])
        end
    end
    
    if !isempty(gates_sol)
        printstyled("\n","========================================================================","\n"; color = :cyan)
        printstyled("Optimal decomposition of the target gate using elementary gates:","\n"; color = :cyan)
        for i=1:length(gates_sol)
            if i != length(gates_sol)
                printstyled(gates_sol[i], " * "; color = :cyan)
            else
                printstyled(gates_sol[i], " = ", data["target_gate"]; color = :cyan)
            end
        end
        printstyled("\n","========================================================================","\n"; color = :cyan)
    end
    if gate_sequence
        return gates_sol
    end
end