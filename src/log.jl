function visualize_QCModel_solution(results::Dict, data::Dict; gate_sequence = false)
    if results["primal_status"] != MOI.FEASIBLE_POINT
        Memento.error(_LOGGER, "Non-feasible primal status. Gate decomposition may not be exact!")
    end

    gates_sol = get_postprocessed_solutions(results, data)

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
    else
        Memento.warn(_LOGGER, "Valid integral feasible solutions could not be found to visualize the solution")
    end

    if gate_sequence
        return gates_sol
    end
end

function get_postprocessed_solutions(results::Dict, data::Dict)
    gates_sol = Array{String,1}()
    id_sequence = Array{Int64,1}()

    for d = 1:data["depth"]
        id = findall(isone.(round.(abs.(results["solution"]["z_onoff_var"][:,d]), digits=3)))[1]
        push!(id_sequence, id)

        if data["M_complex_dict"]["$id"]["type"] != "Identity"
            s1 = data["M_complex_dict"]["$id"]["type"]

            if data["M_complex_dict"]["$id"]["angle"] == "na"
                push!(gates_sol, s1)
            else
                s2 = String[]
                if data["M_complex_dict"]["$id"]["qubit_location"] == "qubit_1"
                    s2 = "1"
                elseif data["M_complex_dict"]["$id"]["qubit_location"] == "qubit_2"
                    s2 = "2"
                end
                s3 = "$(rad2deg(data["M_complex_dict"]["$id"]["angle"]))"
                push!(gates_sol, string(s1,"-",s2," (",s3,")"))
            end

        end

    end
    
    validate_solutions(data, id_sequence)

    return gates_sol  
end

function validate_solutions(data::Dict, id_sequence::Array{Int64,1})
    M = Matrix(LA.I, 2^(data["n_qubits"]+1), 2^(data["n_qubits"]+1))
    
    for i in id_sequence
        M *= QCO.get_complex_to_real_matrix(data["M_complex_dict"]["$i"]["matrix"])
    end

    if !isapprox(M, data["Target_real"])
        Memento.error(_LOGGER, "Decomposition is not valid: Problem may be infeasible")
    end
    
end