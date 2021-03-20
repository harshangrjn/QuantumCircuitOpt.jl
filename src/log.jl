function visualize_QCModel_solution(results::Dict{String, Any}, data::Dict{String, Any}; gate_sequence = false)

    if results["primal_status"] != MOI.FEASIBLE_POINT
        Memento.error(_LOGGER, "Non-feasible primal status. Gate decomposition may not be exact!")
    end

    gates_sol, gates_sol_compressed = get_postprocessed_solutions(results, data)

    if !isempty(gates_sol_compressed)

        printstyled("\n","=============================================================================","\n"; color = :cyan)
        
        printstyled("Problem statistics","\n"; color = :red)
        
        printstyled("\n","  ","Number of qubits: ", data["n_qubits"], "\n"; color = :cyan)
        
        printstyled("  ","Total number of elementary gates (including discretization): ",size(data["M_real"])[3],"\n"; color = :cyan)
        
        printstyled("  ","Maximum depth of decomposition: ", data["depth"],"\n"; color = :cyan)
        
        printstyled("  ","Input elementary gates: ", data["elementary_gates"],"\n"; color = :cyan)
        
        printstyled("  ","Input target gate: ", data["target_gate"],"\n"; color = :cyan)
        
        printstyled("  ","Type of decomposition: ", data["decomposition_type"],"\n"; color = :cyan)

        printstyled("\n","Optimal decomposition","\n","\n"; color = :red)
        
        print("  ")
        
        for i=1:length(gates_sol_compressed)
            
            if i != length(gates_sol_compressed)
                printstyled(gates_sol_compressed[i], " * "; color = :cyan)
            else
                printstyled(gates_sol_compressed[i], " = ", data["target_gate"],"\n"; color = :cyan)
            end

        end

        if data["objective"] == "minimize_depth"
            printstyled("  ","Minimum optimal depth: ", length(gates_sol_compressed),"\n"; color = :cyan)

        elseif data["objective"] == "minimize_cnot"
            printstyled("  ","Minimum number of CNOT gates: ", results["objective"],"\n"; color = :cyan)

        end

        printstyled("  ","Optimizer run time: ", ceil(results["solve_time"], digits=2)," sec.","\n"; color = :cyan)

        printstyled("=============================================================================","\n"; color = :cyan)
        
        @show gates_sol

    else
        Memento.warn(_LOGGER, "Valid integral feasible solutions could not be found to visualize the solution")
    end

    if gate_sequence
        return gates_sol
    end

end

function get_postprocessed_solutions(results::Dict{String, Any}, data::Dict{String, Any})

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
                push!(gates_sol, string(s1," (",s2,", ",s3,")"))
            end

        end

    end
    
    validate_solutions(data, id_sequence)

    gates_sol_compressed = get_compressed_solutions(data, gates_sol)

    return gates_sol, gates_sol_compressed
end

function validate_solutions(data::Dict{String, Any}, id_sequence::Array{Int64,1})
    
    M_sol = Array{Complex{Float64},2}(Matrix(LA.I, 2^(data["n_qubits"]), 2^(data["n_qubits"])))
    
    for i in id_sequence
        M_sol *= data["M_complex_dict"]["$i"]["matrix"]
    end

    # @show M_sol 
    # @show QCO.get_real_to_complex_matrix(data["Target_real"])

    # This tolerance is very important for the final feasiblity check
    if !isapprox(M_sol, QCO.get_real_to_complex_matrix(data["Target_real"]), atol = 1E-4)
        Memento.error(_LOGGER, "Decomposition is not valid: Problem may be infeasible")
    end
    
end

"""
get_compressed_solutions returns sequence of gates after compressing adjacent pair of gates represented on two separate qubits. 
For example, gates H1 and H2 appearing in a sequence will be compressed to H1⊗H2. 
"""
function get_compressed_solutions(data::Dict{String, Any}, gates_sol::Array{String,1})
    gates_sol_compressed = String[]

    if length(gates_sol) == 1 
        return gates_sol
    end

    # This part of the code is a bit ad-hoc. This needs to be updated once the input format for elementary gates gets cleaned up. 
    if isempty(findall(x -> startswith(x, "R") || startswith(x, "U"), data["elementary_gates"]))
        
        status = false

        for i=1:(length(gates_sol))
            if i <= length(gates_sol) - 1
                if status 
                    status = false
                    continue
                else
                    if !(startswith(gates_sol[i], "cnot")) && !(startswith(gates_sol[i+1], "cnot"))
                        if (occursin('1', gates_sol[i]) && occursin('2', gates_sol[i+1])) || (occursin('2', gates_sol[i]) && occursin('1', gates_sol[i+1])) 
                            if occursin('1', gates_sol[i])
                                gate_string = string(gates_sol[i],"⊗",gates_sol[i+1])
                            else 
                                gate_string = string(gates_sol[i+1],"⊗",gates_sol[i])
                            end
                            push!(gates_sol_compressed, gate_string)
                            status = true
                            continue
                        else
                            push!(gates_sol_compressed, gates_sol[i])

                        end
                    else
                        push!(gates_sol_compressed, gates_sol[i])

                    end
                end
            else
                if !status
                    push!(gates_sol_compressed, gates_sol[i])
                end
            end
        end 

    end

    if isempty(gates_sol_compressed)
        Memento.error(_LOGGER, "Compressed gates solution is empty")
    end

    return gates_sol_compressed
end