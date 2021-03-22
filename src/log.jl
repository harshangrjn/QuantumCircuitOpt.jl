function visualize_QCModel_solution(results::Dict{String, Any}, data::Dict{String, Any}; gate_sequence = false)

    if results["primal_status"] != MOI.FEASIBLE_POINT
        Memento.error(_LOGGER, "Non-feasible primal status. Gate decomposition may not be exact!")
    end

    R_gates_ids = findall(x -> startswith(x, "R"), data["elementary_gates"])
    U_gates_ids = findall(x -> startswith(x, "U"), data["elementary_gates"])

    gates_sol, gates_sol_compressed = get_postprocessed_solutions(results, data)

    if !isempty(gates_sol_compressed)

        printstyled("\n","=============================================================================","\n"; color = :cyan)
        
        printstyled("Quantum Circuit Model Data","\n"; color = :red)
        
        printstyled("\n","  ","Number of qubits: ", data["n_qubits"], "\n"; color = :cyan)
        
        if isempty(R_gates_ids) && isempty(U_gates_ids)
            printstyled("  ","Total number of elementary gates: ",size(data["M_real"])[3],"\n"; color = :cyan)
        else
            printstyled("  ","Total number of elementary gates (including discretization): ",size(data["M_real"])[3],"\n"; color = :cyan)
        end
        
        printstyled("  ","Maximum depth of decomposition: ", data["depth"],"\n"; color = :cyan)
        
        printstyled("  ","Input elementary gates: ", data["elementary_gates"],"\n"; color = :cyan)

        if !isempty(R_gates_ids) 
            for i in R_gates_ids
                if data["elementary_gates"][i] == "R_x"
                    printstyled("  ","R_x gate discretization: ", ceil.(rad2deg.(data["discretization"]["R_x"]), digits = 1),"\n"; color = :cyan)
                elseif data["elementary_gates"][i] == "R_y"
                    printstyled("  ","R_y gate discretization: ", ceil.(rad2deg.(data["discretization"]["R_y"]), digits = 1),"\n"; color = :cyan)
                elseif data["elementary_gates"][i] == "R_z"
                    printstyled("  ","R_z gate discretization: ", ceil.(rad2deg.(data["discretization"]["R_z"]), digits = 1),"\n"; color = :cyan)
                end
            end
        end

        if !isempty(U_gates_ids)
            for i in U_gates_ids
                if data["elementary_gates"][i] == "U3"
                    printstyled("  ","U3 gate - θ discretization: ", ceil.(rad2deg.(data["discretization"]["U3_θ"]), digits = 1),"\n"; color = :cyan)
                    printstyled("  ","U3 gate - ϕ discretization: ", ceil.(rad2deg.(data["discretization"]["U3_ϕ"]), digits = 1),"\n"; color = :cyan)
                    printstyled("  ","U3 gate - λ discretization: ", ceil.(rad2deg.(data["discretization"]["U3_λ"]), digits = 1),"\n"; color = :cyan)
                end
            end
        end
        
        printstyled("  ","Input target gate: ", data["target_gate"],"\n"; color = :cyan)
        
        printstyled("  ","Type of decomposition: ", data["decomposition_type"],"\n"; color = :cyan)

        printstyled("\n","Optimal Decomposition","\n","\n"; color = :red)
        
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

        gate_id = data["M_complex_dict"]["$id"]

        if gate_id["type"] != "Identity"
            
            s1 = gate_id["type"]

            if !(startswith(s1, "R") || startswith(s1, "U"))
            
                push!(gates_sol, s1)
            
            else
                
                s2 = String[]
                
                if gate_id["qubit_location"] == "qubit_1"
                    s2 = "1"
                elseif gate_id["qubit_location"] == "qubit_2"
                    s2 = "2"
                end

                if startswith(s1, "R")
                    θ = rad2deg(gate_id["angle"])
                    s3 = "$(θ)"
                    push!(gates_sol, string(s1," (",s2,", ",s3,")"))

                elseif startswith(s1, "U")
                    
                    θ = rad2deg(gate_id["θ"])
                    ϕ = rad2deg(gate_id["ϕ"])
                    λ = rad2deg(gate_id["λ"])
                    s3 = string("(","$(θ)",",","$(ϕ)", ",","$(λ)",")")
                    push!(gates_sol, string(s1," (",s2,", ",s3,")"))

                end
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

    else 
        return gates_sol
    end

    if isempty(gates_sol_compressed)
        Memento.error(_LOGGER, "Compressed gates solution is empty")
    end

    return gates_sol_compressed
end