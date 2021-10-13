"""
    visualize_solution(results::Dict{String, Any}, data::Dict{String, Any}; gate_sequence = false)

Given dictionaries of results and data, and assuming that the optimization model had a feasible solution, 
this function aids in visualizing the optimal circuit decomposition.
"""
function visualize_solution(results::Dict{String, Any}, data::Dict{String, Any}; gate_sequence = false)

    if data["relax_integrality"]
        if results["primal_status"] == MOI.FEASIBLE_POINT 
            Memento.info(_LOGGER, "Integrality-relaxed solutions can be found in the results dictionary")
        else
            Memento.info(_LOGGER, "Infeasible primal status for the integrality-relaxed problem")
        end

        return
    end

    if results["primal_status"] != MOI.FEASIBLE_POINT 
        
        if results["termination_status"] != MOI.TIME_LIMIT
            Memento.warn(_LOGGER, "Infeasible primal status. Gate decomposition may be inaccurate")
        else 
            Memento.warn(_LOGGER, "Optimizer hits time limit with an infeasible primal status. Gate decomposition may be inaccurate")
        end

        return
    else
        gates_sol, gates_sol_compressed = QCO.get_postprocessed_decomposition(results, data)
    end

    if !isempty(gates_sol_compressed)

        printstyled("\n","=============================================================================","\n"; color = :cyan)
        
        printstyled("Quantum Circuit Model Data","\n"; color = :red)
        
        printstyled("\n","  ","Number of qubits: ", data["num_qubits"], "\n"; color = :cyan)
        
        printstyled("  ","Total number of elementary gates (after presolve): ",size(data["gates_real"])[3],"\n"; color = :cyan)
        
        printstyled("  ","Maximum depth of decomposition: ", data["maximum_depth"],"\n"; color = :cyan)
        
        printstyled("  ","Input elementary gates: ", data["elementary_gates"],"\n"; color = :cyan)

        if "discretization" in keys(data)
            for i in keys(data["discretization"])
                printstyled("    ","$i discretization: ", ceil.(rad2deg.(data["discretization"][i]), digits = 1),"\n"; color = :cyan)
            end
        end
        
        # printstyled("  ","Input target gate: ", data["target_gate"]["type"],"\n"; color = :cyan)
        
        printstyled("  ","Type of decomposition: ", data["decomposition_type"],"\n"; color = :cyan)

        printstyled("\n","Optimal Circuit Decomposition","\n","\n"; color = :red)
        
        print("  ")
        
        for i=1:length(gates_sol_compressed)
            
            if i != length(gates_sol_compressed)
                printstyled(gates_sol_compressed[i], " * "; color = :cyan)
            else    
                if data["decomposition_type"] == "exact"
                    printstyled(gates_sol_compressed[i], " = ", "Target gate","\n"; color = :cyan)
                elseif data["decomposition_type"] == "approximate"
                    printstyled(gates_sol_compressed[i], " ≈ ", "Target gate","\n"; color = :cyan)
                end
            end

        end

        if data["decomposition_type"] == "approximate"
            printstyled("  ","||Decomposition error||₂: ", round(LA.norm(results["solution"]["slack_var"]), digits = 10),"\n"; color = :cyan)
        end

        if data["objective"] == "minimize_depth"

            if length(data["identity_idx"]) >= 1 && !(results["termination_status"] == MOI.TIME_LIMIT)
                printstyled("  ","Minimum optimal depth: ", length(gates_sol_compressed),"\n"; color = :cyan)
            else 
                printstyled("  ","Decomposition depth: ", length(gates_sol_compressed),"\n"; color = :cyan)
            end

        elseif data["objective"] == "minimize_cnot"

            if !isempty(data["cnot_idx"])
                
                if data["decomposition_type"] == "exact" 
                    printstyled("  ","Minimum number of CNOT gates: ", round(results["objective"], digits = 6),"\n"; color = :cyan)
                
                elseif data["decomposition_type"] == "approximate"
                    printstyled("  ","Minimum number of CNOT gates: ", round((results["objective"] - data["slack_penalty"]*LA.norm(results["solution"]["slack_var"])^2), digits = 6),"\n"; color = :cyan)
                end
            
            elseif !isempty(data["identity_idx"])
               
                printstyled("  ","Objective function changed to minimize_depth (CNOT gate not found in input)", "\n"; color = :cyan) 
                
                if !(results["termination_status"] == MOI.TIME_LIMIT)
                    printstyled("  ","Minimum optimal depth: ", length(gates_sol_compressed),"\n"; color = :cyan)
                else
                    printstyled("  ","Decomposition depth: ", length(gates_sol_compressed),"\n"; color = :cyan)
                end
            
            else

                printstyled("  ","Objective function changed to feasiblity (CNOT and Identity gates not found in input)", "\n"; color = :cyan) 
                printstyled("  ","Compressed depth for a feasbile decomposition: ", length(gates_sol_compressed),"\n"; color = :cyan)

            end

        end

        printstyled("  ","Optimizer run time: ", ceil(results["solve_time"], digits=2)," sec.","\n"; color = :cyan)
            
        if results["termination_status"] == MOI.TIME_LIMIT
            printstyled("  ","Termination status: TIME_LIMIT", "\n"; color = :cyan)
            printstyled("  ","Decomposition may not be optimal", "\n"; color = :cyan)
        end

        printstyled("=============================================================================","\n"; color = :cyan)      

    else
        Memento.warn(_LOGGER, "Valid integral feasible solutions could not be found to visualize the solution")
    end

    if gate_sequence
        return gates_sol
    end

end

function get_postprocessed_decomposition(results::Dict{String, Any}, data::Dict{String, Any})

    gates_sol = Array{String,1}()
    id_sequence = Array{Int64,1}()

    for d = 1:data["maximum_depth"]
        id = findall(isone.(round.(abs.(results["solution"]["z_onoff_var"][:,d]), digits=3)))[1]
        push!(id_sequence, id)

        gate_id = data["gates_dict"]["$id"]

        if !("Identity" in gate_id["type"])
            
            s1 = gate_id["type"][1]

            if occursin(kron_symbol, s1)
                push!(gates_sol, s1) 

            elseif !(QCO._parse_gate_string(s1, type = true) in union(QCO.ONE_QUBIT_GATES_ANGLE_PARAMETERS, QCO.TWO_QUBIT_GATES_ANGLE_PARAMETERS))
                push!(gates_sol, s1) 

            else
                
                s2 = String[]          
                for i_qu = 1:data["num_qubits"]
                    if gate_id["qubit_loc"] == "qubit_$i_qu"    
                        s2 = "$i_qu"
                    end
                end

                if "angle" in keys(gate_id)

                    if length(keys(gate_id["angle"])) == 1 
                        θ = round(rad2deg(gate_id["angle"]), digits = 3)
                        s3 = "$(θ)"
                        push!(gates_sol, string(s1,"(", s3, ")"))
                    elseif length(keys(gate_id["angle"])) == 3
                        θ = round(rad2deg(gate_id["angle"]["θ"]), digits = 3)
                        ϕ = round(rad2deg(gate_id["angle"]["ϕ"]), digits = 3)
                        λ = round(rad2deg(gate_id["angle"]["λ"]), digits = 3)
                        s3 = string("(","$(θ)",",","$(ϕ)", ",","$(λ)",")")
                        push!(gates_sol, string(s1, s3))
                    end

                end
            end
        end
    end
    
    QCO.validate_circuit_decomposition(data, id_sequence)

    gates_sol_compressed = QCO.get_compressed_decomposition(data["num_qubits"], gates_sol)

    return gates_sol, gates_sol_compressed
end

"""
    validate_circuit_decomposition(data::Dict{String, Any}, id_sequence::Array{Int64,1})

This function validates the circuit decomposition if it is indeed exact with respect to the specified target gate. 
"""
function validate_circuit_decomposition(data::Dict{String, Any}, id_sequence::Array{Int64,1})
    
    M_sol = Array{Complex{Float64},2}(Matrix(LA.I, 2^(data["num_qubits"]), 2^(data["num_qubits"])))
    
    for i in id_sequence
        M_sol *= data["gates_dict"]["$i"]["matrix"]
    end

    # This tolerance is very important for the final feasiblity check
    if data["are_gates_real"]
        target_gate = real(data["target_gate"])
    else
        target_gate = QCO.real_to_complex_gate(data["target_gate"])
    end
    
    if (data["decomposition_type"] == "exact") && (!isapprox(M_sol, target_gate, atol = 1E-4))
        Memento.error(_LOGGER, "Decomposition is not valid: Problem may be infeasible")
    end
    
end

"""
    get_compressed_decomposition(num_qubits::Int64, gates_sol::Array{String,1})

Given the number of qubits and the sequence of gates from the solution, this function returns a 
decomposition of gates after compressing adjacent pair of gates represented on two separate qubits. 
For example, gates H1 and H2 appearing in a sequence will be compressed to H1xH2 (kron(H1,H2)). 
This functionality is currently supported only for two qubit circuits and gates without angle parameters. 
"""
function get_compressed_decomposition(num_qubits::Int64, gates_sol::Array{String,1})
    # This part of the code may be hacky. This needs to be updated once the input format gets cleaned up for elementary gates with U and R gates.     

    if (length(gates_sol) == 1) || (num_qubits > 2)
        return gates_sol
    end
    
    gates_sol_compressed = String[]

    angle_param_gate = false
    for i=1:length(gates_sol)
        if !occursin(kron_symbol, gates_sol[i])
            gates_sol_type = QCO._parse_gate_string(gates_sol[i], type = true)
            
            if gates_sol_type in union(QCO.ONE_QUBIT_GATES_ANGLE_PARAMETERS, QCO.TWO_QUBIT_GATES_ANGLE_PARAMETERS)
                angle_param_gate = true
                break
            end
        end       
    end

    if !angle_param_gate
        
        status = false

        for i=1:(length(gates_sol))
            if i <= length(gates_sol) - 1
                if status 
                    status = false
                    continue
                else
                    gate_i = QCO.is_multi_qubit_gate(gates_sol[i])
                    gate_iplus1 = QCO.is_multi_qubit_gate(gates_sol[i+1])

                    if !(gate_i) && !(gate_iplus1)
                        if (occursin('1', gates_sol[i]) && occursin('2', gates_sol[i+1])) || (occursin('2', gates_sol[i]) && occursin('1', gates_sol[i+1])) 
                            if occursin('1', gates_sol[i])
                                gate_string = string(gates_sol[i],"x",gates_sol[i+1])
                            else 
                                gate_string = string(gates_sol[i+1],"x",gates_sol[i])
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
