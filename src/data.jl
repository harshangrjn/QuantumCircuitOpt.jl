import LinearAlgebra: I

function get_data(params::Dict{String, Any})
    if isempty(params["elementary_gates"])
        Memento.error(_LOGGER, "Input elementary gates are empty. Enter at least two unique unitary gates")
    end

    # Initial gate
    if "initial_gate" in keys(params)
        initial_gate = params["initial_gate"]
    else
        initial_gate = "Identity"
    end

    if initial_gate == "Identity"
        M_initial = Matrix(I, 2^(params["num_qubits"]+1), 2^(params["num_qubits"]+1))
    else
        Memento.error(_LOGGER, "Currently, non-identity gate is not supported as the initial condition")
    end
    # Add code here to support non-identity as an initial condition gate. 

    # Depth
    if params["depth"] < 2 
        Memento.error(_LOGGER, "Minimum depth of 2 is necessary")
    end

    # Decomposition type 
    if "decomposition_type" in keys(params)
        decomposition_type = params["decomposition_type"]
    else
        decomposition_type = "exact"
    end
    
    # Optimizer
    if "optimizer" in keys(params)
        optimizer = params["optimizer"]
    else
        Memento.error(_LOGGER, "Input a valid MIP optimizer")
    end
    
    # Relax Integrality 
    if "relax_integrality" in keys(params)
        relax_integrality = params["relax_integrality"]
    else
        # default value
        relax_integrality = false
    end

    if "slack_penalty" in keys(params)
        slack_penalty = params["relax_integrality"]
    else
        # default value
        slack_penalty = 1E3
    end

    elementary_gates = unique(params["elementary_gates"])
    
    if length(elementary_gates) < length(params["elementary_gates"])
        Memento.warn(_LOGGER, "Eliminating non-unique gates in the input elementary gates")
    end

    gates_dict, target_real = get_quantum_gates(params, elementary_gates)

    gates_dict_unique, M_real_unique, identity_idx, cnot_idx = eliminate_nonunique_gates(gates_dict)
    
    data = Dict{String, Any}("num_qubits" => params["num_qubits"],
                             "depth" => params["depth"],
                             "gates_dict" => gates_dict_unique,
                             "gates_real" => M_real_unique,
                             "initial_gate" => M_initial,
                             "identity_idx" => identity_idx,
                             "cnot_idx" => cnot_idx,
                             "elementary_gates" => elementary_gates,
                             "target_gate" => target_real,
                             "objective" => params["objective"],
                             "slack_penalty" => slack_penalty,
                             "decomposition_type" => decomposition_type,                         
                             "relax_integrality" => relax_integrality
                             )

    R_gates_ids = findall(x -> startswith(x, "R"), data["elementary_gates"])
    U_gates_ids = findall(x -> startswith(x, "U"), data["elementary_gates"])
    
    if !isempty(R_gates_ids) || !isempty(U_gates_ids)
        data["discretization"] = Dict{String, Any}()
    end

    if !isempty(R_gates_ids) 
        for i in R_gates_ids
            if data["elementary_gates"][i] == "R_x"
                data["discretization"]["R_x"] = params["R_x_discretization"]
            elseif data["elementary_gates"][i] == "R_y"
                data["discretization"]["R_y"] = params["R_y_discretization"]
            elseif data["elementary_gates"][i] == "R_z"
                data["discretization"]["R_z"] = params["R_z_discretization"]
            end
        end
    end

    if !isempty(U_gates_ids)
        for i in U_gates_ids
            if data["elementary_gates"][i] == "U3"
                data["discretization"]["U3_θ"] = params["U_θ_discretization"]
                data["discretization"]["U3_ϕ"] = params["U_ϕ_discretization"]
                data["discretization"]["U3_λ"] = params["U_λ_discretization"]
            end
        end
    end
                         
    return data
end

"""
    eliminate_nonunique_gates(gates_dict::Dict{String, Any})

"""
function eliminate_nonunique_gates(gates_dict::Dict{String, Any}; eliminate_gates = false)

    num_gates = length(keys(gates_dict))

    M_real = zeros(2*size(gates_dict["1"]["matrix"])[1], 2*size(gates_dict["1"]["matrix"])[2], num_gates)
    
    for i=1:num_gates
        M_real[:,:,i] = complex_to_real_matrix(gates_dict["$i"]["matrix"])
    end

    M_real_unique = M_real
    M_real_idx = collect(1:size(M_real)[3]) 

    if eliminate_gates
        M_real_unique, M_real_idx = QCO.unique_matrices(M_real)
    end
    
    gates_dict_unique = Dict{String, Any}()

    if size(M_real_unique)[3] < size(M_real)[3]

        Memento.info(_LOGGER, "Eliminating $(size(M_real)[3]-size(M_real_unique)[3]) non-unique gates (after discretization)")

        for i = 1:length(M_real_idx)
            gates_dict_unique["$i"] = gates_dict["$(M_real_idx[i])"]
        end
    
    else
        gates_dict_unique = gates_dict
    end

    identity_idx = _get_identity_idx(M_real_unique)

    if length(identity_idx) > 0
        if !("Identity" in gates_dict_unique["$(identity_idx[1])"]["type"])
            push!(gates_dict_unique["$(identity_idx[1])"]["type"], "Identity")
        end
    end

    cnot_idx = _get_cnot_idx(gates_dict_unique)

    return gates_dict_unique, M_real_unique, identity_idx, cnot_idx

end

function _get_identity_idx(M::Array{Float64,3})
    
    identity_idx = Int64[]
    
    for i=1:size(M)[3]
        if isapprox(M[:,:,i], Matrix(I, size(M)[1], size(M)[2]), atol=1E-5)   
            push!(identity_idx, i)
        end
    end

    # if length(identity_idx) >= 2
    #     Memento.error(_LOGGER, "Detected more than one Identity matrix in the input set of gates (possibly due to discretization)")
    # end

    return identity_idx
end

function _get_cnot_idx(gates_dict::Dict{String, Any})
    
    cnot_idx = Int64[]

    # Note: The below objective minimizes both cnot_12 and cnot_21 in the decomposition
    for i in keys(gates_dict)
        if !isempty(findall(x -> startswith(x, "cnot"), gates_dict[i]["type"]))
            push!(cnot_idx, parse(Int64, i))
        end
    end
    
    return cnot_idx
end

"""
    get_quantum_gates(params::Dict{String, Any}, elementary_gates::Array{String,1})

Given a vector of input with the names of gates (see examples folder), `get_quantum_gates` function 
returns the corresponding elementary gates in the three-dimensional complex matrix form. 
""" 
function get_quantum_gates(params::Dict{String, Any}, elementary_gates::Array{String,1})

    num_qubits = params["num_qubits"]

    gates_dict = get_all_gates_dictionary(params, elementary_gates)

    if !("target_gate" in keys(params)) || isempty(params["target_gate"])
        Memento.error(_LOGGER, "Target gate not found in the input data")
    end 

    if (size(params["target_gate"])[1] != size(params["target_gate"])[2]) || (size(params["target_gate"])[1] != 2^params["num_qubits"])
        Memento.error(_LOGGER, "Dimensions of target gate are incorrect")
    end
 
    return gates_dict, complex_to_real_matrix(params["target_gate"])
end

function get_all_gates_dictionary(params::Dict{String, Any}, elementary_gates::Array{String,1})

    num_qubits = params["num_qubits"]

    R_gates_ids = findall(x -> startswith(x, "R"), elementary_gates)
    U_gates_ids = findall(x -> startswith(x, "U"), elementary_gates)

    R_complex_dict = Dict{}
    if !isempty(R_gates_ids)
        R_complex_dict = get_all_R_gates(params, elementary_gates)
    end
    
    U_complex_dict = Dict{}
    if !isempty(U_gates_ids)
        U_complex_dict = get_all_U_gates(params, elementary_gates)
    end

    gates_dict = Dict{String, Any}()

    counter = 1

    for i=1:length(elementary_gates)

        if startswith(elementary_gates[i], "R") || startswith(elementary_gates[i], "U")
            M_elementary_dict = Dict{}

            if startswith(elementary_gates[i], "R")
                M_elementary_dict = R_complex_dict
            elseif startswith(elementary_gates[i], "U")
                M_elementary_dict = U_complex_dict
            end

            for j in keys(M_elementary_dict)
                if (j == elementary_gates[i])
                    for k in keys(M_elementary_dict[j])
                        for l in keys(M_elementary_dict[j][k]["2qubit_rep"])
                            
                            gates_dict["$counter"] = Dict{String, Any}("type" => [j],
                                                                       "angle" => Any,
                                                                       "qubit_location" => l,
                                                                       "matrix" => M_elementary_dict[j][k]["2qubit_rep"][l])

                            if startswith(elementary_gates[i], "R")
                                gates_dict["$counter"]["angle"] = M_elementary_dict[j][k]["angle"]

                            elseif startswith(elementary_gates[i], "U")
                                gates_dict["$counter"]["angle"] = Dict{String, Any}("θ" => U_complex_dict[j][k]["θ"],
                                                                                    "ϕ" => U_complex_dict[j][k]["ϕ"],
                                                                                    "λ" => U_complex_dict[j][k]["λ"],)
                            end

                            counter += 1

                        end
                    end
                end
            end
            
        else 

            gates_dict["$counter"] = Dict{String, Any}("type" => [elementary_gates[i]],
                                                       "matrix" => get_full_sized_gate(elementary_gates[i], num_qubits))
            counter += 1

        end

    end

    return gates_dict
    
end

function get_all_R_gates(params::Dict{String, Any}, elementary_gates::Array{String,1})

    R_gates_ids = findall(x -> startswith(x, "R"), elementary_gates)

    R_complex = Dict{String, Any}()

    if length(R_gates_ids) >= 1 
        for i=1:length(R_gates_ids)

            if (elementary_gates[R_gates_ids[i]] == "R_x")
                if isempty(params["R_x_discretization"])
                    Memento.error(_LOGGER, "No discretized angle was specified for R_x gate. Input at least one angle")
                end        

                R_complex["R_x"] = Dict{String, Any}()    
                R_complex["R_x"] = get_discretized_R_gates("R_x", R_complex["R_x"], collect(params["R_x_discretization"]), params["num_qubits"])
            end

            if (elementary_gates[R_gates_ids[i]] == "R_y")
                if isempty(params["R_y_discretization"])
                    Memento.error(_LOGGER, "No discretized angle was specified for R_y gate. Input at least one angle")
                end    

                R_complex["R_y"] = Dict{String, Any}()             
                R_complex["R_y"] = get_discretized_R_gates("R_y", R_complex["R_y"], collect(params["R_y_discretization"]), params["num_qubits"])
            end

            if (elementary_gates[R_gates_ids[i]] == "R_z")
                if isempty(params["R_z_discretization"])
                    Memento.error(_LOGGER, "No discretized angle was specified for R_z gate. Input at least one angle")
                end         

                R_complex["R_z"] = Dict{String, Any}()             
                R_complex["R_z"] = get_discretized_R_gates("R_z", R_complex["R_z"], collect(params["R_z_discretization"]), params["num_qubits"])
            end

        end
    end
    
    return R_complex    
end

function get_discretized_R_gates(R_type::String, R::Dict{String, Any}, discretization::Array{Float64,1}, num_qubits::Int64)

    if length(discretization) >= 1

        for i=1:length(discretization)
            R_discrete = get_pauli_rotation_gates(discretization[i])[R_type]
            R["angle_$i"] = Dict{String, Any}("angle" => discretization[i],
                                             "1qubit_rep" => R_discrete,
                                             "2qubit_rep" => Dict{String, Any}("qubit_1" => get_full_sized_gate(R_type, num_qubits, matrix=R_discrete, qubit_location = "qubit_1"),
                                                                               "qubit_2" => get_full_sized_gate(R_type, num_qubits, matrix=R_discrete, qubit_location = "qubit_2")
                                                                              ) 
                                            )            
        end
    end 

    return R 
end

function get_all_U_gates(params::Dict{String, Any}, elementary_gates::Array{String,1})

    U_gates_ids = findall(x -> startswith(x, "U"), elementary_gates)

    U_complex = Dict{String, Any}()

    if length(U_gates_ids) >= 1 
        for i=1:length(U_gates_ids)

            if (elementary_gates[U_gates_ids[i]] == "U3")        
                U_complex["U3"] = Dict{String, Any}()    
                U_complex["U3"] = get_discretized_U3_gates("U3", U_complex["U3"], collect(float(params["U_θ_discretization"])), collect(float(params["U_ϕ_discretization"])), collect(float(params["U_λ_discretization"])), params["num_qubits"])
            end

            # Add support for U1 and U2 universal gates here

        end
    end
    
    return U_complex    
end

function get_discretized_U3_gates(U_type::String, U::Dict{String, Any}, θ_discretization::Array{Float64,1}, ϕ_discretization::Array{Float64,1}, λ_discretization::Array{Float64,1}, num_qubits::Int64)
    
    counter = 1 

    for i=1:length(θ_discretization)
        for j=1:length(ϕ_discretization)
            for k=1:length(λ_discretization)
                
                U_discrete = QCO.U3Gate(θ_discretization[i], ϕ_discretization[j], λ_discretization[k])
                
                U["angle_$(counter)"] = Dict{String, Any}("θ" => θ_discretization[i],
                                                          "ϕ" => ϕ_discretization[j],
                                                          "λ" => λ_discretization[k],
                                                          "1qubit_rep" => U_discrete,
                                                          "2qubit_rep" => Dict{String, Any}("qubit_1" => get_full_sized_gate(U_type, num_qubits, matrix = U_discrete, qubit_location = "qubit_1"),
                                                                                            "qubit_2" => get_full_sized_gate(U_type, num_qubits, matrix = U_discrete, qubit_location = "qubit_2")
                                                                                           ) 
                                                         )
                counter += 1
            end
        end
    end
    
    return U
end

"""
    get_full_sized_gate(input::String, num_qubits::Int64; matrix = nothing, qubit_location = nothing)

For a given string and number of qubits in the input specified input, this function returns a full 
sized gate with respect to the input number of qubits. 
"""
function get_full_sized_gate(input::String, num_qubits::Int64; matrix = nothing, qubit_location = nothing)

    if input == "Identity"
        return QCO.IGate(num_qubits)
    end
    
    # All 2-qubit full-sized gates
    if num_qubits == 2

        if input == "H1"
            return kron(QCO.HGate(), QCO.IGate(1))

        elseif input == "H2"
            return kron(QCO.IGate(1), QCO.HGate())

        elseif input == "cnot_12"
            return QCO.CNotGate()

        elseif input == "cnot_21"
            return QCO.CNotRevGate()

        elseif input == "cnot_swap"
            return QCO.CNotGate() * QCO.CNotRevGate()

        elseif input == "H1⊗H2"
            return kron(QCO.HGate(), QCO.HGate())

        elseif input == "T1"
            return kron(QCO.TGate(), QCO.IGate(1)) 

        elseif input == "T2"
            return kron(QCO.IGate(1), QCO.TGate())   

        elseif input == "T1_dagger"
            return kron(QCO.TdaggerGate(), QCO.IGate(1)) 

        elseif input == "T2_dagger"
            return kron(QCO.IGate(1), QCO.TdaggerGate())   

        elseif input == "S1"
            return kron(QCO.SGate(), QCO.IGate(1)) 

        elseif input == "S2"
            return kron(QCO.IGate(1), QCO.SGate())  

        elseif input == "X1"
            return kron(QCO.XGate(), QCO.IGate(1)) 

        elseif input == "X2"
            return kron(QCO.IGate(1), QCO.XGate())   

        elseif input == "Y1"
            return kron(QCO.YGate(), QCO.IGate(1)) 

        elseif input == "Y2"
            return kron(QCO.IGate(1), QCO.YGate())   

        elseif input == "Z1"
            return kron(QCO.ZGate(), QCO.IGate(1)) 

        elseif input == "Z2"
            return kron(QCO.IGate(1), QCO.ZGate())        
            
        elseif input == "CZ_12"
            return QCO.CZGate()

        elseif input == "CH_12"
            return QCO.CHGate()

        elseif input == "CV_12"
            return QCO.CVGate()

        elseif input == "swap"
            return QCO.SwapGate()

        elseif input == "M_12"
            return QCO.MGate()

        elseif input == "QFT_12"
            return QCO.QFT2Gate()

        elseif input == "CSX_12"
            return QCO.C2SXGate()
        
        elseif input == "W_12"
            return QCO.WGate()
        
        # Gates with continuous angle parameters
        elseif input in ["R_x", "R_y", "R_z", "U3"] 

            if qubit_location == "qubit_1"
                return kron(matrix, QCO.IGate(1))
            elseif qubit_location == "qubit_2"
                return kron(QCO.IGate(1), matrix)
            end

        else
            
            Memento.error(_LOGGER, "Specified input elementary gates or the target gate does not exist in the predefined set of gates.")
        end
    end

    # All 3-qubit full-sized gates
    if num_qubits == 3

        if input == "toffoli"
            return QCO.ToffoliGate()

        elseif input == "CSwap"
            return QCO.CSwapGate()

        elseif input == "CCZ"
            return QCO.CCZGate()

        elseif input == "peres"
            return QCO.PeresGate()

        elseif input == "cnot_13"
            # |0⟩⟨0| ⊗ I ⊗ I 
            control_0 = kron(Array{Complex{Float64},2}([1 0; 0 0]) , kron(QCO.IGate(1), QCO.IGate(1)))
            # |1⟩⟨1| ⊗ I ⊗ X 
            control_1 = kron(Array{Complex{Float64},2}([0 0; 0 1]) , kron(QCO.IGate(1), QCO.XGate()))
            return control_0 + control_1

        elseif input == "cnot_31"
            # I ⊗ I ⊗ |0⟩⟨0| 
            control_0 = kron(QCO.IGate(1), kron(QCO.IGate(1), Array{Complex{Float64},2}([1 0; 0 0])))
            # X ⊗ I ⊗ |1⟩⟨1|
            control_1 = kron(kron(QCO.XGate(), QCO.IGate(1)), Array{Complex{Float64},2}([0 0; 0 1]))
            return control_0 + control_1

        else
            Memento.error(_LOGGER, "Specified input elementary gates or the target gate does not exist in the predefined set of gates.")
        end
    
    end

end

function num_input_gates(params::Dict{String, Any}, elementary_gates::Array{String,1})

    # Update this if the naming convention for gates changes
    R_gates = findall(x -> startswith(x, "R"), (elementary_gates))
    U_gates = findall(x -> startswith(x, "U"), (elementary_gates))

    num_gates = 0
    
    if isempty(R_gates) && isempty(U_gates)

        num_gates = length(elementary_gates)

    elseif !isempty(R_gates) && isempty(U_gates)
        
        n_discretized_R_gates = get_number_of_R_gates(params, elementary_gates)

        num_gates = length(elementary_gates) - length(R_gates) + n_discretized_R_gates

    elseif !isempty(U_gates) && isempty(R_gates)
        
        n_discretized_U_gates = _num_U3_gates(params, elementary_gates)

        num_gates = length(elementary_gates) - length(U_gates) + n_discretized_U_gates

    elseif !isempty(R_gates) && !isempty(U_gates)
        
        n_discretized_R_gates = get_number_of_R_gates(params, elementary_gates)

        n_discretized_U_gates = _num_U3_gates(params, elementary_gates)

        num_gates = length(elementary_gates) - length(R_gates) - length(U_gates) + n_discretized_R_gates + n_discretized_U_gates

    end
    
    if num_gates == 1
        Memento.error(_LOGGER, "Input at least two unique elementary gates for non-trivial solutions")
    else
        return num_gates
    end

end

function get_number_of_R_gates(params::Dict{String, Any}, elementary_gates::Array{String,1})
    
    num_R_x = 0; num_R_y = 0; num_R_z = 0;
    
    # A factor of two to account for R gates on both qubits
    if ("R_x" in elementary_gates)
        if !("R_x_discretization" in keys(params))
            Memento.error(_LOGGER, "Discretization angles for the R_x gate are unspecified")
        end

        num_R_x = 2*length(params["R_x_discretization"])
        if num_R_x == 0
            Memento.error(_LOGGER, "Dicsretization not specified for R_x gate")
        end
    end
    
    if ("R_y" in elementary_gates)
        if !("R_y_discretization" in keys(params))
            Memento.error(_LOGGER, "Discretization angles for the R_y gate are unspecified")
        end

        num_R_y = 2*length(params["R_y_discretization"])
        if num_R_y == 0 
            Memento.error(_LOGGER, "Dicsretization not specified for R_y gate")
        end 
    end

    if ("R_z" in elementary_gates)
        if !("R_z_discretization" in keys(params))
            Memento.error(_LOGGER, "Discretization angles for the R_z gate are unspecified")
        end

        num_R_z = 2*length(params["R_z_discretization"])
        if num_R_z == 0
            Memento.error(_LOGGER, "Dicsretization not specified for R_z gate")
        end
    end

    return (num_R_x + num_R_y + num_R_z)
end

function _num_U3_gates(params::Dict{String, Any}, elementary_gates::Array{String,1})

    num_U3 = 0

        if ("U3" in elementary_gates)    
            
            if !("U_θ_discretization" in keys(params)) || !("U_ϕ_discretization" in keys(params)) || !("U_λ_discretization" in keys(params))
                Memento.error(_LOGGER, "Discretization angles for the U3 gate are unspecified")
            end

            if (isempty(params["U_θ_discretization"])) || (isempty(params["U_ϕ_discretization"])) || (isempty(params["U_λ_discretization"]))
                
                Memento.error(_LOGGER, "Input atleast one discretization angle for every angle (θ,ϕ,λ)")

            else
                # A factor of two to account for U gates on both qubits                
                num_U3 = 2 * length(params["U_θ_discretization"]) * length(params["U_ϕ_discretization"]) * length(params["U_λ_discretization"])

            end
        end

    return num_U3
end

"""
    get_pauli_rotation_gates(θ::Number)

For a given angle in radiaons, this function returns standard rotation gates those define rotations around the Pauli axis {X,Y,Z}.
Note that R_x(θ) = u3(θ, -π/2, π/2), R_y(θ) = u3(θ, 0, 0), R_z(λ) = exp((-λ/2)im)*u1(λ). 
"""
function get_pauli_rotation_gates(θ::Number)
    #input angles in radians
    if !(-2*π <= θ <= 2*π)
        Memento.error(_LOGGER, "θ angle in Pauli rotation gate is not within valid bounds")
    end

    R_x = Array{Complex{Float64},2}([cos(θ/2) -(sin(θ/2))im; -(sin(θ/2))im cos(θ/2)])
    R_y = Array{Complex{Float64},2}([cos(θ/2) -(sin(θ/2)); (sin(θ/2)) cos(θ/2)])
    R_z = Array{Complex{Float64},2}([(cos(θ/2) - (sin(θ/2))im) 0; 0 (cos(θ/2) + (sin(θ/2))im)])
    pauli_rotation_gates = Dict{String, Any}("R_x" => round_complex_values(R_x), 
                                             "R_y" => round_complex_values(R_y), 
                                             "R_z" => round_complex_values(R_z)
                                            )

    return pauli_rotation_gates
end