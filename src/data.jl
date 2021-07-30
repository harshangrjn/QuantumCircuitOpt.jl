import LinearAlgebra: I

function get_data(params::Dict{String, Any}; eliminate_identical_gates = true)
    
    # Number of qubits
    if "num_qubits" in keys(params)
        if params["num_qubits"] < 2 
            Memento.error(_LOGGER, "Minimum of 2 qubits is necessary")
        end
        num_qubits = params["num_qubits"]
    else
        Memento.error(_LOGGER, "Number of qubits has to be specified by the user")
    end

    # Depth
    if "depth" in keys(params)
        if params["depth"] < 2 
            Memento.error(_LOGGER, "Minimum depth of 2 is necessary")
        end
        depth = params["depth"]
    else
        Memento.error(_LOGGER, "Depth of decomposition has to be specified by the user")
    end

    # Elementary gates
    if !("elementary_gates" in keys(params)) || isempty(params["elementary_gates"])
        Memento.error(_LOGGER, "Input elementary gates is empty. Enter at least two unique unitary gates")
    end

    # Initial gate
    if "initial_gate" in keys(params)
        if params["initial_gate"] == "Identity"
            initial_gate = QCO.complex_to_real_matrix(QCO.IGate(num_qubits))
        else 
            Memento.error(_LOGGER, "Currently, only \"Identity\" is supported as an initial gate")
            # Add code here to support non-identity as an initial gate. 
        end
    else
        initial_gate = QCO.complex_to_real_matrix(QCO.IGate(num_qubits))
    end

    # Input Circuit
    if "input_circuit" in keys(params)
        input_circuit = params["input_circuit"]
    else
        # default value
        input_circuit = []
    end
    
    input_circuit_dict = Dict{String,Any}()

    if length(input_circuit) > 0 && (length(input_circuit) <= params["depth"])
        
        input_circuit_dict = QCO.get_input_circuit_dict(input_circuit, params)

    else
        (length(input_circuit) > 0) && (Memento.warn(_LOGGER, "Neglecting the input circuit as it's depth is greater than the allowable depth"))
    end

    # Decomposition type 
    if "decomposition_type" in keys(params)
        decomposition_type = params["decomposition_type"]
    else
        decomposition_type = "exact"
    end

    # Objective function
    if "objective" in keys(params)
        objective = params["objective"]
    else
        objective = "minimize_depth"
    end
    
    # Relax Integrality 
    if "relax_integrality" in keys(params)
        relax_integrality = params["relax_integrality"]
    else
        # default value
        relax_integrality = false
    end

    # Optimizer time limit (in seconds)
    if "time_limit" in keys(params) && params["time_limit"] isa Number
        time_limit = params["time_limit"]
    else
        # default value
        time_limit = 10800
    end

    elementary_gates = unique(params["elementary_gates"])
    
    if length(elementary_gates) < length(params["elementary_gates"])
        Memento.warn(_LOGGER, "Eliminating non-unique gates in the input elementary gates")
    end

    gates_dict, target_real = QCO.get_quantum_gates(params, elementary_gates)

    gates_dict_unique, M_real_unique, identity_idx, cnot_idx = eliminate_nonunique_gates(gates_dict, eliminate_identical_gates = eliminate_identical_gates)
    
    data = Dict{String, Any}("num_qubits" => num_qubits,
                             "depth" => depth,
                             "gates_dict" => gates_dict_unique,
                             "gates_real" => M_real_unique,
                             "initial_gate" => initial_gate,
                             "identity_idx" => identity_idx,
                             "cnot_idx" => cnot_idx,
                             "elementary_gates" => elementary_gates,
                             "target_gate" => target_real,
                             "objective" => objective,
                             "decomposition_type" => decomposition_type,                         
                             "relax_integrality" => relax_integrality,
                             "time_limit" => time_limit
                             )

    R_gates_ids = findall(x -> startswith(x, "R"), data["elementary_gates"])
    U_gates_ids = findall(x -> startswith(x, "U"), data["elementary_gates"])
    
    if !isempty(R_gates_ids) || !isempty(U_gates_ids)
        data["discretization"] = Dict{String, Any}()
    end

    if !isempty(R_gates_ids) 
        for i in R_gates_ids
            if data["elementary_gates"][i] == "RX"
                data["discretization"]["RX"] = params["RX_discretization"]
            elseif data["elementary_gates"][i] == "RY"
                data["discretization"]["RY"] = params["RY_discretization"]
            elseif data["elementary_gates"][i] == "RZ"
                data["discretization"]["RZ"] = params["RZ_discretization"]
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

    # Input circuit
    if length(keys(input_circuit_dict)) > 0
        data["input_circuit"] = input_circuit_dict
    end

    # Slack Penalty
    if decomposition_type == "approximate"
        if "slack_penalty" in keys(params)
            data["slack_penalty"] = params["slack_penalty"]
        else
            # default value
            data["slack_penalty"] = 1E3
        end
    end
                         
    return data
end

"""
    eliminate_nonunique_gates(gates_dict::Dict{String, Any})

"""
function eliminate_nonunique_gates(gates_dict::Dict{String, Any}; eliminate_identical_gates = false)

    num_gates = length(keys(gates_dict))

    M_real = zeros(2*size(gates_dict["1"]["matrix"])[1], 2*size(gates_dict["1"]["matrix"])[2], num_gates)
    
    for i=1:num_gates
        M_real[:,:,i] = complex_to_real_matrix(gates_dict["$i"]["matrix"])
    end

    M_real_unique = M_real
    M_real_idx = collect(1:size(M_real)[3]) 

    if eliminate_identical_gates
        M_real_unique, M_real_idx = QCO.unique_matrices(M_real)
    end
    
    gates_dict_unique = Dict{String, Any}()

    if size(M_real_unique)[3] < size(M_real)[3]

        Memento.info(_LOGGER, "Detected $(size(M_real)[3]-size(M_real_unique)[3]) non-unique gates (after discretization)")

        for i = 1:length(M_real_idx)
            gates_dict_unique["$i"] = gates_dict["$(M_real_idx[i])"]
        end
    
    else
        gates_dict_unique = gates_dict
    end

    identity_idx = QCO._get_identity_idx(M_real_unique)

    for i_id = 1:length(identity_idx)
        if !("Identity" in gates_dict_unique["$(identity_idx[i_id])"]["type"])
            push!(gates_dict_unique["$(identity_idx[i_id])"]["type"], "Identity")
        end
    end

    cnot_idx = QCO._get_cnot_idx(gates_dict_unique)

    return gates_dict_unique, M_real_unique, identity_idx, cnot_idx

end

function _get_identity_idx(M::Array{Float64,3})
    
    identity_idx = Int64[]
    
    for i=1:size(M)[3]
        if isapprox(M[:,:,i], Matrix(I, size(M)[1], size(M)[2]), atol=1E-5)   
            push!(identity_idx, i)
        end
    end

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

    gates_dict = QCO.get_all_gates_dictionary(params, elementary_gates)

    if !("target_gate" in keys(params)) || isempty(params["target_gate"])
        Memento.error(_LOGGER, "Target gate not found in the input data")
    end 

    if (size(params["target_gate"])[1] != size(params["target_gate"])[2]) || (size(params["target_gate"])[1] != 2^params["num_qubits"])
        Memento.error(_LOGGER, "Dimensions of target gate do not match the input num_qubits")
    end
 
    return gates_dict, complex_to_real_matrix(params["target_gate"])
end

function get_all_gates_dictionary(params::Dict{String, Any}, elementary_gates::Array{String,1})

    num_qubits = params["num_qubits"]

    R_gates_ids = findall(x -> startswith(x, "R"), elementary_gates)
    U_gates_ids = findall(x -> startswith(x, "U"), elementary_gates)

    R_complex_dict = Dict{}
    if !isempty(R_gates_ids)
        R_complex_dict = QCO.get_all_R_gates(params, elementary_gates)
    end
    
    U_complex_dict = Dict{}
    if !isempty(U_gates_ids)
        U_complex_dict = QCO.get_all_U_gates(params, elementary_gates)
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
                        for l in keys(M_elementary_dict[j][k]["$(num_qubits)qubit_rep"])
                            
                            M_sqrd = M_elementary_dict[j][k]["$(num_qubits)qubit_rep"][l]^2

                            gates_dict["$counter"] = Dict{String, Any}("type" => [j],
                                                                       "angle" => Any,
                                                                       "qubit_loc" => l,
                                                                       "matrix" => M_elementary_dict[j][k]["$(num_qubits)qubit_rep"][l],
                                                                       "isInvolutory" => isapprox(M_sqrd, Matrix(LA.I, size(M_sqrd)[1], size(M_sqrd)[2]), atol=1E-6))

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
            
            M = QCO.get_full_sized_gate(elementary_gates[i], num_qubits)
            M_sqrd = M^2

            gates_dict["$counter"] = Dict{String, Any}("type" => [elementary_gates[i]],
                                                       "matrix" => M,
                                                       "isInvolutory" => isapprox(M_sqrd, Matrix(LA.I, size(M_sqrd)[1], size(M_sqrd)[2]), atol=1E-6))
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

            gate_type = elementary_gates[R_gates_ids[i]]

            if !(gate_type in ["RX", "RY", "RZ"])
                Memento.error(_LOGGER, "Input R gate type ($(gate_type)) is not supported.")
            end

            if isempty(params[string(gate_type,"_discretization")])
                Memento.error(_LOGGER, "No discretized angle was specified for $(gate_type) gate. Input at least one angle")
            end        

            R_complex["$(gate_type)"] = Dict{String, Any}()    
            R_complex["$(gate_type)"] = QCO.get_discretized_R_gates(gate_type, R_complex[gate_type], collect(params[string(gate_type,"_discretization")]), params["num_qubits"])

        end
    end
    
    return R_complex    
end

function get_discretized_R_gates(R_type::String, R::Dict{String, Any}, discretization::Array{Float64,1}, num_qubits::Int64)

    if length(discretization) >= 1

        for i=1:length(discretization)
            if     R_type == "RX"
                R_discrete = QCO.RXGate(discretization[i])
            elseif R_type == "RY"
                R_discrete = QCO.RYGate(discretization[i])
            elseif R_type == "RZ"
                R_discrete = QCO.RZGate(discretization[i])
            end

            R["angle_$i"] = Dict{String, Any}("angle" => discretization[i],
                                             "1qubit_rep" => R_discrete,
                                             "$(num_qubits)qubit_rep" => Dict{String, Any}()
                                            )
            
            for i_qu = 1:num_qubits
                R["angle_$i"]["$(num_qubits)qubit_rep"]["qubit_$i_qu"] = QCO.get_full_sized_gate(R_type, num_qubits, matrix = R_discrete, qubit_location = "q$i_qu")
            end

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
                U_complex["U3"] = QCO.get_discretized_U3_gates("U3", U_complex["U3"], collect(float(params["U_θ_discretization"])), collect(float(params["U_ϕ_discretization"])), collect(float(params["U_λ_discretization"])), params["num_qubits"])
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
                                                          "$(num_qubits)qubit_rep" => Dict{String, Any}()
                                                         )
                
                for i_qu=1:num_qubits
                    U["angle_$(counter)"]["$(num_qubits)qubit_rep"]["qubit_$i_qu"] = QCO.get_full_sized_gate(U_type, num_qubits, matrix = U_discrete, qubit_location = "q$i_qu")
                end

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

    if num_qubits > 3
        Memento.error(_LOGGER, "Gates with greater than 3 qubits are not currently supported")
    end
    
    if input == "H_1"
        return QCO.kron_single_gate(num_qubits, QCO.HGate(), "q1")

    elseif input == "H_2"
        return QCO.kron_single_gate(num_qubits, QCO.HGate(), "q2")

    elseif input == "T_1"
        return QCO.kron_single_gate(num_qubits, QCO.TGate(), "q1")

    elseif input == "T_2"
        return QCO.kron_single_gate(num_qubits, QCO.TGate(), "q2")

    elseif input == "Tdagger_1"
        return QCO.kron_single_gate(num_qubits, QCO.TdaggerGate(), "q1")

    elseif input == "Tdagger_2"
        return QCO.kron_single_gate(num_qubits, QCO.TdaggerGate(), "q2")   

    elseif input == "S_1"
        return QCO.kron_single_gate(num_qubits, QCO.SGate(), "q1") 

    elseif input == "S_2"
        return QCO.kron_single_gate(num_qubits, QCO.SGate(), "q2")  

    elseif input == "Sdagger_1"
        return QCO.kron_single_gate(num_qubits, QCO.SdaggerGate(), "q1") 

    elseif input == "Sdagger_2"
        return QCO.kron_single_gate(num_qubits, QCO.SdaggerGate(), "q2")   

    elseif input == "SX_1"
        return QCO.kron_single_gate(num_qubits, QCO.SXGate(), "q1") 

    elseif input == "SX_2"
        return QCO.kron_single_gate(num_qubits, QCO.SXGate(), "q2")  

    elseif input == "SXdagger_1"
        return QCO.kron_single_gate(num_qubits, QCO.SXdaggerGate(), "q1")

    elseif input == "SXdagger_2"
        return QCO.kron_single_gate(num_qubits, QCO.SXdaggerGate(), "q2")   

    elseif input == "X_1"
        return QCO.kron_single_gate(num_qubits, QCO.XGate(), "q1") 

    elseif input == "X_2"
        return QCO.kron_single_gate(num_qubits, QCO.XGate(), "q2")    

    elseif input == "Y_1"
        return QCO.kron_single_gate(num_qubits, QCO.YGate(), "q1")  

    elseif input == "Y_2"
        return QCO.kron_single_gate(num_qubits, QCO.YGate(), "q2")   

    elseif input == "Z_1"
        return QCO.kron_single_gate(num_qubits, QCO.ZGate(), "q1") 

    elseif input == "Z_2"
        return QCO.kron_single_gate(num_qubits, QCO.ZGate(), "q2")     

    # Gates with continuous angle parameters
    elseif input in ["RX", "RY", "RZ", "U3"] 
        return QCO.kron_single_gate(num_qubits, matrix, qubit_location)
    
    end

    # All 2-qubit full-sized gates
    if num_qubits == 2

        if input == "cnot_12"
            return QCO.CNotGate()

        elseif input == "cnot_21"
            return QCO.CNotRevGate()

        elseif input == "cnot_swap"
            return QCO.CNotGate() * QCO.CNotRevGate()

        elseif input == "H_1⊗H_2"
            return kron(QCO.HGate(), QCO.HGate())   
            
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
        
        elseif input == "HCoin"
            return QCO.HCoinGate()
        
        else
            
            Memento.error(_LOGGER, "Specified input elementary gates or the target gate does not exist in the predefined set of gates.")
        end
    end

    # All 3-qubit full-sized gates
    if num_qubits == 3
        
        if input == "H_3"
            return QCO.kron_single_gate(num_qubits, QCO.HGate(), "q3")
    
        elseif input == "T_3"
            return QCO.kron_single_gate(num_qubits, QCO.TGate(), "q3")
    
        elseif input == "Tdagger_3"
            return QCO.kron_single_gate(num_qubits, QCO.TdaggerGate(), "q3")   
    
        elseif input == "S_3"
            return QCO.kron_single_gate(num_qubits, QCO.SGate(), "q3") 
    
        elseif input == "Sdagger_3"
            return QCO.kron_single_gate(num_qubits, QCO.SdaggerGate(), "q3") 
    
        elseif input == "SX_3"
            return QCO.kron_single_gate(num_qubits, QCO.SXGate(), "q3") 
    
        elseif input == "SXdagger_3"
            return QCO.kron_single_gate(num_qubits, QCO.SXdaggerGate(), "q3")
    
        elseif input == "X_3"
            return QCO.kron_single_gate(num_qubits, QCO.XGate(), "q3") 
    
        elseif input == "Y_3"
            return QCO.kron_single_gate(num_qubits, QCO.YGate(), "q3")  
     
        elseif input == "Z_3"
            return QCO.kron_single_gate(num_qubits, QCO.ZGate(), "q3") 
     
        # Gates with continuous angle parameters
        elseif input in ["RX", "RY", "RZ", "U3"] 
            return QCO.kron_single_gate(num_qubits, matrix, qubit_location)
        
        elseif input == "toffoli"
            return QCO.ToffoliGate()

        elseif input == "CSwap"
            return QCO.CSwapGate()

        elseif input == "CCZ"
            return QCO.CCZGate()

        elseif input == "peres"
            return QCO.PeresGate()

        elseif input == "cnot_12"
            return kron(QCO.CNotGate(), QCO.IGate(1))

        elseif input == "cnot_23"
            return kron(QCO.IGate(1), QCO.CNotGate())

        elseif input == "cnot_21"
            return kron(QCO.CNotRevGate(), QCO.IGate(1))

        elseif input == "cnot_32"
            return kron(QCO.IGate(1), QCO.CNotRevGate())

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

        elseif input == "CV_12"
            return kron(QCO.CVGate(), QCO.IGate(1))

        elseif input == "CV_23"
            return kron(QCO.IGate(1), QCO.CVGate())

        elseif input == "CV_21"
            return kron(QCO.CVRevGate(), QCO.IGate(1))

        elseif input == "CV_32"
            return kron(QCO.IGate(1), QCO.CVRevGate())

        elseif input == "CV_13"
            # |0⟩⟨0| ⊗ I ⊗ I 
            control_0 = kron(Array{Complex{Float64},2}([1 0; 0 0]) , kron(QCO.IGate(1), QCO.IGate(1)))
            # |1⟩⟨1| ⊗ I ⊗ V 
            control_1 = kron(Array{Complex{Float64},2}([0 0; 0 1]) , kron(QCO.IGate(1), QCO.SXGate()))
            return control_0 + control_1

        elseif input == "CV_31"
            # I ⊗ I ⊗ |0⟩⟨0| 
            control_0 = kron(QCO.IGate(1), kron(QCO.IGate(1), Array{Complex{Float64},2}([1 0; 0 0])))
            # V ⊗ I ⊗ |1⟩⟨1|
            control_1 = kron(kron(QCO.SXGate(), QCO.IGate(1)), Array{Complex{Float64},2}([0 0; 0 1]))
            return control_0 + control_1

        elseif input == "CVdagger_12"
            return kron(QCO.CVdaggerGate(), QCO.IGate(1))

        elseif input == "CVdagger_23"
            return kron(QCO.IGate(1), QCO.CVdaggerGate())

        elseif input == "CVdagger_21"
            return kron(QCO.CVRevdaggerGate(), QCO.IGate(1))

        elseif input == "CVdagger_32"
            return kron(QCO.IGate(1), QCO.CVRevdaggerGate())

        elseif input == "CVdagger_13"
            # |0⟩⟨0| ⊗ I ⊗ I 
            control_0 = kron(Array{Complex{Float64},2}([1 0; 0 0]) , kron(QCO.IGate(1), QCO.IGate(1)))
            # |1⟩⟨1| ⊗ I ⊗ V 
            control_1 = kron(Array{Complex{Float64},2}([0 0; 0 1]) , kron(QCO.IGate(1), QCO.SXdaggerGate()))
            return control_0 + control_1

        elseif input == "CVdagger_31"
            # I ⊗ I ⊗ |0⟩⟨0| 
            control_0 = kron(QCO.IGate(1), kron(QCO.IGate(1), Array{Complex{Float64},2}([1 0; 0 0])))
            # V ⊗ I ⊗ |1⟩⟨1|
            control_1 = kron(kron(QCO.SXdaggerGate(), QCO.IGate(1)), Array{Complex{Float64},2}([0 0; 0 1]))
            return control_0 + control_1

        else
            Memento.error(_LOGGER, "Specified input elementary gates or the target gate does not exist in the predefined set of gates.")
        end
    
    end

end

function get_input_circuit_dict(input_circuit::Vector{Tuple{Int64,String}}, params::Dict{String,Any})

    input_circuit_dict = Dict{String, Any}()

    status = true
    gate_type = []
    
    for i = 1:length(input_circuit)
        if !(input_circuit[i][2] in params["elementary_gates"])
            status = false
            gate_type = input_circuit[i][2]
        end
    end

    if status  
        for i = 1:length(input_circuit)
    
            if i == input_circuit[i][1]
                input_circuit_dict["$i"] = Dict{String, Any}("depth" => input_circuit[i][1],
                                                             "gate" => input_circuit[i][2])
                # Later: add support for universal and rotation gates here
            else
                input_circuit_dict = Dict{String, Any}()          
                Memento.warn(_LOGGER, "Neglecting the input circuit as multiple gates cannot be input at the same depth")
                break
            end

        end
    else
        Memento.warn(_LOGGER, "Neglecting the input circuit as gate $gate_type is not in input elementary gates")
    end
    
    return input_circuit_dict
end 
