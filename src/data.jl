"""
Given a vector of input with the names of gates (see examples folder), `get_quantum_gates` function 
returns the corresponding elementary gates in the three-dimensional complex matrix form. 
""" 
function get_quantum_gates(params::Dict{String, Any}, n_gates::Int64, elementary_gates::Array{String,1})

    n_qubits = params["n_qubits"]
    
    R_gates_ids = findall(x -> startswith(x, "R"), elementary_gates)

    R_complex_dict = Dict{}
    if !isempty(R_gates_ids)
        R_complex_dict = get_all_R_gates(params, elementary_gates)
    end
    
    M_complex_dict = Dict{String, Any}()

    counter = 1
    for i=1:length(elementary_gates)
        if startswith(elementary_gates[i], "R")
            for j in keys(R_complex_dict)
                if (j == elementary_gates[i])
                    for k in keys(R_complex_dict[j])
                        for l in keys(R_complex_dict[j][k]["2qubit_rep"])
                            
                            M_complex_dict["$counter"] = Dict{String, Any}("type" => j,
                                                                           "angle" => R_complex_dict[j][k]["angle"],
                                                                           "qubit_location" => l,
                                                                           "matrix" => R_complex_dict[j][k]["2qubit_rep"][l])
                            counter += 1

                        end
                    end
                end
            end
        else     
    
            M_complex_dict["$counter"] = Dict{String, Any}("type" => elementary_gates[i],
                                                           "angle" => "na",
                                                           "qubit_location" => "na",
                                                           "matrix" => get_full_sized_gate(elementary_gates[i], n_qubits))
            counter += 1

        end
    end
    
    M_complex = Array{Complex{Float64},3}(zeros(2^n_qubits, 2^n_qubits, n_gates))
    
    for i=1:n_gates
        M_complex[:,:,i] = M_complex_dict["$i"]["matrix"]
    end

    T_complex = get_full_sized_gate(params["target_gate"], n_qubits)

    if ((size(M_complex[:,:,1])[1] != size(T_complex)[1]) || (size(M_complex[:,:,1])[2] != size(T_complex)[2]))
        Memento.error(_LOGGER, "Dimension mis-match for n_gates elementary gates vs. the target gate.")
    end
 
    return M_complex_dict, M_complex, T_complex
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
                R_complex["R_x"] = get_discretized_R_gates("R_x", R_complex["R_x"], collect(params["R_x_discretization"]), params["n_qubits"])
            end

            if (elementary_gates[R_gates_ids[i]] == "R_y")
                if isempty(params["R_y_discretization"])
                    Memento.error(_LOGGER, "No discretized angle was specified for R_y gate. Input at least one angle")
                end    

                R_complex["R_y"] = Dict{String, Any}()             
                R_complex["R_y"] = get_discretized_R_gates("R_y", R_complex["R_y"], collect(params["R_y_discretization"]), params["n_qubits"])
            end

            if (elementary_gates[R_gates_ids[i]] == "R_z")
                if isempty(params["R_z_discretization"])
                    Memento.error(_LOGGER, "No discretized angle was specified for R_z gate. Input at least one angle")
                end         

                R_complex["R_z"] = Dict{String, Any}()             
                R_complex["R_z"] = get_discretized_R_gates("R_z", R_complex["R_z"], collect(params["R_z_discretization"]), params["n_qubits"])
            end

        end
    end
    
    return R_complex    
end

function get_discretized_R_gates(R_type::String, R::Dict{String, Any}, discretization::Array{Float64,1}, n_qubits::Int64)

    if length(discretization) >= 1

        for i=1:length(discretization)
            R_discrete = get_pauli_rotation_gates(discretization[i])[R_type]
            R["angle_$i"] = Dict{String, Any}("angle" => discretization[i],
                                             "1qubit_rep" => R_discrete,
                                             "2qubit_rep" => Dict{String, Any}("qubit_1" => get_full_sized_gate(R_type, n_qubits, M=R_discrete, qubit_location = "qubit_1"),
                                                                               "qubit_2" => get_full_sized_gate(R_type, n_qubits, M=R_discrete, qubit_location = "qubit_2")
                                                                              ) 
                                            )            
        end
    end 

    return R 
end

function get_discretized_U3_gates(θ_dicretization::Array{Number,1}, ϕ_dicretization::Array{Number,1}, λ_dicretization::Array{Number,1}, n_qubits::Int64)
    n_gates_U = length(θ_dicretization) * length(ϕ_dicretization) * length(λ_dicretization)
    
    U = Array{Complex{Float64},3}(zeros(2^n_qubits, 2^n_qubits, n_gates_U))

    counter = 1
    for i=1:length(θ_dicretization)
        for j=1:length(ϕ_dicretization)
            for k=1:length(λ_dicretization)
                U[:,:,counter] = get_u3_gate(θ_dicretization[i], ϕ_dicretization[j], λ_dicretization[k])
                counter += 1
            end
        end
    end
    
    return U
end

# function get_full_sized_gate(input::String, n_qubits::Int64; M::Array{Complex{Float64},2} = nothing, qubit_location::String = nothing)
function get_full_sized_gate(input::String, n_qubits::Int64; M = nothing, qubit_location = nothing)

    gates = get_elementary_gates(n_qubits)
    # All 2-qubit full-sized gates
    if n_qubits == 2
        
        if input == "Identity"
            return kron(gates["I_2"], gates["I_2"])

        elseif input == "H1"
            return kron(gates["hadamard_H"], gates["I_2"])

        elseif input == "H2"
            return kron(gates["I_2"], gates["hadamard_H"])

        elseif input == "cnot_12"
            return gates["cnot_12"]

        elseif input == "cnot_21"
            return gates["cnot_21"]

        elseif input == "cnot_swap"
            return gates["cnot_12"] * gates["cnot_21"]

        elseif input == "H1⊗H2"
            return kron(gates["hadamard_H"], gates["hadamard_H"])

        elseif input == "Z1"
            return kron(gates["ph_shift_Z"], gates["I_2"]) 

        elseif input == "Z2"
            return kron(gates["I_2"], gates["ph_shift_Z"])        

        elseif input == "T1"
            return kron(gates["ph_shift_T"], gates["I_2"]) 

        elseif input == "T2"
            return kron(gates["I_2"], gates["ph_shift_T"])   

        elseif input == "T1_conjugate"
            return kron(gates["ph_shift_T_conj"], gates["I_2"]) 

        elseif input == "T2_conjugate"
            return kron(gates["I_2"], gates["ph_shift_T_conj"])   

        elseif input == "S1"
            return kron(gates["ph_shift_S"], gates["I_2"]) 

        elseif input == "S2"
            return kron(gates["I_2"], gates["ph_shift_S"])  

        elseif input == "X1"
            return kron(gates["pauli_X"], gates["I_2"]) 

        elseif input == "X2"
            return kron(gates["I_2"], gates["pauli_X"])   

        elseif input == "Y1"
            return kron(gates["pauli_Y"], gates["I_2"]) 

        elseif input == "Y2"
            return kron(gates["I_2"], gates["pauli_Y"])   
            
        elseif input == "controlled_Z"
            return gates["controlled_Z"]

        elseif input == "controlled_H_12"
            return gates["controlled_H_12"]

        elseif input == "controlled_V"
            return gates["controlled_V"]

        elseif input == "swap"
            return gates["swap"]

        elseif input == "magic_M"
            return gates["magic_M"]

        elseif input == "qft2"
            return gates["qft2"]

        elseif input == "test_R_x_1"
            return kron(gates["test_R_x"], gates["I_2"])

        elseif input == "test_R_x_2"
            return kron(gates["I_2"], gates["test_R_x"])

        elseif input == "test_R_y_1"
            return kron(gates["test_R_y"], gates["I_2"])

        elseif input == "test_R_y_2"
            return kron(gates["I_2"], gates["test_R_y"])

        elseif input == "test_R_z_1"
            return kron(gates["test_R_z"], gates["I_2"])

        elseif input == "test_R_z_2"
            return kron(gates["I_2"], gates["test_R_z"])

        elseif input in ["R_x", "R_y", "R_z"] 

            if qubit_location == "qubit_1"
                return kron(M, gates["I_2"])
            elseif qubit_location == "qubit_2"
                return kron(gates["I_2"], M)
            end

        # Add here the other gates from get_elementary_gates
        else
            Memento.error(_LOGGER, "Specified input elementary gates or the target gate does not exist in the predefined set of gates.")
        end
    end
    # All 3-qubit full-sized gates
    if n_qubits == 3
        if input == "toffoli"
            return gates["toffoli"]

        elseif input == "cnot_13"
            return gates["cnot_13"]

        elseif input == "cnot_31"
            return gates["cnot_31"]

        else
            Memento.error(_LOGGER, "Specified input elementary gates or the target gate does not exist in the predefined set of gates.")
        end
    end
end

function get_total_number_of_input_gates(params::Dict{String, Any}, elementary_gates::Array{String,1})

    R_gates = findall(x -> startswith(x, "R"), (elementary_gates))
    U_gates = findall(x -> startswith(x, "U") || startswith(x, "u"), (elementary_gates))

    n_gates = 0
    
    if isempty(R_gates) && isempty(U_gates)
        n_gates = length(elementary_gates)

    elseif !isempty(R_gates)
        # Twice to handle R_x gates on both qubits

        num_R_x = 0; num_R_y = 0; num_R_z = 0;
        
        if ("R_x" in elementary_gates)
            num_R_x = 2*length(params["R_x_discretization"])
            if num_R_x == 0
                Memento.error(_LOGGER, "Dicsretization not specified for R_x gate")
            end
        end
        
        if ("R_y" in elementary_gates)
            num_R_y = 2*length(params["R_y_discretization"])
            if num_R_y == 0 
                Memento.error(_LOGGER, "Dicsretization not specified for R_y gate")
            end 
        end

        if ("R_z" in elementary_gates)
            num_R_z = 2*length(params["R_z_discretization"])
            if num_R_z == 0
                Memento.error(_LOGGER, "Dicsretization not specified for R_z gate")
            end
        end

        n_gates = length(elementary_gates) - length(R_gates) + num_R_x + num_R_y + num_R_z

    elseif !isempty(U_gates)
        # Twice to handle U gates on both qubits
        num_U_θ = 0; num_U_ϕ = 0; num_U_λ = 0;
        
        if ("U_θ" in elementary_gates)
            num_U_θ = 2*length(params["U_θ_discretization"])
            if num_U_θ == 0
                Memento.error(_LOGGER, "Dicsretization not specified for U_θ gate")
            end
        end
        if ("U_ϕ" in elementary_gates)
            num_U_ϕ = 2*length(params["U_ϕ_discretization"])
            if num_U_ϕ == 0
                Memento.error(_LOGGER, "Dicsretization not specified for U_ϕ gate")
            end
        end
        if ("U_λ" in elementary_gates)
            num_U_λ = 2*length(params["U_λ_discretization"])
            if num_U_λ == 0
                Memento.error(_LOGGER, "Dicsretization not specified for U_λ gate")
            end
        end

        n_gates = length(elementary_gates) - length(U_gates) + (num_U_θ * num_U_ϕ * num_U_λ)

    elseif !isempty(R_gates) && !isempty(U_gates)
        num_R_x = 0; num_R_y = 0; num_R_z = 0;
        
        if ("R_x" in elementary_gates)
            num_R_x = 2*length(params["R_x_discretization"])
            if (num_R_x == 0 )
                Memento.error(_LOGGER, "Dicsretization not specified for R_x gate")
            end
        end
        
        if ("R_y" in elementary_gates)
            num_R_y = 2*length(params["R_y_discretization"])
            if (num_R_y == 0 )
                Memento.error(_LOGGER, "Dicsretization not specified for R_y gate")
            end
        end

        if ("R_z" in elementary_gates)
            num_R_z = 2*length(params["R_z_discretization"])
            if (num_R_z == 0 )
                Memento.error(_LOGGER, "Dicsretization not specified for R_z gate")
            end
        end

        num_U_θ = 0; num_U_ϕ = 0; num_U_λ = 0;
        
        if ("U_θ" in elementary_gates)
            num_U_θ = 2*length(params["U_θ_discretization"])
            if num_U_θ == 0
                Memento.error(_LOGGER, "Dicsretization not specified for U_θ gate")
            end
        end
        if ("U_ϕ" in elementary_gates)
            num_U_ϕ = 2*length(params["U_ϕ_discretization"])
            if num_U_ϕ == 0
                Memento.error(_LOGGER, "Dicsretization not specified for U_ϕ gate")
            end
        end
        if ("U_λ" in elementary_gates)
            num_U_λ = 2*length(params["U_λ_discretization"])
            if num_U_λ == 0
                Memento.error(_LOGGER, "Dicsretization not specified for U_λ gate")
            end
        end

        n_gates = length(elementary_gates) - length(R_gates) - length(U_gates) + num_R_x + num_R_y + num_R_z + (num_U_θ * num_U_ϕ * num_U_λ)

    end
    
    if n_gates == 1
        Memento.error(_LOGGER, "Input at least two unique elementary gates for non-trivial solutions")
    else
        return n_gates
    end
end

function get_data(params::Dict{String, Any})
    if isempty(params["elementary_gates"])
        Memento.error(_LOGGER, "Input elementary gates are empty. Enter at least two unique unitary gates")
    end

    elementary_gates = unique(params["elementary_gates"])
    
    if length(elementary_gates) < length(params["elementary_gates"])
        Memento.warn(_LOGGER, "Eliminating non-unique gates in the input elementary gates")
    end
    
    n_gates = get_total_number_of_input_gates(params, elementary_gates)

    M_complex_dict, M_complex, T_complex = get_quantum_gates(params, n_gates, elementary_gates)
   
    M_real = zeros(2*size(M_complex)[1], 2*size(M_complex)[2], size(M_complex)[3])

    for d=1:n_gates
        M_real[:,:,d] = get_complex_to_real_matrix(M_complex[:,:,d])
    end

    if params["initial_gate"] == "Identity"
        M_initial = Matrix(LA.I, 2^(params["n_qubits"]+1), 2^(params["n_qubits"]+1))
    else
        Memento.error(_LOGGER, "Currently non-identity gate is not supported as the initial condition")
    end
    # Add code here to support non-identity as an initial condition gate. 
    
    data = Dict{String, Any}("n_qubits" => params["n_qubits"],
                             "depth" => params["D"],
                             "M_complex_dict" => M_complex_dict,
                             "M_real" => M_real,
                             "M_initial" => M_initial,
                             "Target_real" => get_complex_to_real_matrix(T_complex),
                             "elementary_gates" => elementary_gates,
                             "target_gate" => params["target_gate"],
                             "objective" => params["objective"],
                             "decomposition_type" => params["decomposition_type"],
                             "optimizer" => params["optimizer"],                         
                             "relax_integrality" => params["relax_integrality"]
                             )
    return data
end

"""
Given the number of qubits (`n_qubits`), `get_elementary_gates` function 
returns all the elementary gates in the basic `2⨉2` form, without applying 
kronecker tensor product operations. 
""" 
function get_elementary_gates(n_qubits::Int64)
    # 1 and 2-qubit gates 
    I_2 = Array{Complex{Float64},2}([1 0; 0 1])
    pauli_X = Array{Complex{Float64},2}([0 1; 1 0])
    pauli_Y = Array{Complex{Float64},2}([0 -im; im 0])  
    
    hadamard_H = Array{Complex{Float64},2}(1/sqrt(2)*[1 1; 1 -1])
    ph_shift_Z = Array{Complex{Float64},2}([1 0; 0 -1]) # Also called pauli-Z gate (π/8 gate)
    ph_shift_S = Array{Complex{Float64},2}([1 0; 0 im])
    ph_shift_T = Array{Complex{Float64},2}([1 0; 0 (1/sqrt(2)) + (1/sqrt(2))im])
    ph_shift_T_conj = Array{Complex{Float64},2}([1 0; 0 (1/sqrt(2)) - (1/sqrt(2))im])
 
    cnot_12 = Array{Complex{Float64},2}([1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]) 
    cnot_21 = Array{Complex{Float64},2}([1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0]) 
    controlled_Z = Array{Complex{Float64},2}([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1])
    # controlled_H_12 = Array{Complex{Float64},2}([1 0 0 0; 0 1/sqrt(2) 0 1/sqrt(2); 0 0 1 0; 0 1/sqrt(2) 0 -1/sqrt(2)])
    controlled_H_12 = Array{Complex{Float64},2}([1 0 0 0; 0 1 0 0; 0 0 1/sqrt(2) 1/sqrt(2); 0  0 1/sqrt(2) -1/sqrt(2)])
    controlled_V = Array{Complex{Float64},2}([1 0 0 0; 0 1 0 0; 0 0 0.5+(0.5)im 0.5-(0.5)im; 0 0 0.5-(0.5)im 0.5+(0.5)im])  #Also called sqrt(CNOT) gate
    swap = Array{Complex{Float64},2}([1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1])
    magic_M = Array{Complex{Float64},2}(1/sqrt(2)*[1 im 0 0; 0 0 im 1; 0 0 im -1; 1 -im 0 0])
    qft2 = Array{Complex{Float64},2}(0.5*[1 1 1 1; 1 im -1 -im; 1 -1 1 -1; 1 -im -1 im])

    test_R_x = Array{Complex{Float64},2}([ 0.92388+0.0im   0.0-0.382683im
                                            0.0-0.382683im  0.92388+0.0im])

    test_R_y = Array{Complex{Float64}, 2}([  0.92388+0.0im  -0.382683+0.0im
                                            0.382683+0.0im    0.92388+0.0im
                                        ])
    test_R_z = Array{Complex{Float64}, 2}([  0.92388-0.382683im      0.0+0.0im
                                                0.0+0.0im       0.92388+0.382683im])

    elementary_gates = Dict{String, Any}("I_2" => I_2,
                                         "pauli_X" => pauli_X, 
                                         "pauli_Y" => pauli_Y,
                                         "hadamard_H" => hadamard_H,
                                         "ph_shift_S" => ph_shift_S,
                                         "ph_shift_Z" => ph_shift_Z,
                                         "ph_shift_T" => ph_shift_T,
                                         "ph_shift_T_conj" => ph_shift_T_conj,
                                         "cnot_12" => cnot_12,
                                         "cnot_21" => cnot_21,
                                         "controlled_Z" => controlled_Z,
                                         "controlled_H_12" => controlled_H_12,
                                         "controlled_V" => controlled_V,
                                         "swap" => swap,
                                         "magic_M" => magic_M,
                                         "qft2" => qft2,
                                         "test_R_x" => test_R_x,
                                         "test_R_y" => test_R_y,
                                         "test_R_z" => test_R_z
                                         )
    
    # 3-qubit gates 
    if n_qubits == 3
        I_3 = sparse(Array{Complex{Float64},2}(Matrix(LA.I, 2^3, 2^3)))

        toffoli = sparse(Array{Complex{Float64},2}([1  0  0  0  0  0  0  0
                                                    0  1  0  0  0  0  0  0
                                                    0  0  1  0  0  0  0  0
                                                    0  0  0  1  0  0  0  0
                                                    0  0  0  0  1  0  0  0
                                                    0  0  0  0  0  1  0  0
                                                    0  0  0  0  0  0  0  1
                                                    0  0  0  0  0  0  1  0]))
        
        cnot_13 = sparse(Array{Complex{Float64},2}([1	0	0	0	0	0	0	0
                                                    0	1	0	0	0	0	0	0
                                                    0	0	1	0	0	0	0	0
                                                    0	0	0	1	0	0	0	0
                                                    0	0	0	0	0	1	0	0
                                                    0	0	0	0	1	0	0	0
                                                    0	0	0	0	0	0	0	1
                                                    0	0	0	0	0	0	1	0]))
        
        cnot_31 = sparse(Array{Complex{Float64},2}([1	0	0	0	0	0	0	0
                                                    0	0	0	0	0	1	0	0
                                                    0	0	1	0	0	0	0	0
                                                    0	0	0	0	0	0	0	1
                                                    0	0	0	0	1	0	0	0
                                                    0	1	0	0	0	0	0	0
                                                    0	0	0	0	0	0	1	0
                                                    0	0	0	1	0	0	0	0]))

        elementary_gates["I_3"] = I_3
        elementary_gates["toffoli"] = toffoli
        elementary_gates["cnot_13"] = cnot_13
        elementary_gates["cnot_31"] = cnot_31
    end
    return elementary_gates
end

"""
    get_pauli_rotation_gates

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
    pauli_rotation_gates = Dict{String, Any}("R_x" => verify_tolerances_complex_values(R_x), 
                                             "R_y" => verify_tolerances_complex_values(R_y), 
                                             "R_z" => verify_tolerances_complex_values(R_z)
                                            )

    return pauli_rotation_gates
end

function get_u1_gate(λ::Number)
    #input angles in radians
    θ = 0
    ϕ = 0

    if !(-2*π <= λ <= 2*π)
        Memento.error(_LOGGER, "λ angle in universal u1 gate is not within valid bounds")
    end

    u1 = get_u3_gate(θ, ϕ, λ)
    return u1
end

function get_u2_gate(ϕ::Number, λ::Number)
    #input angles in radians
    θ = π/2

    if !(-2*π <= λ <= 2*π)
        Memento.error(_LOGGER, "λ angle in universal u2 gate is not within valid bounds")
    end
    if !(-2*π <= ϕ <= 2*π)
        Memento.error(_LOGGER, "ϕ angle in universal u2 gate is not within valid bounds")
    end

    u2 = get_u3_gate(θ, ϕ, λ)
    return u2
end

"""
    get_u3_gate

Given three angles (θ,ϕ,λ),  this function returns the most general form of a single qubit unitar gate.
"""
function get_u3_gate(θ::Number, ϕ::Number, λ::Number)
    #input angles in radians

    if !(-π <= θ <= π)
        Memento.error(_LOGGER, "θ angle in universal u3 gate is not within valid bounds")
    end
    if !(-2*π <= ϕ <= 2*π)
        Memento.error(_LOGGER, "ϕ angle in universal u3 gate is not within valid bounds")
    end
    if !(-2*π <= λ <= 2*π)
        Memento.error(_LOGGER, "λ angle in universal u3 gate is not within valid bounds")
    end

    u3 = Array{Complex{Float64},2}([           cos(θ/2)               -(cos(λ) + (sin(λ))im)*sin(θ/2) 
                                    (cos(ϕ) + (sin(ϕ))im)*sin(θ/2)  (cos(λ+ϕ) + (sin(λ+ϕ))im)*cos(θ/2)])

    return verify_tolerances_complex_values(u3)
end