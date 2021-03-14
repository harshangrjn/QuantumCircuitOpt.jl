"""
Given a vector of input with the names of gates (see examples folder), `get_quantum_gates` function 
returns the corresponding elementary gates in the three-dimensional complex matrix form. 
""" 
function get_quantum_gates(params::Dict{String, Any})
    n_qubits = params["n_qubits"]
    n_gates = params["elementary_gates"]

    M_complex = Array{Complex{Float64},3}(zeros(2^n_qubits, 2^n_qubits, length(n_gates)))
    for i=1:length(n_gates)
        M_complex[:,:,i] = get_full_sized_gate(n_gates[i], n_qubits)
    end
    T_complex = get_full_sized_gate(params["target_gate"], n_qubits)

    if ((size(M_complex[:,:,1])[1] != size(T_complex)[1]) || (size(M_complex[:,:,1])[2] != size(T_complex)[2]))
        Memento.error(_LOGGER, "Dimension mis-match for n_gates elementary gates vs. the target gate.")
    end

    return M_complex, T_complex
end

function get_full_sized_gate(input::String, n_qubits::Int64)
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
        elseif input == "H⊗H"
            return kron(gates["hadamard_H"], gates["hadamard_H"])
        elseif input == "Z1"
            return kron(gates["ph_shift_Z"], gates[I_2]) 
        elseif input == "Z2"
            return kron(gates[I_2], gates["ph_shift_Z"])        
        elseif input == "T1"
            return kron(gates["ph_shift_T"], gates[I_2]) 
        elseif input == "T2"
            return kron(gates[I_2], gates["ph_shift_T"])   
        elseif input == "T1_conjugate"
            return kron(gates["ph_shift_T_conj"], gates[I_2]) 
        elseif input == "T2_conjugate"
            return kron(gates[I_2], gates["ph_shift_T_conj"])   
        elseif input == "S1"
            return kron(gates["ph_shift_S"], gates[I_2]) 
        elseif input == "S2"
            return kron(gates[I_2], gates["ph_shift_S"])  
        elseif input == "X1"
            return kron(gates["pauli_X"], gates[I_2]) 
        elseif input == "X2"
            return kron(gates[I_2], gates["pauli_X"])   
        elseif input == "Y1"
            return kron(gates["pauli_Y"], gates[I_2]) 
        elseif input == "Y2"
            return kron(gates[I_2], gates["pauli_Y"])   
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

function get_data(params::Dict{String, Any})
    n_gates = length(params["elementary_gates"])

    M_complex, T_complex = get_quantum_gates(params)
    M_real = zeros(2*size(M_complex)[1], 2*size(M_complex)[2], size(M_complex)[3])
    for d=1:n_gates
        M_real[:,:,d] = get_complex_to_real_matrix(M_complex[:,:,d])
    end

    if params["initial_gate"] == "Identity"
        M_initial = Matrix(LA.I, 2^(params["n_qubits"]+1), 2^(params["n_qubits"]+1))
    end
    # Add code here to support non-identity as an initial condition gate. 
    
    data = Dict{String, Any}("n_qubits" => params["n_qubits"],
                             "depth" => params["D"],
                             "M_real" => M_real,
                             "M_initial" => M_initial,
                             "Target_real" => get_complex_to_real_matrix(T_complex),
                             "elementary_gates" => params["elementary_gates"],
                             "objective" => params["objective"],
                             "optimizer" => params["optimizer"],
                             "presolve" => params["presolve"],
                             "optimizer_log" => params["optimizer_log"],                           
                             "binary_relax" => params["binary_relax"]
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
    elementary_gates = Dict{String, Any}("I_2" => I_2,
                                         "pauli_X" => pauli_X, 
                                         "pauli_Y" => pauli_Y,
                                         "hadamard_H" => hadamard_H,
                                         "ph_shift_Z" => ph_shift_Z,
                                         "ph_shift_T" => ph_shift_T,
                                         "cnot_12" => cnot_12,
                                         "cnot_21" => cnot_21
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

function get_pauli_rotation_gates(θ::Number)
    #input angles in radians
    @assert 0 <= θ <= π
    R_x = Array{Complex{Float64},2}([cos(θ/2) -(sin(θ/2))im; -(sin(θ/2))im cos(θ/2)])
    R_y = Array{Complex{Float64},2}([cos(θ/2) -(sin(θ/2)); (sin(θ/2)) cos(θ/2)])
    R_z = Array{Complex{Float64},2}([(cos(θ/2) - (sin(θ/2))im) 0; 0 (cos(θ/2) + (sin(θ/2))im)])
    pauli_rotation_gates = Dict{String, Any}("R_x" => verify_tolerances_complex_values(R_x), 
                                             "R_y" => verify_tolerances_complex_values(R_y), 
                                             "R_z" => verify_tolerances_complex_values(R_z)
                                            )

    return pauli_rotation_gates
end

function get_universal_gate(θ::Number, ϕ::Number, λ::Number)
    #input angles in radians
    @assert 0 <= θ <= π
    @assert 0 <= ϕ <= 2*π
    @assert 0 <= λ <= 2*π

    U3_11 = cos(θ/2)   
    U3_12 = (cos(λ) - (sin(λ))im)*sin(θ/2)
    U3_21 = (cos(ϕ) + (sin(ϕ))im)*sin(θ/2)
    U3_22 = (cos(λ + ϕ) + (sin(λ + ϕ))im)*cos(θ/2)
    U3 = Array{Complex{Float64},2}([U3_11 U3_12; U3_21 U3_22])

    return verify_tolerances_complex_values(U3)
end
