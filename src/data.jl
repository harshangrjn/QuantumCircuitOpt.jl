"""
Given a vector of input with the names of gates (see examples folder), `get_quantum_gates` function 
returns the corresponding elementary gates in the three-dimensional complex matrix form. 
""" 
function get_quantum_gates(params::Dict{String, Any})
    n_qubits = params["n_qubits"]
    input = params["elementary_gates"]
    M = zeros(Complex{Float64}, 2^n_qubits, 2^n_qubits, length(input))
    for i=1:length(input)
        M[:,:,i] = get_full_sized_gate(input[i])
    end
    T = get_full_sized_gate(params["target_gate"])

    return M, T
end

function get_full_sized_gate(input::String)
    gates = get_elementary_gates()
    if input == "H1"
        return kron(gates["hadamard_H"], gates["I_2"])
    elseif input == "H2"
        return kron(gates["I_2"], gates["hadamard_H"])
    elseif input == "cnot_12"
        return gates["cnot_12"]
    elseif input == "cnot_21"
        return gates["cnot_21"]
    elseif input == "H⊗H"
        return kron(gates["hadamard_H"], gates["hadamard_H"])
    # Add here the other gates from get_elementary_gates
    else
        Memento.error(_LOGGER, "Specified input elementary gates or the target gate does not exist in the predefined set of gates.")
    end
end

function get_data(params::Dict{String, Any})
    M, T = get_quantum_gates(params)
    data = Dict{String, Any}("n_qubits" => params["n_qubits"],
                             "depth" => params["D"],
                             "M_complex" => M,
                             "T_complex" => T,
                             "elementary_gates" => params["elementary_gates"],
                             "optimizer" => params["optimizer"],
                             "presolve" => params["presolve"],
                             "optimizer_log" => params["optimizer_log"],                           
                             "lp_relax" => params["lp_relax"])
    return data
end

function get_elementary_gates()
    # 1-qubit gates 
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
    
    # 2-qubit gates 
    I_3 = Array{Complex{Float64},2}(Matrix(LA.I, 2^3, 2^3)) 

    toffoli = Array{Complex{Float64},2}([1  0  0  0  0  0  0  0
                                         0  1  0  0  0  0  0  0
                                         0  0  1  0  0  0  0  0
                                         0  0  0  1  0  0  0  0
                                         0  0  0  0  1  0  0  0
                                         0  0  0  0  0  1  0  0
                                         0  0  0  0  0  0  0  1
                                         0  0  0  0  0  0  1  0])
    
    cnot_13 = Array{Complex{Float64},2}([1	0	0	0	0	0	0	0
                                         0	1	0	0	0	0	0	0
                                         0	0	1	0	0	0	0	0
                                         0	0	0	1	0	0	0	0
                                         0	0	0	0	0	1	0	0
                                         0	0	0	0	1	0	0	0
                                         0	0	0	0	0	0	0	1
                                         0	0	0	0	0	0	1	0]) 
    
    cnot_31 = Array{Complex{Float64},2}([1	0	0	0	0	0	0	0
                                         0	0	0	0	0	1	0	0
                                         0	0	1	0	0	0	0	0
                                         0	0	0	0	0	0	0	1
                                         0	0	0	0	1	0	0	0
                                         0	1	0	0	0	0	0	0
                                         0	0	0	0	0	0	1	0
                                         0	0	0	1	0	0	0	0]) 

    elementary_gates = Dict{String, Any}("I_2" => I_2,
                                         "pauli_X" => pauli_X, 
                                         "pauli_Y" => pauli_Y,
                                         "hadamard_H" => hadamard_H,
                                         "ph_shift_Z" => ph_shift_Z,
                                         "ph_shift_T" => ph_shift_T,
                                         "cnot_12" => cnot_12,
                                         "cnot_21" => cnot_21,
                                         "I_3" => I_3,
                                         "toffoli" => toffoli,
                                         "cnot_13" => cnot_13,
                                         "cnot_31" => cnot_31
                                        )
    return elementary_gates
end

function get_pauli_rotation_gates(θ::Float64)
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

function get_universal_gate(θ::Float64, ϕ::Float64, λ::Float64)
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
