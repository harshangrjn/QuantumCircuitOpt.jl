"""
    auxiliary_variable_bounds(v::Array{JuMP.VariableRef,1})

Given a vector of JuMP variables (maximum 4 variables), this function returns the worst-case 
bounds, the product of these input variables can admit.  
"""
function auxiliary_variable_bounds(v::Array{JuMP.VariableRef,1}) 

    v_l = [JuMP.lower_bound(v[1]), JuMP.upper_bound(v[1])]
    v_u = [JuMP.lower_bound(v[2]), JuMP.upper_bound(v[2])]
    M = v_l * v_u'
 
    if length(v) == 2 
       # Bilinear
        return minimum(M), maximum(M)
 
    elseif (length(v) == 3) || (length(v) == 4) 
       M1 = zeros(Float64, (2, 2, 2))
       M1[:,:,1] = M * JuMP.lower_bound(v[3])
       M1[:,:,2] = M * JuMP.upper_bound(v[3])
       if length(v) == 3
          # Trilinear
          return minimum(M1), maximum(M1)
       elseif length(v) == 4
          # Quadrilinear
          M2 = zeros(Float64, (2, 2, 4))
          M2[:,:,1:2] = M1 * JuMP.lower_bound(v[4])
          M2[:,:,3:4] = M1 * JuMP.upper_bound(v[4])
          return minimum(M2), maximum(M2)
       end
    end
end
 
"""
    gate_element_bounds(M::Array{Float64,3})

Given a set of elementary gates, `{G1, G2, ... ,Gn}`, this function evaluates 
the range of every co-ordinate of the superimposed gates, over all possible gates.  
"""
function gate_element_bounds(M::Array{Float64,3}) 

    M_l = zeros(size(M)[1], size(M)[2])
    M_u = zeros(size(M)[1], size(M)[2])

    for i = 1:size(M)[1]
        for j = 1:size(M)[2]

            M_l[i,j] = minimum(M[i,j,:])
            M_u[i,j] = maximum(M[i,j,:])

            if M_l[i,j] > M_u[i,j]
                Memento.error(_LOGGER, "Lower and upper bound conflict in the elements of input elementary gates")
            elseif (M_l[i,j] in [-Inf, Inf]) || (M_u[i,j] in [-Inf, Inf])
                Memento.error(_LOGGER, "Unbounded entry detected in the input elementary gates")
            end
        end
    end

    return M_l, M_u
end

"""
    get_commutative_gate_pairs(M::Dict{String,Any}; identity_pairs = true)

Given a dictionary of elementary quantum gates, this function returns all pairs of commuting 
gates. Optional argument, `identity_pairs` can be set to `false` if identity matrix need not be part of the commuting pairs. 
"""
function get_commutative_gate_pairs(M::Dict{String,Any}; identity_in_pairs = true)
    
    num_gates = length(keys(M))
    commute_pairs = Array{Tuple{Int64,Int64},1}()
    commute_pairs_prodIdentity = Array{Tuple{Int64,Int64},1}()

    for i = 1:(num_gates-1)
        for j = (i+1):num_gates
            M_i = M["$i"]["matrix"]
            M_j = M["$j"]["matrix"]

            if ("Identity" in M["$i"]["type"]) || ("Identity" in M["$j"]["type"])
                continue
            end
            
            M_ij = M_i*M_j
            M_ji = M_j*M_i
            Id = Matrix(LA.I, size(M_ij)[1], size(M_ij)[2])

            # Commuting pairs == Identity 
            if isapprox(M_ij, Id, atol=1E-6)
                push!(commute_pairs_prodIdentity, (i,j))
            
            # Commuting pairs != Identity 
            elseif isapprox(M_ij, M_ji, atol = 1E-4)
                push!(commute_pairs, (i, j))

            end

        end
    end

    if identity_in_pairs
        # Commuting pairs involving Identity
        identity_idx = []
        for i in keys(M)
            if "Identity" in M[i]["type"]
                push!(identity_idx, parse(Int64, i))
            end
        end

        if length(identity_idx) > 0
            for i = 1:length(identity_idx)
                for j = 1:num_gates 
                    if j != identity_idx[i]
                        push!(commute_pairs, (j, identity_idx[i]))
                    end
                end
            end
        end

    end 

    return commute_pairs, commute_pairs_prodIdentity
end

"""
    get_redundant_gate_product_pairs(M::Dict{String,Any})

Given a dictionary of elementary quantum gates, this function returns all pairs of gates whose product is 
one of the input elementary gates. For example, let `G_basis = {G1, G2, G3}` be the elementary gates. If `G1*G2 ∈ G_basis`, 
then `(1,2)` is considered as a redundant pair. 
"""
function get_redundant_gate_product_pairs(M::Dict{String,Any})
    num_gates = length(keys(M))
    redundant_pairs_idx = Array{Tuple{Int64,Int64},1}()

    # Non-Identity redundant pairs
    for i = 1:(num_gates-1)
        for j = (i+1):num_gates
            M_i = M["$i"]["matrix"]
            M_j = M["$j"]["matrix"]

            if ("Identity" in M["$i"]["type"]) || ("Identity" in M["$j"]["type"])
                continue
            end
            
            for k = 1:num_gates 
                if (k != i) && (k != j) && !("Identity" in M["$k"]["type"])                
                    M_k = M["$k"]["matrix"]
                
                    if isapprox(M_i*M_j, M_k, atol = 1E-4)
                        push!(redundant_pairs_idx, (i, j))
                    end

                    if isapprox(M_j*M_i, M_k, atol = 1E-4)
                        push!(redundant_pairs_idx, (j, i))
                    end
                end
            end

        end
    end

    return redundant_pairs_idx
end

"""
    get_idempotent_gates(M::Dict{String,Any})

Given the dictionary of complex quantum gates, this function returns the indices of matrices which are self-idempotent 
or idempotent with other set of input gates, excluding the Identity gate.
"""
function get_idempotent_gates(M::Dict{String,Any})
    num_gates = length(keys(M))
    idempotent_gates_idx = Vector{Int64}()

    # Excluding Identity gate in input 
    for i = 1:num_gates
        M_i = M["$i"]["matrix"]
        for j = 1:num_gates
            M_j = M["$j"]["matrix"]
            if ("Identity" in M["$i"]["type"]) || ("Identity" in M["$j"]["type"])
                continue
            end

            if isapprox(M_i^2, M_j, atol=1E-4)
                push!(idempotent_gates_idx, i)
            end

        end
    end

    return idempotent_gates_idx
end

"""
    get_involutory_gates(M::Dict{String,Any})

Given the dictionary of complex gates `G_1, G_2, ..., G_n`, this function returns the indices of these gates 
which are involutory, i.e, `G_i^2 = Identity`, excluding the Identity gate. 
"""
function get_involutory_gates(M::Dict{String,Any})
    num_gates = length(keys(M))
    involutory_gates_idx = Vector{Int64}()

    if num_gates > 0
        n_r = size(M["1"]["matrix"])[1]
        n_c = size(M["1"]["matrix"])[2]
    end
    Id = Matrix{ComplexF64}(LA.I, n_r, n_c)

    # Excluding Identity gate in input 
    for i=1:num_gates
        if !("Identity" in M["$i"]["type"]) && (isapprox((M["$i"]["matrix"])^2, Id, atol=1E-5))
            push!(involutory_gates_idx, i)
        end
    end

    return involutory_gates_idx
end

"""
    complex_to_real_gate(M::Array{Complex{Float64},2})

Given a complex-valued two-dimensional quantum gate of size NxN, this function returns a real-valued gate 
of dimensions 2Nx2N. 
"""
function complex_to_real_gate(M::Array{Complex{Float64},2})

    n = size(M)[1]
    M_real = zeros(2*n, 2*n)
  
    ii = 1; jj = 1;
    for i = collect(1:2:2*n)
        for j = collect(1:2:2*n)

            if isapprox(real(M[ii,jj]), 0, atol=1E-6)
                M_real[i,j] = 0
                M_real[i+1,j+1] = 0
            else
                M_real[i,j] = real(M[ii,jj])
                M_real[i+1,j+1] = real(M[ii,jj])
            end

            if isapprox(imag(M[ii,jj]), 0, atol=1E-6)
                M_real[i,j+1] = 0
                M_real[i+1,j] = 0
            else
                M_real[i,j+1] = imag(M[ii,jj])
                M_real[i+1,j] = -imag(M[ii,jj])
            end
            
            jj += 1
        end
        jj = 1
        ii += 1
    end
  
    return M_real
end

"""
    real_to_complex_gate(M::Array{Complex{Float64},2})

Given a real-valued two-dimensional quantum gate of size 2Nx2N, this function returns a complex-valued gate 
of size NxN, if the input gate is in a valid complex form. 
"""
function real_to_complex_gate(M::Array{Float64,2})
    
    n = size(M)[1]

    if !iseven(n)
        Memento.error(_LOGGER, "Input gate can admit only even numbered columns and rows")
    end
    
    M_complex = zeros(Complex{Float64}, (Int(n/2), Int(n/2)))
  
    ii = 1; jj = 1;
    for i = collect(1:2:n)
        for j = collect(1:2:n)

            if !isapprox(M[i,j], M[i+1, j+1], atol = 1E-5) || !isapprox(M[i+1,j], -M[i,j+1], atol = 1E-5)
                Memento.error(_LOGGER, "Input real form of the complex gate is invalid")
            end

            M_re = M[i,j]
            M_im = M[i,j+1]
            (isapprox(M_re, 0, atol=1E-6)) && (M_re = 0)
            (isapprox(M_im, 0, atol=1E-6)) && (M_im = 0)
            
            M_complex[ii,jj] = complex(M_re, M_im)
            jj += 1
        end
        jj = 1
        ii += 1
    end

    return M_complex
end

"""
    round_complex_values(M::Array{Complex{Float64},2})

Given a complex-valued gate, this function returns a complex-valued gate which 
rounds the values closest to 0 and 1. This is useful to avoid numerical issues. 
"""
function round_complex_values(M::Array{Complex{Float64},2})   
    if size(M)[1] == 0
        Memento.error(_LOGGER, "Input cannot be a scalar")
    end

    if length(size(M)) == 2
        n_r = size(M)[1]
        n_c = size(M)[2]

        M_round = Array{Complex{Float64},2}(zeros(n_r,n_c))
        
        for i=1:n_r
            for j=1:n_c 
                M_round[i,j] = complex(QCO.round_real_value(real(M[i,j])), QCO.round_real_value(imag(M[i,j])))
            end
        end
        return M_round
    end
    
end

"""
    round_real_value(x::Float64)

Given a real-valued number, this function returns a real-value which rounds the values closest to 0 and 1. 
"""
function round_real_value(x::Float64)
    if isapprox(abs(x), 0, atol=1E-6)
        x = 0
    elseif isapprox(x, 1, atol=1E-6)
        x = 1
    end  

    return x
end

"""
    unique_idx(x::AbstractArray{T})

This function returns the indices of unique elements in a given array of scalar or vector inputs. Overall, 
this function computes faster than Julia's built-in `findfirst` command. 
"""
function unique_idx(x::AbstractArray{T}) where T
    uniqueset = Set{T}()
    ex = eachindex(x)
    idxs = Vector{eltype(ex)}()
    for i in ex
        xi = x[i]
        if !(xi in uniqueset)
            push!(idxs, i)
            push!(uniqueset, xi)
        end
    end
    idxs
end

"""
    unique_matrices(M::Array{Float64, 3})

This function returns the unique set of matrices and the corresponding indices 
of unique matrices from the given set of matrices.  
"""
function unique_matrices(M::Array{Float64, 3})
    M[isapprox.(M, 0, atol=1E-6)] .= 0

    M_reshape = [];
    for i=1:size(M)[3]
        push!(M_reshape, round.(reshape(M[:,:,i], size(M)[1]*size(M)[2]), digits=5))
    end

    idx = QCO.unique_idx(M_reshape)

    return M[:,:,idx], idx
end

"""
    kron_single_qubit_gate(num_qubits::Int64, M::Array{Complex{Float64},2}, qubit_loc::String)

Given number of qubits of the circuit, the complex-valued one-qubit gate and the qubit location ("q1","q2',"q3",...),
this function returns a full-sized gate after applying appropriate kronecker products. This function supports any number 
integer-valued qubits.  
"""
function kron_single_qubit_gate(num_qubits::Int64, M::Array{Complex{Float64},2}, qubit_loc::String)
    
    if size(M)[1] != 2
        Memento.error(_LOGGER, "Input should be an one-qubit gate")
    end

    qubit = parse(Int, qubit_loc[2:end])

    if !(qubit in 1:num_qubits)
        Memento.error(_LOGGER, "Input qubit location, $qubit, has to be ∈ [q1,...,q$num_qubits]")
    end

    I = QCO.IGate(1)
    M_kron = 1

    for i = 1:num_qubits
        M_iter = I
        if i == qubit 
            M_iter = M
        end
        M_kron = kron(M_kron, M_iter)
    end

    if size(M_kron)[1] == 2^(num_qubits)
        return M_kron 
    else 
        Memento.error(_LOGGER, "Dimensions mismatch in evaluation of Kronecker of single-qubit gate")
    end

end


"""
    kron_two_qubit_gate(num_qubits::Int64, M::Array{Complex{Float64},2}, c_qubit_loc::String, t_qubit_loc::String)

Given number of qubits of the circuit, the complex-valued two-qubit gate and the control and target qubit locations ("q1","q2',"q3",...),
this function returns a full-sized gate after applying appropriate kronecker products. This function supports any number integer-valued qubits.  
"""
function kron_two_qubit_gate(num_qubits::Int64, M::Array{Complex{Float64},2}, c_qubit_loc::String, t_qubit_loc::String)
    
    if size(M)[1] != 4
        Memento.error(_LOGGER, "Input should be a two-qubit gate")
    end

    c_qubit = parse(Int, c_qubit_loc[2:end])
    t_qubit = parse(Int, t_qubit_loc[2:end])

    if !(c_qubit in 1:num_qubits) || !(t_qubit in 1:num_qubits)
        Memento.error(_LOGGER, "Input control and target qubit locations have to be ∈ [q1,...,q$num_qubits]")
    elseif isapprox(c_qubit, t_qubit, atol = 1E-6)
        Memento.error(_LOGGER, "Control and target qubits cannot be identical for a multi-qubit elementary gate")
    end

    I = QCO.IGate(1)
    Swap = QCO.SwapGate()

    # If on adjacent qubits
    if abs(c_qubit - t_qubit) == 1

        M_kron = 1
        for i = 1:(num_qubits-1)
            M_iter = I
            if (i in [c_qubit, t_qubit]) && (i+1 in [c_qubit, t_qubit])
                M_iter = M
            end
            M_kron = kron(M_kron, M_iter)
        end

    # If on non-adjacent qubits
    # Assuming the qubit numbering starts at 1, and not 0
    elseif abs(c_qubit - t_qubit) >= 2
        ct = abs(c_qubit - t_qubit)
        
        M_sub_depth = 2*ct - 1 
        M_sub_loc = ct
        M_sub_swap_id = ct
        M_sub = Matrix{Complex{Float64}}(LA.I, 2^(ct+1),2^(ct+1))
        
        # Idea is to represent gate_14 = swap_34 * swap_23 * gate_12 * swap_23 * swap_34
        for d=1:M_sub_depth
            M_sub_kron = 1 
            swap_qubits = true

            for i=1:ct
                M_sub_iter = I 

                if (d == M_sub_loc) && (i == 1)
                    M_sub_iter = M 
                elseif !(d == M_sub_loc) && (i == M_sub_swap_id) && swap_qubits
                    M_sub_iter = Swap
                    if (d < M_sub_loc) && (M_sub_swap_id > 2)
                        M_sub_swap_id -= 1
                    elseif (d > M_sub_loc) 
                        M_sub_swap_id += 1
                    end
                    swap_qubits = false
                end

                M_sub_kron = kron(M_sub_kron, M_sub_iter)
            end

            M_sub *= M_sub_kron
        end

        M_kron = 1
        i = 1 
        while i <= num_qubits
            M_iter = I 
            if (i in [c_qubit, t_qubit])
                M_iter = M_sub
                i += (ct+1)
            else 
                i += 1
            end
            M_kron = kron(M_kron, M_iter)
        end 

    end

    if size(M_kron)[1] == 2^(num_qubits)
        return M_kron 
    else 
        Memento.error(_LOGGER, "Dimensions mismatch in evaluation of Kronecker of two-qubit gate")
    end

end

"""
    _parse_gates_with_kron_symbol(s::String)

Given a string with gates separated by kronecker symbols (x), this function parses and returns the vector of gates. For 
example, if the input string is `H_1xCNot_2_3xT_4`, the output will be `Vector{String}(["H_1", "CNot_2_3", "T_4"])`.
"""
function _parse_gates_with_kron_symbol(s::String)

    gates = Vector{String}()
    gate_id = string()
 
    for i = 1:length(s)
       if s[i] != kron_symbol
          gate_id = gate_id * s[i]
       else
          push!(gates, gate_id)
          (i != length(s)) && (gate_id = string())
       end
 
       if i == length(s) 
          push!(gates, gate_id)
       end
    end
 
    return gates
 end

"""
    _parse_gate_string(s::String)

Given a string representing a single gate with qubit numbers separated by symbol `_`, this function parses and returns the vector of qubits on
which the input gate is located. For example, if the input string is `CRX_2_3`, the output will be `Vector{Int64}([2,3])`.
"""
function _parse_gate_string(s::String; type=false, qubits=false)

    if occursin(kron_symbol, s)
        Memento.error(_LOGGER, "Kron symbol is not supported for parsing qubit numbers in $s")
    end
    
    gates = Vector{String}()
    gate_id = string()
 
    for i = 1:length(s)
       if s[i] != qubit_separator
          gate_id = gate_id * s[i]
       else
          push!(gates, gate_id)
          (i != length(s)) && (gate_id = string())
       end
 
       if i == length(s) 
          push!(gates, gate_id)
       end
    end
    
    if type && qubits
        return gates[1], parse.(Int, gates[2:end]) # Assuming 1st element is the gate type/name
    elseif type 
        return gates[1]
    elseif qubits 
        return parse.(Int, gates[2:end])
    end

 end
 
"""
    is_gate_real(M::Array{Complex{Float64},2})

Given a complex-valued quantum gate, M, this function returns if M has purely real parts or not as it's elements. 
"""
function is_gate_real(M::Array{Complex{Float64},2})
    M_imag = imag(M)
    n_r = size(M_imag)[1]
    n_c = size(M_imag)[2]

    if sum(isapprox.(M_imag, zeros(n_r, n_c), atol=1E-6)) == n_r*n_c
        return true
    else
        return false
    end
 end

"""
    _get_constraint_slope_intercept(vertex1::Vector{Float64}, vertex2::Vector{Float64})

Given co-ordinates of two points in a plane, this function returns the slope (m) and intercept (c) of the 
line joining these two points. 
"""
function _get_constraint_slope_intercept(vertex1::Vector{Float64}, vertex2::Vector{Float64})
    
    if isapprox.(vertex1, vertex2, atol=1E-6) == [true, true]
        Memento.warn(_LOGGER, "Invalid slope and intercept for two identical vertices")
        return
    end

    if isapprox(vertex1[1], vertex2[1], atol = 1E-6)
        return Inf, Inf 
    else
        m = QCO.round_real_value((vertex2[2] - vertex1[2]) / (vertex2[1] - vertex1[1]))
        c = QCO.round_real_value(vertex1[2] - (m * vertex1[1]))

        return m,c
    end

 end

"""
    is_multi_qubit_gate(gate::String)

Given the input gate string, this function returns a boolean if the input gate is a multi qubit gate or not. 
For example, for a 2-qubit gate `CRZ_1_2`, output is `true`. 
"""
 function is_multi_qubit_gate(gate::String)
    
    if occursin(kron_symbol, gate)
        return true
    end

    qubit_loc = QCO._parse_gate_string(gate, qubits = true)
    
    if length(qubit_loc) > 1 
        return true
    elseif length(qubit_loc) == 1 
        return false 
    else 
        Memento.error(_LOGGER, "Atleast one qubit has to be specified for an input gate")
    end

end