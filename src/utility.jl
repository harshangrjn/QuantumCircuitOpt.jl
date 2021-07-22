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

Given a set of elementary gates, {G_1, G_2, ... G_n}, `gate_element_bounds` function evaluates 
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
                Memento.error(_LOGGER, "lower bound cannot be greater than upper bound in the elements of input elementary gates")
            end
            if (M_l[i,j] in [-Inf, Inf]) || (M_u[i,j] in [-Inf, Inf])
                Memento.error(_LOGGER, "one of the input elementary gates has at least one unbounded entry, which is invalid")
            end
        end
    end
    return M_l, M_u
end

"""
    get_commutative_gate_pairs(M::Dict{String,Any}; identity_pairs = true)

Given a dictionary of gates, as processed in `src/data.jl`, `get_commutative_gate_pairs` outputs all pairs of commuting 
matrices. Optional argument, `identity_pairs` can be set to `false` if identity matrix need not be part of the commuting pairs. 
"""
function get_commutative_gate_pairs(M::Dict{String,Any}; identity_in_pairs = true)
    
    depth = length(keys(M))
    commute_pairs = Array{Tuple{Int64,Int64},1}()
    commute_pairs_prodIdentity = Array{Tuple{Int64,Int64},1}()

    for i = 1:(depth-1)
        for j = (i+1):depth
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
                for j = 1:depth
                    if j != identity_idx[i]
                        push!(commute_pairs, (j, identity_idx[i]))
                    end
                end
            end
        end

    end 

    return commute_pairs, commute_pairs_prodIdentity

end

function get_redundant_gate_product_pairs(M::Dict{String,Any})
    depth = length(keys(M))
    redundant_pairs = Array{Tuple{Int64,Int64},1}()

    # Non-Identity redundant pairs
    for i = 1:(depth-1)
        for j = (i+1):depth
            M_i = M["$i"]["matrix"]
            M_j = M["$j"]["matrix"]

            if ("Identity" in M["$i"]["type"]) || ("Identity" in M["$j"]["type"])
                continue
            end
            
            for k = 1:depth
                if (k != i) && (k != j) && !("Identity" in M["$k"]["type"])                
                    M_k = M["$k"]["matrix"]
                
                    if isapprox(M_i*M_j, M_k, atol = 1E-4)
                        push!(redundant_pairs, (i, j))
                    end

                    if isapprox(M_j*M_i, M_k, atol = 1E-4)
                        push!(redundant_pairs, (j, i))
                    end
                end
            end

        end
    end

    return redundant_pairs
end

"""
    get_idempotent_gates(M::Dict{String,Any})

Given the dictionary of complex gates, this function returns the indices of matrices which are self-idempotent 
or idempotent with other set of input gates, excluding the Identity gate.
"""
function get_idempotent_gates(M::Dict{String,Any})
    depth = length(keys(M))
    idempotent_gates = Vector{Int64}()

    # Excluding Identity gates in input 
    for i=1:depth
        M_i = M["$i"]["matrix"]
        for j=1:depth
            M_j = M["$j"]["matrix"]
            if ("Identity" in M["$i"]["type"]) || ("Identity" in M["$j"]["type"])
                continue
            end

            if isapprox(M_i^2, M_j, atol=1E-4)
                push!(idempotent_gates, i)
            end

        end
    end

    return idempotent_gates
end

"""
    complex_to_real_matrix(M::Array{Complex{Float64},2})

Given a complex-valued 2D matrix of size NxN, complex_to_real_matrix function returns a real-valued matrix 
of size 2Nx2N. 
"""
function complex_to_real_matrix(M::Array{Complex{Float64},2})

    n = size(M)[1]
    M_real = zeros(2*n, 2*n)
  
    ii = 1; jj = 1;
    for i = collect(1:2:2*n)
        for j = collect(1:2:2*n)
            M_real[i,j] = real(M[ii,jj])
            M_real[i+1,j+1] = real(M[ii,jj])
            if imag(M[ii,jj]) == 0
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
    real_to_complex_matrix(M::Array{Complex{Float64},2})

Given a real-valued 2D matrix of size 2Nx2N, real_to_complex_matrix function returns a complex-valued matrix 
of size NxN, if the input matrix is in a valid complex matrix form. 
"""
function real_to_complex_matrix(M::Array{Float64,2})
    
    n = size(M)[1]

    if !iseven(n)
        Memento.error(_LOGGER, "Input real matrix can admit only even numbered columns and rows")
    end
    
    M_complex = zeros(Complex{Float64}, (Int(n/2), Int(n/2)))
  
    ii = 1; jj = 1;
    for i = collect(1:2:n)
        for j = collect(1:2:n)

            if !isapprox(M[i,j], M[i+1, j+1], atol = 1E-5) || !isapprox(M[i+1,j], -M[i,j+1], atol = 1E-5)
                Memento.error(_LOGGER, "Input real matrix cannot be converted into a valid complex matrix form")
            end

            M_complex[ii,jj] = complex(M[i,j], M[i,j+1])
            jj += 1
        end
        jj = 1
        ii += 1
    end

    return M_complex
end

"""
    round_complex_values(M::Array{Complex{Float64},2})

Given a complex-valued 2D matrix, round_complex_values function returns a complex-valued matrix which 
rounds the valuest closest to 0 and 1. This is useful to avoid numerical issues. 
"""
function round_complex_values(M::Array{Complex{Float64},2})
    # round values close to 0 (within toleranes) for both real and imaginary values
    # Input can be a vector (>= 1 element) or a matrix of complex values
   
    if size(M)[1] == 0
        Memento.error(_LOGGER, "Input cannot be a scalar")
    end

    if (length(size(M)) == 2)
        n_r = size(M)[1]
        n_c = size(M)[2]

        M_round = Array{Complex{Float64},2}(zeros(n_r,n_c))
        
        for i=1:n_r
            for j=1:n_c 

                # Real components
                if isapprox(real(M[i,j]), 0, atol=1E-6)
                    Mij_r = 0
                elseif isapprox(real(M[i,j]), 1, atol=1E-6)
                    Mij_r = 1
                else 
                    Mij_r = real(M[i,j])
                end
                
                # Imaginary components
                if isapprox(imag(M[i,j]), 0, atol=1E-6)
                    Mij_i = 0
                elseif isapprox(imag(M[i,j]), 1, atol=1E-6)
                    Mij_i = 1
                else 
                    Mij_i = imag(M[i,j])
                end

                M_round[i,j] = Mij_r + (Mij_i)im
                
            end
        end
    end
    
    return M_round
end

"""
    unique_idx(x::AbstractArray{T})

unique_idx returns the indices of unique elements in a given array of scalar or vector inputs. Overall, 
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

unique_matrices returns the unique set of matrices and the corresponding indices 
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
    kron_single_gate(num_qubits::Int64, M::Array{Complex{Float64},2}, qubit_loc::String)

Given number of qubits of the circuit, the complex-valued gate and the qubit location ("q1","q2',"q3",...),
kron_single_gate function returns a full-sized gate after applying appropriate kronecker products. 
"""
function kron_single_gate(num_qubits::Int64, M::Array{Complex{Float64},2}, qubit_loc::String)
    
    if size(M)[1] >= 2^num_qubits
        Memento.warn(_LOGGER, "Input gate is already in $num_qubits qubits")
        return M
    end

    I = QCO.IGate(1)

    if num_qubits == 2     

        if qubit_loc == "q1" 
            return kron(M,I)
        elseif qubit_loc == "q2"
            return kron(I,M)
        else
            Memento.error(_LOGGER, "For num_qubits = $num_qubits, qubit location has to be ∈ [q1, q2]") 
        end

    elseif num_qubits == 3 
        
        if qubit_loc == "q1" 
            return kron(kron(M,I),I)
        elseif qubit_loc == "q2" 
            return kron(kron(I,M),I)
        elseif qubit_loc == "q3" 
            return kron(kron(I,I),M)
        else
            Memento.error(_LOGGER, "For num_qubits = $num_qubits, qubit location has to be ∈ [q1, q2, q3]")
        end

    elseif num_qubits == 4 
        
        if qubit_loc == "q1" 
            return kron(kron(kron(M,I),I),I)
        elseif qubit_loc == "q2" 
            return kron(kron(kron(I,M),I),I)
        elseif qubit_loc == "q3" 
            return kron(kron(kron(I,I),M),I)
        elseif qubit_loc == "q4" 
            return kron(kron(kron(I,I),I),M)
        else
            Memento.error(_LOGGER, "For num_qubits = $num_qubits, qubit location has to be ∈ [q1, q2, q3, q4]")
        end
    
    # Larger qubit circuits can be supported here.

    end

end