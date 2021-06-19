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

function get_commutative_gates(M::Array{Float64,3})

    depth = size(M)[3]

    M_commute_2 = Array{Tuple{Int64,Int64},1}()
    M_commute_3 = Array{Tuple{Int64,Int64,Int64},1}()

    # Check for commutative matrix pairs
    for d1 = 1:depth
        for d2 = (d1+1):depth
            M_d1 = real_to_complex_matrix(M[:,:,d1])
            M_d2 = real_to_complex_matrix(M[:,:,d2])

            if isapprox(M_d1*M_d2, M_d2*M_d1, atol = 1E-4)
                push!(M_commute_2, (d1, d2))
            end

            for d3 = (d2+1):depth
                M_d3 = real_to_complex_matrix(M[:,:,d3])

                M_product = Array{Complex{Float64},3}(zeros(size(M_d3)[1], size(M_d3)[2], 6))

                # All commuting permutations
                M_product[:,:,1] = M_d1 * M_d2 * M_d3
                M_product[:,:,2] = M_d2 * M_d1 * M_d3
                M_product[:,:,3] = M_d1 * M_d3 * M_d2
                M_product[:,:,4] = M_d3 * M_d1 * M_d2
                M_product[:,:,5] = M_d2 * M_d3 * M_d1
                M_product[:,:,6] = M_d3 * M_d2 * M_d1
                
                triplet_commutes = true

                for i=1:6
                    for j=(i+1):6
                        if !isapprox(M_product[:,:,i], M_product[:,:,j], atol = 1E-4)
                            triplet_commutes = false
                        end
                    end
                end

                if triplet_commutes
                    push!(M_commute_3, (d1, d2, d3))
                end

            end
        end
    end
 
    return M_commute_2, M_commute_3
end

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

function real_to_complex_matrix(M::Array{Float64,2})
    
    n = size(M)[1]
    if !iseven(n)
        Memento.error(_LOGGER, "Input real matrix can admit only even numbered columns and rows")
    end
    
    M_complex = zeros(Complex{Float64}, (Int(n/2), Int(n/2)))
  
    ii = 1; jj = 1;
    for i = collect(1:2:n)
        for j = collect(1:2:n)

            if !isapprox(M[i,j], M[i+1, j+1], atol = 1E-4) || !isapprox(M[i+1,j], -M[i,j+1], atol = 1E-4)
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

function round_complex_values(M::Array{Complex{Float64},2})
    # round values close to 0 and 1 (within toleranes) for both real and imaginary values
    # Input can be a vector (>= 1 element) or a matrix of complex values
   
    if size(M)[1] == 0
        Memento.error(_LOGGER, "Input cannot be a scalar")
    end

    tol = 1E-6
    n_r = size(M)[1]

    if (length(size(M)) == 1)
        M_round = Array{Complex{Float64},2}(zeros(n_r))
        for i=1:n_r
            M_r = 0.0; M_i = 0.0;
            if abs(real(M[i])) > tol
                M_r = real(M)
            end
            if abs(imag(M[i])) > tol
                M_i = imag(M[i])
            end
            if (abs(M_r) > 0) || (abs(M_i) > 0)
                M_round[i] = M_r + (M_i)im
            end    
        end
    end

    if (length(size(M)) == 2)
        n_c = size(M)[2]
        M_round = Array{Complex{Float64},2}(zeros(n_r,n_c))
        for i=1:n_r
            for j=1:n_c 
                Mij_r = 0.0; Mij_i = 0.0;
                if abs(real(M[i,j])) > tol
                    Mij_r = real(M[i,j])
                end
                if abs(imag(M[i,j])) > tol
                    Mij_i = imag(M[i,j])
                end
                if (abs(Mij_r) > 0) || (abs(Mij_i) > 0)
                    M_round[i,j] = Mij_r + (Mij_i)im
                end    
            end
        end
    end

    if (length(size(M)) == 3)
        n_c = size(M)[2]
        n_d = size(M)[3]
        M_round = Array{Complex{Float64},3}(zeros(n_r, n_c, n_d))
        for i=1:n_r
            for j=1:n_c 
                for k=1:n_d
                    Mijk_r = 0.0; Mijk_i = 0.0;
                    if abs(real(M[i,j,k])) > tol
                        Mijk_r = real(M[i,j,k])
                    end
                    if abs(imag(M[i,j])) > tol
                        Mijk_i = imag(M[i,j,k])
                    end
                    if (abs(Mijk_r) > 0) || (abs(Mijk_i) > 0)
                        M_round[i,j,k] = Mijk_r + (Mijk_i)im
                    end    
                end
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

unique_matrices returns the unique set of matrices and the corresponding indices of unique matrices from the given set of matrices.  
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
