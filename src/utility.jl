function get_auxiliary_var_bounds(v::Array{JuMP.VariableRef,1}) 
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
 
function get_matrix_bounds(M) 
    tol_0 = 1E-6
    n_r = length(M[:,1,1])
    n_c = length(M[1,:,1])
    M_l = zeros(n_r, n_c)
    M_u = zeros(n_r, n_c)
    for i = 1:n_r
        for j = 1:n_c
            M_l[i,j] = minimum(M[i,j,:])
            M_u[i,j] = maximum(M[i,j,:])
        end
    end
    
    k = 0
    for i = 1:n_c
        for j = 1:n_r
            @assert M_l[i,j] <= M_u[i,j]
            ((abs(M_l[i,j] - M_u[i,j])) <= tol_0) && (k += 1)
        end
    end
    (k >= 1) && (println(">>> $k (of $(n_c * n_r)) number of entries in the given set of matrices have equal lower and upper bounds."))
    return M_l, M_u
    
end

function get_commutative_matrices(M)
    tol_0 = 1E-6
    n_r = length(M[:,1,1])
    n_c = length(M[1,:,1])
    n = length(M[1,1,:])
    Id = Matrix(I, n_r, n_c)
    Id_idx = [];
    for i = 1:n
        (all(isapprox.(M[:,:,i], Id))) && (push!(Id_idx, i))
    end
    M_idx = setdiff(1:n, Id_idx)
    M_commute = Any[];
    # Check for pairs of commutative matrices
    for i in M_idx
        for j in M_idx
            if (i < j) && (i != j)
             # Evaluate the commutator 
                (all(isapprox.(M[:,:,i] * M[:,:,j], M[:,:,j] * M[:,:,i]))) && push!(M_commute, (i, j))
            end
        end
    end
 
    return M_commute, Id_idx
end

function complex_to_real_matrix(M)
    @assert ((typeof(M) == Array{Complex{Int64},2}) || (typeof(M) == Array{Complex{Float64},2}))
    n = size(M)[1]
    M_real = zeros(2*n, 2*n)
  
    ii = 1; jj = 1;
    for i = collect(1:2:2 * n)
        for j = collect(1:2:2 * n)
            M_real[i,j] = real(M[ii,jj])
            M_real[i+1,j+1] = real(M[ii,jj])
            if imag(M[ii,jj]) == 0
                M_real[i,j+1] = 0
                M_real[i+1,j] = 0
            else
                M_real[i,j+1] = -imag(M[ii,jj])
                M_real[i+1,j] = imag(M[ii,jj])
            end
            jj += 1
        end
        jj = 1
        ii += 1
    end
  
    return M_real
end

function real_to_complex_matrix(M)
    @assert ((typeof(M) == Array{Int64,2}) || (typeof(M) == Array{Float64,2}))
    n = size(M)[1]
    @assert iseven(n)
    M_complex = zeros(Complex{Float64}, (Int(n/2), Int(n/2)))
  
    ii = 1; jj = 1;
    for i = collect(1:2:n)
        for j = collect(1:2:n)
            @assert M[i,j] == M[i+1,j+1]
            @assert M[i+1,j] == -M[i,j+1]
            M_complex[ii,jj] = complex(M[i,j], M[i+1,j])
            jj += 1
        end
        jj = 1
        ii += 1
    end

    return M_complex
end

function verify_tolerances_complex_values(M::Array{Complex{Float64},2})
    # round values close to 0 and 1 (within toleranes) for both real and imaginary values
    # Input can be a vector (>= 1 element) or a matrix of complex values
    @assert size(M)[1] != 0
    tol = 1E-6
    n_r = size(M)[1]
    if (length(size(M)) == 1)
        M_round = zeros(Complex, n_r)
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
        M_round = zeros(Complex, n_r, n_c)
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
        M_round = zeros(Complex, n_r, n_c, n_d)
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