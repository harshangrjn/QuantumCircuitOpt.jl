# Unit tests for functions in utility.jl

M_c = [complex(1,0)       complex(0,-1)
       complex(-0.5,0.5)  complex(0,0)]

M_r = [1.0   0.0   0.0  -1.0
       0.0   1.0   1.0  0.0
      -0.5  0.5    0.0  0.0
      -0.5  -0.5   0.0  0.0]

# Variable bounds to test auxiliary_variable_bounds function
l = [-2, -1, 0, 2]
u = [-1,  2, 3, 2]

@testset "complex to real matrix function tests" begin
    test_M_r = QCO.complex_to_real_matrix(M_c)
    @test test_M_r == M_r
end

@testset "real to complex matrix function tests" begin
    test_M_c = QCO.real_to_complex_matrix(M_r)
    @test test_M_c == M_c
end

@testset "auxiliary variable product bounds tests" begin
    m = Model()
    @variable(m, l[i] <= x[i=1:length(l)] <= u[i])
    # x1*x2 
    @test QCO.auxiliary_variable_bounds([x[1], x[2]]) == (-4, 2)
    # x1*x3 
    @test QCO.auxiliary_variable_bounds([x[1], x[3]]) == (-6, 0)
    # x1*x4
    @test QCO.auxiliary_variable_bounds([x[1], x[4]]) == (-4, -2)
    # x4*x4
    @test QCO.auxiliary_variable_bounds([x[4], x[4]]) == (4, 4)
    # x1*x2*x3 
    @test QCO.auxiliary_variable_bounds([x[1], x[2], x[3]]) == (-12, 6)
    @test QCO.auxiliary_variable_bounds([x[3], x[1], x[2]]) == (-12, 6)
    @test QCO.auxiliary_variable_bounds([x[2], x[1], x[3]]) == (-12, 6)
    # x1*x2*x3*x4
    @test QCO.auxiliary_variable_bounds([x[1], x[2], x[3], x[4]]) == (-24, 12)
    @test QCO.auxiliary_variable_bounds([x[4], x[1], x[3], x[2]]) == (-24, 12)
    @test QCO.auxiliary_variable_bounds([x[4], x[3], x[2], x[1]]) == (-24, 12)
end

@testset "unique_matrices and unique_idx tests" begin
    
    tol = 1E-6
    D = zeros(3,3,4)
    D[:,:,1] = Matrix(LA.I,3,3)
    D[:,:,2] = rand(3,3)
    D[:,:,3] = Matrix(LA.I,3,3)
    D[:,:,4] = rand(3,3)

    D[1,1,3] = 1 + 2*tol
    D[3,3,3] = 1 - tol
    D[1,2,3] = 0 - tol
    D[1,3,3] = 0 + tol

    D[isapprox.(D, 0, atol=1E-6)] .= 0

    D_unique, D_unique_idx = QCO.unique_matrices(D)

    @test D_unique_idx == [1,2,4]
    @test D_unique[:,:,1] == D[:,:,1]
    @test D_unique[:,:,2] == D[:,:,2]
    @test D_unique[:,:,3] == D[:,:,4]

end

#TODO: Add test to verify verify_tolerances_complex_matrix function