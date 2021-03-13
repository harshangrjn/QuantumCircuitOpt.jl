# Unit tests for functions in utility.jl

M_c = [complex(1,0)       complex(0,-1)
       complex(-0.5,0.5)  complex(0,0)]

M_r = [1.0   0.0   0.0  1.0
       0.0   1.0  -1.0  0.0
      -0.5  -0.5   0.0  0.0
       0.5  -0.5   0.0  0.0]

# Variable bounds to test get_auxiliary_var_bounds function
l = [-2, -1, 0, 2]
u = [-1,  2, 3, 2]

@testset "complex to real matrix function tests" begin
    test_M_r = QCO.get_complex_to_real_matrix(M_c)
    @test test_M_r == M_r
end

@testset "real to complex matrix function tests" begin
    test_M_c = QCO.get_real_to_complex_matrix(M_r)
    @test test_M_c == M_c
end

@testset "auxiliary variable product bounds tests" begin
    m = Model()
    @variable(m, l[i] <= x[i=1:length(l)] <= u[i])
    # x1*x2 
    @test QCO.get_auxiliary_var_bounds([x[1], x[2]]) == (-4, 2)
    # x1*x3 
    @test QCO.get_auxiliary_var_bounds([x[1], x[3]]) == (-6, 0)
    # x1*x4
    @test QCO.get_auxiliary_var_bounds([x[1], x[4]]) == (-4, -2)
    # x4*x4
    @test QCO.get_auxiliary_var_bounds([x[4], x[4]]) == (4, 4)
    # x1*x2*x3 
    @test QCO.get_auxiliary_var_bounds([x[1], x[2], x[3]]) == (-12, 6)
    @test QCO.get_auxiliary_var_bounds([x[3], x[1], x[2]]) == (-12, 6)
    @test QCO.get_auxiliary_var_bounds([x[2], x[1], x[3]]) == (-12, 6)
    # x1*x2*x3*x4
    @test QCO.get_auxiliary_var_bounds([x[1], x[2], x[3], x[4]]) == (-24, 12)
    @test QCO.get_auxiliary_var_bounds([x[4], x[1], x[3], x[2]]) == (-24, 12)
    @test QCO.get_auxiliary_var_bounds([x[4], x[3], x[2], x[1]]) == (-24, 12)
end

#TODO: Add test to verify verify_tolerances_complex_matrix function