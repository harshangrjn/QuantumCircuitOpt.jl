M_c = [complex(1,0)       complex(0,-1)
       complex(-0.5,0.5)  complex(0,0)]

M_r = [1.0   0.0   0.0  1.0
       0.0   1.0  -1.0  0.0
      -0.5  -0.5   0.0  0.0
       0.5  -0.5   0.0  0.0]

@testset "complex to real matrix function tests" begin
    test_M_r = QCO.complex_to_real_matrix(M_c)
    @test test_M_r == M_r
end

@testset "real to complex matrix function tests" begin
    test_M_c = QCO.real_to_complex_matrix(M_r)
    @test test_M_c == M_c
end