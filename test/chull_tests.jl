@testset "Tests: convex_hull" begin
    points_1a = [(0.5, 0.0), (0.0, 0.0), (1, 0), (0,1), (-0.5, 0), (0.0, 0.5)]
    points_1b = [(0.0, 0.0), (1, 0), (0,1), (-0.5, 0), (0.5, 0.0), (0.0, 0.5)]
    points_2a = [(0.5, 0.0), (0.0, 0.0), (1, 0), (0,1), (-0.5, 0), (0.0, 0.5), (0.5, -0.5)]
    points_2b = [(0.5, 0.0), (0.5, -0.5), (0.0, 0.0), (0,1), (-0.5, 0), (0.0, 0.5), (1, 0)]

    points_3a = [(0,0)]
    points_3b = [(-0.5,-0.5), (0.5,0.5)]
    points_3c = [(-0.5, -0.5), (0.5, 0.5), (0, 0), (0.25, 0.25)]

    @test QCO.convex_hull(points_1a) == QCO.convex_hull(points_1b)
    @test QCO.convex_hull(points_2a) == QCO.convex_hull(points_2b)
    @test QCO.convex_hull(points_3a) == points_3a
    @test QCO.convex_hull(points_3b) == points_3b
    @test QCO.convex_hull(points_3c) == points_3b
    
end