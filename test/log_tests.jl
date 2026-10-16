@testset "Logging tests" begin
    io = IOBuffer()
    QCO._LOGGER[] =
        Logging.ConsoleLogger(io, Logging.Debug; meta_formatter = QCO._meta_formatter)

    QCO.@_debug("debug message")
    QCO.@_info("info message")
    QCO.@_warn("warning message")
    @test_throws ErrorException QCO.@_error("error message")

    output = String(take!(io))
    @test occursin("debug message", output)
    @test occursin("info message", output)
    @test occursin("warning message", output)
    @test occursin("error message", output)

    QCO.logger_config!(Logging.Error)
end
