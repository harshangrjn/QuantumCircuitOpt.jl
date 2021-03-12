using QuantumCircuitOpt
using Memento
using JuMP
using LinearAlgebra
using Test
using Cbc

const QCO = QuantumCircuitOpt

# Suppress warnings during testing
QCO.logger_config!("error")

const CBC = JuMP.optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)

@testset "QuantumCircuitOpt" begin

    include("utility_tests.jl")
    include("data_tests.jl")

end