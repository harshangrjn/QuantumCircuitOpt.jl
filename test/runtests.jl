using QuantumCircuitOptimization
import Memento
using JuMP
import LinearAlgebra

const QCO = QuantumCircuitOptimization

# Suppress warnings during testing
QCO.logger_config!("error")

using Test
using Cbc

const CBC = JuMP.optimizer_with_attributes(Cbc.Optimizer, "logLevel" => 0)

@testset "QuantumCircuitOptimization" begin

    include("utility_tests.jl")

end