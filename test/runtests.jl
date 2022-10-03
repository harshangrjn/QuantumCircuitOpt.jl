using JuMP
using Test 

import QuantumCircuitOpt as QCO
import Memento
import LinearAlgebra as LA
import Cbc
import HiGHS

# Suppress warnings during testing
QCO.logger_config!("error")

const HIGHS = MOI.OptimizerWithAttributes(
    HiGHS.Optimizer,
    "presolve" => "on",
    "log_to_console" => false,
)

const MIP_SOLVER = HIGHS

tol_0 = 1E-6

@testset "QuantumCircuitOpt" begin

    include("utility_tests.jl")
    include("chull_tests.jl")
    include("gates_tests.jl")
    include("types_tests.jl")
    include("data_tests.jl")
    include("qc_model_tests.jl")
    include("relaxations_tests.jl")

end