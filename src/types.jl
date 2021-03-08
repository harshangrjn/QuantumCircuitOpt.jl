mutable struct OptimizerOptions
    log_level :: Int
    time_limit :: Float64
    tol :: Float64
    silent :: Bool
    binary_relax :: Bool
    mip_gap :: Float64

    mip_solver :: Any 
    nlp_solver :: Any 
    minlp_solver :: Any
end

function get_default_options()
    log_level = 1
    time_limit = 10800
    tol = 1E-6
    silent = true
    binary_relax = false
    mip_gap = 1E-4

    mip_solver = nothing
    nlp_solver = nothing 
    minlp_solver = nothing

    return OptimizerOptions(log_level, time_limit, tol, silent, binary_relax, 
                            mip_gap, mip_solver, nlp_solver, minlp_solver)
end

mutable struct QCOoptimizer <: MOI.AbstractOptimizer
    options::OptimizerOptions    
    Qmodel :: JuMP.Model
    objval :: Float64
    best_bound :: Float64
    objective :: Union{Nothing, MOI.ScalarAffineFunction{Float64}, MOI.ScalarQuadraticFunction{Float64}}
    solution :: Any
    sol_time :: Float64

    mip_solver :: Any 
    nlp_solver :: Any 
    minlp_solver :: Any 

    Qmodel_status::MOI.TerminationStatusCode

    # constructor
    function QCOoptimizer()
        m = new()
        m.options = get_default_options()
        MOI.empty!(m)
        return m
    end
end

mutable struct QCOdata 
    Q_gates :: Array{Complex,3}
    target_gate :: Array{Complex,2}
    n_r :: Int
    n_c :: Int
    Q_gates_lb :: Array{Complex,2}
    Q_gates_ub :: Array{Complex,2}
end