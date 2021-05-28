export QuantumCircuitModel

"""
    QuantumCircuitModel
The composite mutable struct, `QuantumCircuitModel`, holds dictionaries for input data, abstract JuMP model for optimization,
variable references and result from solving the JuMP model.
"""
mutable struct QuantumCircuitModel 
    data::Dict{String,Any}
    model::JuMP.Model
    variables::Dict{Symbol,Any}
    #constraints::Dict{Symbol,Any}
    result::Dict{String,Any}

    "Contructor for struct `QuantumCircuitModel`"
    function QuantumCircuitModel(data::Dict{String,Any})
        
        data = data
        model = JuMP.Model() 
        variables = Dict{Symbol,Any}()
        result = Dict{String,Any}()
        qcm = new(data, model, variables, result)

        return qcm
    end

end

"""
    GateData
"""
mutable struct GateData 
    type::String
    complex::Array{Complex{Float64},2}
    real::Array{Float64,2}
    inverse::Array{Float64,2}
    isallreal::Bool

    "Contructor for struct `Gate`"
    function Gate(gate_type::String, n_qubits::Int64)
        type = gate_type

        complex  = getfield(Main, Symbol(gate_type))(n_qubits) # String to function name 
        real = QCO.get_complex_to_real_matrix(complex)
        inverse = inv(real)
        isallreal = false
        gate = new(type, complex, real, inverse, isallreal)

        return gate
    end

end

# mutable struct OptimizerOptions
#     log_level :: Int
#     time_limit :: Float64
#     tol :: Float64
#     silent :: Bool
#     binary_relax :: Bool
#     mip_gap :: Float64

#     mip_solver :: Any 
#     nlp_solver :: Any 
#     minlp_solver :: Any
# end

# function get_default_options()
#     log_level = 1
#     time_limit = 10800
#     tol = 1E-6
#     silent = true
#     binary_relax = false
#     mip_gap = 1E-4

#     mip_solver = nothing
#     nlp_solver = nothing 
#     minlp_solver = nothing

#     return OptimizerOptions(log_level, time_limit, tol, silent, binary_relax, 
#                             mip_gap, mip_solver, nlp_solver, minlp_solver)
# end
