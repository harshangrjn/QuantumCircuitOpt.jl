export QuantumCircuitModel

"""
The composite mutable struct, `QuantumCircuitModel`, holds dictionaries for input data, abstract JuMP model for optimization,
variable references and result from solving the JuMP model.
"""
mutable struct QuantumCircuitModel 
    data::Dict{String,Any}
    model::JuMP.Model
    variables::Dict{Symbol,Any}
    #constraints::Dict{Symbol,Any}
    result::Dict{String,Any}
end

"Contructor for struct `QuantumCircuitModel`"
function QuantumCircuitModel(data)
    qcm = new(data)
    qcm.data = data
    qcm.model = JuMP.Model() 
    qcm.variables = Dict{Symbol,Any}()
    qcm.result = Dict{String,Any}()
    return qcm
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
