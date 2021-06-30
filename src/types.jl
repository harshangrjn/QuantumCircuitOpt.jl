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
The composite mutable struct, `GateData`, type of the gate, the complex matrix form 
    of the gate, full sized real form of the gate, inverse of the gate and a boolean which
    states if the gate has all real entries.
"""
mutable struct GateData 
    type::String
    complex::Array{Complex{Float64},2}
    real::Array{Float64,2}
    inverse::Array{Float64,2}
    isreal::Bool

    "Contructor for struct `Gate`"
    function GateData(gate_type::String, num_qubits::Int64)
        type = gate_type
        complex = QCO.get_full_sized_gate(type, num_qubits)
        real = QCO.complex_to_real_matrix(complex)
        inverse = inv(real)
        isreal = iszero(imag(complex))
        gate = new(type, complex, real, inverse, isreal)

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
