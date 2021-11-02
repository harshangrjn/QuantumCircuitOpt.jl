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
        real = QCO.complex_to_real_gate(complex)
        inverse = inv(real)
        isreal = iszero(imag(complex))
        gate = new(type, complex, real, inverse, isreal)

        return gate
    end

end