export QuantumCircuitModel

"""
    QCModelOptions
The composite mutable struct, `QCModelOptions`, holds various optimization model options for enhancements 
with defualt options set to the values provided by `get_default_options` function.
"""
mutable struct QCModelOptions
    model_type                         :: String
    
    all_valid_constraints              :: Int64
    commute_gate_constraints           :: Bool
    involutory_gate_constraints        :: Bool
    redundant_gate_pair_constraints    :: Bool 
    identity_gate_symmetry_constraints :: Bool
    visualize_solution                 :: Bool
    idempotent_gate_constraints        :: Bool
    convex_hull_gate_constraints       :: Bool

    time_limit                         :: Float64
    relax_integrality                  :: Bool
    optimizer_log                      :: Bool
    objective_slack_penalty            :: Float64
end

"""
    get_default_options()
This function returns the default options for building the struct `QCModelOptions`.
"""
function get_default_options()
    model_type                         = "compact_formulation" # compact_formulation, balas_formulation
    
    all_valid_constraints              = 0    # -1, 0, 1
    commute_gate_constraints           = true # true, false
    involutory_gate_constraints        = true # true, false
    redundant_gate_pair_constraints    = true # true, false
    identity_gate_symmetry_constraints = true # true, false
    visualize_solution                 = true # true, false

    idempotent_gate_constraints        = false # true, false
    convex_hull_gate_constraints       = false # true, false

    time_limit                         = 10800 # float value
    relax_integrality                  = false # true, false
    optimizer_log                      = true  # true, false
    objective_slack_penalty            = 1E3   # > 0 value

    return QCModelOptions(model_type,
                          all_valid_constraints,
                          commute_gate_constraints,
                          involutory_gate_constraints,
                          redundant_gate_pair_constraints,
                          identity_gate_symmetry_constraints,
                          visualize_solution,
                          idempotent_gate_constraints,
                          convex_hull_gate_constraints,
                          time_limit,
                          relax_integrality,
                          optimizer_log,
                          objective_slack_penalty)
end

"""
    QuantumCircuitModel
The composite mutable struct, `QuantumCircuitModel`, holds dictionaries for input data, abstract JuMP model for optimization,
variable references and result from solving the JuMP model.
"""
mutable struct QuantumCircuitModel 
    data      :: Dict{String,Any}
    model     :: JuMP.Model
    options   :: QCModelOptions
    variables :: Dict{Symbol,Any}
    result    :: Dict{String,Any}

    "Constructor for struct `QuantumCircuitModel`"
    function QuantumCircuitModel(data::Dict{String,Any})     
        data      = data
        model     = JuMP.Model() 
        options   = QCO.get_default_options()
        variables = Dict{Symbol,Any}()
        result    = Dict{String,Any}()
        
        qcm = new(data, model, options, variables, result)

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
    type    :: String
    complex :: Array{Complex{Float64},2}
    real    :: Array{Float64,2}
    inverse :: Array{Float64,2}
    isreal  :: Bool

    "Contructor for struct `Gate`"
    function GateData(gate_type::String, num_qubits::Int64)
        type    = gate_type
        complex = QCO.get_full_sized_gate(type, num_qubits)
        real    = QCO.complex_to_real_gate(complex)
        inverse = inv(real)
        isreal  = iszero(imag(complex))
        
        gate = new(type, complex, real, inverse, isreal)

        return gate
    end

end
