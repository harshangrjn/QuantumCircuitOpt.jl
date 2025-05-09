##############################################################################
#  Imports and lightweight fallbacks
##############################################################################
using Printf
import QuantumCircuitOpt as QCO
import LinearAlgebra as LA

##############################################################################
#  Basic data model
##############################################################################
# Gate struct for post-optimization
struct Gate
    label  :: String                       
    qubits :: BitSet                       
    mat    :: Matrix{ComplexF64}  # full 2^n × 2^n unitary
end


##############################################################################
#  String helpers: single gate, plus “x”‐tensor expansion
##############################################################################
function parse_gate(s::AbstractString)::Gate
    head, params = occursin('(', s) ? split(s, '(', limit = 2) : (s,"")
    parts = split(head, '_')
    op    = parts[1]
    qs    = BitSet(parse.(Int, parts[2:end]))
    label = params=="" ? op : op * "(" * params
    Gate(label, qs, Matrix{ComplexF64}(I,1,1))   # dummy mat (overwritten later)
end

explode_kron(seq::Vector{String}) = [parse_gate(term) for item in seq for term in split(item,'x')]


##############################################################################
#  Greedy ASAP / slide-left scheduler  (gates, n_qubits) → layers
##############################################################################
function compress_layers(gates::Vector{Gate}, n::Int)
    layers, masks = Vector{Vector{Gate}}(), BitSet[]
    full = BitSet(1:n)

    for g in gates
        @assert g.qubits ⊆ full "gate $(g.label) uses invalid qubit (n=$n)"
        pos = length(layers)+1
        for l = length(layers):-1:1
            if !isempty(intersect(g.qubits, masks[l])); break; end
            pos = l
        end
        if pos > length(layers)
            push!(layers, [g]); push!(masks, copy(g.qubits))
        else
            push!(layers[pos], g); union!(masks[pos], g.qubits)
        end
    end
    return layers
end


##############################################################################
#  Kronecker-style pretty printer
##############################################################################
function kron_layer(layer::Vector{Gate}, n::Int)::String
    sorted = sort(layer; by = g -> minimum(g.qubits))
    facs, qptr = String[], 1
    for g in sorted
        firstq, lastq = minimum(g.qubits), maximum(g.qubits)
        for _ in qptr:firstq-1; push!(facs,"I"); end
        if length(g.qubits)==1
            push!(facs, g.label)
        else
            push!(facs, @sprintf("%s_{%s}", g.label,
                                 join(sort(collect(g.qubits)),",")))
        end
        qptr = lastq+1
    end
    for _ in qptr:n; push!(facs,"I"); end
    return join(facs, " ⊗ ")
end


##############################################################################
#  Multiply all layers (full matrices already provided)
##############################################################################
function circuit_unitary(layers::Vector{Vector{Gate}})
    dim = size(layers[1][1].mat,1)        # 2^n
    U = Matrix{ComplexF64}(LA.I, dim, dim)
    for layer in layers
        G = Matrix{ComplexF64}(LA.I, dim, dim)
        for g in layer
            G = g.mat * G                 # gates in same layer commute (disjoint qubits)
        end
        U = G * U
    end
    U
end


##############################################################################
#  Build layers from `id_sequence` + `gates_dict`
##############################################################################
function build_layers(id_seq::Vector{Int}, gdict::Dict{String,Any})
    gates = Gate[]
    for id in id_seq
        info   = gdict[string(id)]
        typestr = info["type"][1]
        typestr == "Identity" && continue          # skip

        head, params = occursin('(', typestr) ? split(typestr,'(',limit=2) : (typestr,"")
        parts = split(head, '_'); op = parts[1]
        qs    = BitSet(parse.(Int, parts[2:end]))
        mat   = convert(Array{ComplexF64,2}, info["matrix"])

        # round angles if present
        label = op
        if haskey(info,"angle")
            a = info["angle"]; order = ["θ","ϕ","λ"]
            vals = [round(rad2deg(a[k]), digits=3) for k in order if haskey(a,k)]
            label *= "(" * join(vals,",") * ")"
        elseif params!=""; label *= "(" * params; end

        push!(gates, Gate(label, qs, mat))
    end

    n_qubits = round(Int, log2(size(gates[1].mat,1)))
    layers_before = [[g] for g in gates]
    layers_after  = compress_layers(gates, n_qubits)
    return layers_before, layers_after
end


##############################################################################
#  Validation helper
##############################################################################
function validate_circuit(data::Dict{String,Any},
                          layers::Vector{Vector{Gate}};
                          atol = 1e-4,
                          error_message = true)

    U_sol = circuit_unitary(layers)

    tgt = data["are_gates_real"] ?
          real(data["target_gate"]) :
          QCO.real_to_complex_gate(data["target_gate"])

    target_gate = convert(Array{ComplexF64,2}, tgt)

    valid = if data["decomposition_type"] in ["exact_optimal", "exact_feasible"]
        QCO.isapprox(U_sol, target_gate; atol=atol)
        @show U_sol
        @show target_gate
    elseif data["decomposition_type"] == "optimal_global_phase"
        QCO.isapprox_global_phase(U_sol, target_gate)
    else
        false
    end

    if !valid && error_message
        # Memento.error(_LOGGER,
        #     "Decomposition is not valid: Problem may be infeasible")
    end
    return valid
end


##############################################################################
#  ---------------------------  DEMO  ----------------------------------------
##############################################################################

include("2qubit_gates.jl")

params = minimize_T_gate()
data = QCO.get_data(params)

id_sequence = [8,1,2,3]

# ----------------------------------------------------------------------------
layers_before, layers_after = build_layers(id_sequence, data["gates_dict"])

println("=== BEFORE compression ===")
for (d,lay) in enumerate(layers_before)
    println("depth $d : ", kron_layer(lay, 2))
end

println("\n=== AFTER compression ===")
for (d,lay) in enumerate(layers_after)
    println("depth $d : ", kron_layer(lay, 2))
end

# ----------------------------------------------------------------------------
println("\nValidation:")
println("  before = ", validate_circuit(data, layers_before))
println("  after  = ", validate_circuit(data, layers_after))
