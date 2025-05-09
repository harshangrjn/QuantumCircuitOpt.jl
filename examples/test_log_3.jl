###############################################################################
# Depth-compression for a gate list carried in `id_seq` / `gates_dict`
###############################################################################

"""
Gate  ─ A minimal record the scheduler needs and nothing more.
`label`  – for humans/debugging, e.g. `"3:U3_1"`  
`qubits` – **BitSet** of the qubit indices this gate touches  
`mat`    – the full unitary kept for any later optimisation or simulation
"""
struct Gate
    label  :: String
    qubits :: BitSet
    mat    :: Matrix{ComplexF64}
end


"""
    build_circuit_layers(id_seq, gates_dict)
           -> (cirq_precompress, cirq_postcompress)

* **`id_seq`**        – Ordered vector of integer keys that define the *original*
                        (uncompressed) execution order.
* **`gates_dict`**    – Dictionary exactly as supplied (`Dict{String,Any}`).

Returns two vectors of layers, each layer being a `Vector{Gate}`:

1. **`cirq_precompress`** – one gate per layer in the incoming order.
2. **`cirq_postcompress`** – layers produced by a greedy
   graph-colouring heuristic where no qubit appears in more than one gate
   per layer.
"""
function build_circuit_layers(id_seq::Vector{Int},
                              gates_dict::Dict{String,Any})

    # ————————————————————————————————————————————————————————————————
    # helper: extract every integer that appears in a string
    # e.g. "CNot_2_3"  →  [2,3]   |   "qubit_1" → [1]
    # ————————————————————————————————————————————————————————————————
    get_ints(s::AbstractString) =
        parse.(Int, collect(m.match for m in eachmatch(r"\d+", s)))

    # ————————————————————————————————————————————————————————————————
    # 1. Build the Gate objects in the *given* order
    # ————————————————————————————————————————————————————————————————
    gates = Gate[]
    for id in id_seq
        gd      = gates_dict[string(id)]
        # ─ which symbolic gate name shall we log?
        tfield  = gd["type"]
        tvec    = isa(tfield, AbstractVector) ? tfield : [tfield]
        typstr  = isempty(tvec) ? "Identity" :
                  (findfirst(!=("Identity"), tvec) === nothing ?
                   "Identity" : tvec[findfirst(!=("Identity"), tvec)])

        # ─ which qubits does it touch?
        qubits  = if haskey(gd, "qubit_loc")
                     get_ints(gd["qubit_loc"])            # single-qubit form
                  else
                     get_ints(typstr)                     # parse the type name
                  end

        push!(gates, Gate(
            "$(id):$typstr",            # label
            BitSet(qubits),             # active qubits
            gd["matrix"],               # full unitary
        ))
    end

    # uncompressed circuit: one gate per depth
    cirq_precompress = [[g] for g in gates]

    # ————————————————————————————————————————————————————————————————
    # 2. Build the undirected conflict adjacency list
    #    edge ⇔ the two gates share ≥ 1 qubit
    # ————————————————————————————————————————————————————————————————
    n   = length(gates)
    adj = [Int[] for _ in 1:n]
    for i in 1:n-1, j in i+1:n
        if !isempty(intersect(gates[i].qubits, gates[j].qubits))
            push!(adj[i], j);  push!(adj[j], i)
        end
    end

    # ————————————————————————————————————————————————————————————————
    # 3. Greedy DSATUR-style colouring  → depth labels
    # ————————————————————————————————————————————————————————————————
    order = sortperm(length.(adj); rev = true)     # visit high-degree first
    depth = zeros(Int, n)

    for v in order
        forbidden = depth[adj[v]]
        d = 1
        while d in forbidden; d += 1; end
        depth[v] = d
    end

    # ————————————————————————————————————————————————————————————————
    # 4. Collect gates by their colour (= depth)
    # ————————————————————————————————————————————————————————————————
    maxd              = maximum(depth)
    cirq_postcompress = [Gate[] for _ in 1:maxd]
    for (g, d) in zip(gates, depth)
        push!(cirq_postcompress[d], g)
    end

    return cirq_precompress, cirq_postcompress
end


import QuantumCircuitOpt as QCOpt
include("4qubit_gates.jl")
params = double_peres()
data = QCOpt.get_data(params)

id_seq = [5, 6, 7, 8, 3, 4, 9]

cirq_precompress, cirq_postcompress = QCOpt.build_circuit_layers(id_seq, data["gates_dict"])

# QCOpt.validate_circuit(data, cirq_precompress)
println("--------------------------------")
QCOpt.validate_circuit(data, cirq_postcompress)
println("--------------------------------")