###############################################################################
# depth-compression without an external graph package
###############################################################################

struct Gate
    name::Symbol
    qubits::Vector{Int}
end

Base.show(io::IO, g::Gate) = print(io, "$(g.name)$(g.qubits)")

# ————————————————————————————————————————————————————————————————
# 1.  Build adjacency lists  (vector of vectors)
# ————————————————————————————————————————————————————————————————
"""
    adjacency(gates) -> Vector{Vector{Int}}

`adj[v]` contains the indices of gates that **share a qubit** with gate `v`.
"""
function adjacency(gates::Vector{Gate})
    n      = length(gates)
    adj    = [Int[] for _ in 1:n]
    qubits = [Set(g.qubits) for g in gates]           # cache for speed

    for v in 1:n-1, w in v+1:n
        if !isempty(intersect(qubits[v], qubits[w]))
            push!(adj[v], w)
            push!(adj[w], v)
        end
    end
    return adj
end

# ————————————————————————————————————————————————————————————————
# 2.  Greedy DSATUR-style colouring using only the adjacency list
# ————————————————————————————————————————————————————————————————
"""
    greedy_depths(gates) -> Vector{Int}

Return a depth label for each gate (1-based colours).
"""
function greedy_depths(gates::Vector{Gate})
    adj   = adjacency(gates)
    n     = length(adj)
    order = sortperm(length.(adj); rev = true)        # high degree first
    depth = zeros(Int, n)

    for v in order
        forbidden = depth[adj[v]]
        d = 1
        while d in forbidden
            d += 1
        end
        depth[v] = d
    end
    return depth
end

# ————————————————————————————————————————————————————————————————
# 3.  Group gates by depth
# ————————————————————————————————————————————————————————————————
"""
    compress_circuit(gates) -> layers::Vector{Vector{Gate}}
"""
function compress_circuit(gates::Vector{Gate})
    depths = greedy_depths(gates)
    maxd   = maximum(depths)
    layers = [Gate[] for _ in 1:maxd]
    for (g,d) in zip(gates, depths)
        push!(layers[d], g)
    end
    return layers
end

# ————————————————————————————————————————————————————————————————
# 4.  Demo run
# ————————————————————————————————————————————————————————————————
circ = Gate[
    Gate(:CVdg_24, [2,4]),
    Gate(:CVdg_34, [3,4]),
    Gate(:CNOT_12, [1,2]),
    Gate(:CNOT_23, [2,3]),
    Gate(:CV_34,   [3,4]),
    Gate(:CVdg_14, [1,4])
]

layers = compress_circuit(circ)
println("Compressed depth = $(length(layers)) layers")
for (d, layer) in pairs(layers)
    println("  depth $d : ", join(layer, "  ⊗  "))
end
