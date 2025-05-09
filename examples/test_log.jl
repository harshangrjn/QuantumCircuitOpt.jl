###############################################################
#  Minimal data model + gate parser           (UNCHANGED)
###############################################################
using Printf                       # <— BitSets is required

struct Gate
    label  :: String
    qubits :: BitSet
end

function parse_gate(s::AbstractString)::Gate
    head, params = occursin('(', s) ? split(s, '(', limit = 2) : (s, "")
    parts = split(head, '_')
    name  = parts[1] * (params == "" ? "" : "(" * params)
    qs    = BitSet(parse.(Int, parts[2:end]))
    Gate(name, qs)
end

###############################################################
#  NEW: explode "…x…" kron-notation into individual Gate items
###############################################################
"""
    explode_kron(seq::Vector{String}) -> Vector{Gate}

Splits any item containing `x` into its tensor-product factors,
then parses each factor with `parse_gate`.
"""
function explode_kron(seq::Vector{String})::Vector{Gate}
    out = Gate[]
    for item in seq
        for term in split(item, 'x')
            push!(out, parse_gate(term))
        end
    end
    return out
end


###############################################################
#  Greedy ASAP / slide-left scheduler        (UNCHANGED)
###############################################################
function compress(gates::Vector{Gate}, n::Int)
    layers = Vector{Vector{Gate}}()
    masks  = Vector{BitSet}()
    full   = BitSet(1:n)

    for g in gates
        @assert g.qubits ⊆ full "gate $(g.label) touches an invalid qubit (n = $n)"

        pos = length(layers) + 1
        for l = length(layers):-1:1
            if !isempty(intersect(g.qubits, masks[l]))
                break
            end
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


###############################################################
#  Accurate Kronecker-style printer          (UNCHANGED)
###############################################################
function kron_layer(layer::Vector{Gate}, n::Int)::String
    sorted = sort(layer; by = g -> minimum(g.qubits))
    facs   = String[]; qptr = 1
    for g in sorted
        firstq, lastq = minimum(g.qubits), maximum(g.qubits)
        for _ in qptr:firstq-1
            push!(facs, "I")
        end
        if length(g.qubits) == 1
            push!(facs, g.label)
        else
            push!(facs, @sprintf("%s_{%s}", g.label,
                                 join(sort(collect(g.qubits)), ",")))
        end
        qptr = lastq + 1
    end
    for _ in qptr:n
        push!(facs, "I")
    end
    return join(facs, " ⊗ ")
end


###############################################################
#  Demo that includes the new “x” syntax
###############################################################
function demo()
    println("=== Example A (original list, unchanged) ===")
    gatesA_raw = [
        "X_1", "X_2", "H_1",
        "CNot_1_2", "H_1",
        "CNot_2_3", "H_3" ]
    layersA = compress(explode_kron(gatesA_raw), 3)
    for (d, lay) in enumerate(layersA)
        println(@sprintf("depth %d : %s", d, kron_layer(lay, 3)))
    end

    println("\n=== Example B (with kron syntax) ===")
    gatesB_raw = [
        "X_1xX_2",                       # kron(X1, X2)
        "H_3",
        "CNot_2_3xZ_1",                  # kron(CNot{2,3}, Z1)
        "U3_2(45,0,0)" ]
    layersB = compress(explode_kron(gatesB_raw), 3)

    for (d, lay) in enumerate(layersB)
        println(@sprintf("depth %d : %s", d, kron_layer(lay, 3)))
    end

    gates_sol = [
    "H_1",
    "CNot_1_2",
    "H_1",
    "X_1xX_2"
    ]

    layers = compress(explode_kron(gates_sol), 2)   # 2-qubit register
    @show layers

    println("\n=== Example C (with kron syntax) ===")
    for (d, lay) in enumerate(layers)
        println(@sprintf("depth %d : %s", d, kron_layer(lay, 2)))
    end

end

demo()
