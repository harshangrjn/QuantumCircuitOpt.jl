# Examples for QuantumCircuitOpt

`run_examples.jl` is the main driver for decomposing a target gate into a
sequence of user-selected elementary gates. Example definitions are grouped in
`2qubit_gates.jl`, `3qubit_gates.jl`, `4qubit_gates.jl`, and
`5qubit_gates.jl`. The driver runs `hadamard` with HiGHS by default. Edit
`decompose_gates` to select additional examples.

From the repository root, prepare the examples environment and run the driver:

```sh
julia --project=examples -e 'using Pkg; VERSION < v"1.11" && Pkg.develop(path="."); Pkg.instantiate()'
julia --project=examples examples/run_examples.jl
```

Commercial optimizer integrations are available but are not installed by
default. Add Gurobi or CPLEX to the examples environment, import it in
`run_examples.jl`, and replace `get_highs` with the corresponding constructor
from `examples/optimizers.jl`.