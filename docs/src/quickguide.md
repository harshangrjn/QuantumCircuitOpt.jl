# Quick Start Guide

## Getting started

To get started, install [QuantumCircuitOpt](https://github.com/harshangrjn/QuantumCircuitOpt.jl) and [JuMP](https://github.com/jump-dev/JuMP.jl), a modeling language layer for optimization. QuantumCircuitOpt also needs a MIP solver such as [CPLEX](https://github.com/jump-dev/CPLEX.jl) or [Gurobi](https://github.com/jump-dev/Gurobi.jl). If you prefer an open-source MIP solver, install [CBC](https://github.com/jump-dev/Cbc.jl) or [GLPK](https://github.com/jump-dev/GLPK.jl) from the Julia package manager, though be warned that the run times of QuantumCircuitOpt can be subastantially slower using these open-source MIP solvers. 

# User inputs
| Necessary Inputs  | Description |
| -----------: | :----------- |
| `num_qubits`      | Number of qubits of the circuit (≥ 2).  |
| `depth`   | Maximum allowable depth for decomposition of the circuit (≥ 2)   |
| `elementary_gates` | Vector of all one and two qubit elementary gates. Comprehensive list of gates currently supported in QuantumCircuitOpt can be found in [gates.jl](https://github.com/harshangrjn/QuantumCircuitOpt.jl/blob/main/src/gates.jl). |
| `target_gate` | Target gate which you wish to decompose using the above-mentioned `elementary_gates`.|
| `R_x_discretization` | Vector of discretization angles (in radians) for `RXGate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `R_y_discretization` | Vector of discretization angles (in radians) for `RYGate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `R_z_discretization` | Vector of discretization angles (in radians) for `RZGate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `U_θ_discretization` | Vector of discretization angles (in radians) for θ parameter in `U3Gate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `U_ϕ_discretization` | Vector of discretization angles (in radians) for ϕ parameter in `U3Gate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `U_λ_discretization` | Vector of discretization angles (in radians) for λ parameter in `U3Gate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `objective` | Choose one of the following: (a) `"minimize_depth"`, which minimizes the total depth of decomposition. For this option, include `"Identity"` matrix in the above-mentioned `elementary_gates`, (b) `"minimize_cnot"`, which minimizes the number of CNOT gates in the decomposition. |
| `decomposition_type` | Choose one of the following: (a) `"exact"`, which finds an exact decomposition if it exists, (b) `"approximate"`, which finds an approximate decomposition if an exact one does not exist; otherwise it will return an exact solution. |
| `optimizer` | Mixed-integer programming (MIP) optimizer. For various MIP solver options, check [solve.jl](https://github.com/harshangrjn/QuantumCircuitOpt.jl/blob/main/examples/solver.jl). |

| Optional Inputs  | Description |
| -----------: | :----------- |
| `initial_gate` |  Intitial condition gate to the decomposition (default: `"Identity"`).  | 
| `relax_integrality` | This option transforms integer variables into continuous variables (default: `false`).  |
| `optimizer_presolve` | This option enables or disables the presolve option in the chosen `optimizer` (default: `true`). Turning it off can lead to slower run times.|
| `optimizer_log` | This option enables or disables console logging for the `optimizer` (default: `true`).|
| `slack_penalty` | This option allows the user to set the penalty for minimizing the slack term in the objective, when `decomposition_type` is set to `"approximate"` (default: `1E3`).  |



# Sample circuit decomposition
Using some of the user input options as described above, an optimization model to minimize the total depth of decomposition for a 2-qubit controlled-Z gate can be executed as follows:

```julia
using QuantumCircuitOpt
using JuMP
using CPLEX

params = Dict{String, Any}(
"num_qubits" => 2, 
"depth" => 4,    
"elementary_gates" => ["U3", "cnot_12", "Identity"], 
"target_gate" => "controlled_Z",
       
"U_θ_discretization" => [-π/2, 0, π/2],
"U_ϕ_discretization" => [0, π/2],
"U_λ_discretization" => [0, π/2],

"objective" => "minimize_depth", 
"decomposition_type" => "exact",
"optimizer" => "cplex"
)

qcm_optimizer = JuMP.optimizer_with_attributes(CPLEX.Optimizer) 
QuantumCircuitOpt.run_QCModel(params, qcm_optimizer)
```

# Extracting results
The run commands (for example, `run_QCModel`) in QuantumCircuitOpt return detailed results in the form of a dictionary. This dictionary can be saved for further processing as follows,

```julia
results = QuantumCircuitOpt.run_QCModel(params, qcm_optimizer)
```
For example, for decomposing the above controlled-Z gate, the QuantumCircuitOpt's runtime and the optimal objective value (minimum depth) can be accessed using,
```julia
results["solve_time"]
results["objective"]
```
Also,  `results["solution"]` contains detailed information about the solution produced by the optimization model.

# Visualizing results
QuantumCircuitOpt also currently supports the visualization of optimal circuit decompositions obtained from the results dictionary (from above), which can be executed using,
```julia
data = QuantumCircuitOpt.get_data(params)
QuantumCircuitOpt.visualize_solution(results, data)
```
For example, for the above controlled-Z gate decomposition, the processed output of QuantumCircuitOpt is as follows: 
```
=============================================================================
Quantum Circuit Model Data

  Number of qubits: 2
  Total number of elementary gates (including discretization): 26
  Maximum depth of decomposition: 4
  Input elementary gates: ["U3", "cnot_12", "Identity"]
    U3 gate - θ discretization: [-90.0, 0.0, 90.0]
    U3 gate - ϕ discretization: [0.0, 90.0]
    U3 gate - λ discretization: [0.0, 90.0]
  Input target gate: controlled_Z
  Type of decomposition: exact

Optimal Circuit Decomposition

  U3 (2, (-90.0,0.0,0.0)) * cnot_12 * U3 (2, (90.0,0.0,0.0)) = controlled_Z
  Minimum optimal depth: 3
  Optimizer run time: 10.03 sec.
=============================================================================
```
 
