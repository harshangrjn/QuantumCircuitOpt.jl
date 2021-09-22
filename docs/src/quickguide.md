# Quick Start Guide

## Getting started

To get started, install [QuantumCircuitOpt](https://github.com/harshangrjn/QuantumCircuitOpt.jl) and [JuMP](https://github.com/jump-dev/JuMP.jl), a modeling language layer for optimization. QuantumCircuitOpt also needs a MIP solver such as [CPLEX](https://github.com/jump-dev/CPLEX.jl) or [Gurobi](https://github.com/jump-dev/Gurobi.jl). If you prefer an open-source MIP solver, install [CBC](https://github.com/jump-dev/Cbc.jl) or [GLPK](https://github.com/jump-dev/GLPK.jl) from the Julia package manager, though be warned that the run times of QuantumCircuitOpt can be substantially slower using these open-source MIP solvers. 

# User inputs
| Necessary Inputs  | Description |
| -----------: | :----------- |
| `num_qubits`      | Number of qubits of the circuit (≥ 2).  |
| `depth`   | Maximum allowable depth for decomposition of the circuit (≥ 2)   |
| `elementary_gates` | Vector of all one and two qubit elementary gates. The menagerie of quantum gates currently supported in QuantumCircuitOpt can be found in [gates.jl](https://github.com/harshangrjn/QuantumCircuitOpt.jl/blob/master/src/gates.jl). |
| `target_gate` | Target gate which you wish to decompose using the above-mentioned `elementary_gates`.|
| `RX_discretization` | Vector of discretization angles (in radians) for `RXGate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `RY_discretization` | Vector of discretization angles (in radians) for `RYGate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `RZ_discretization` | Vector of discretization angles (in radians) for `RZGate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `Phase_discretization` | Vector of discretization angles (in radians) for `PhaseGate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `U3_θ_discretization` | Vector of discretization angles (in radians) for θ parameter in `U3Gate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `U3_ϕ_discretization` | Vector of discretization angles (in radians) for ϕ parameter in `U3Gate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `U3_λ_discretization` | Vector of discretization angles (in radians) for λ parameter in `U3Gate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `CRX_discretization` | Vector of discretization angles (in radians) for `CRXGate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `CRY_discretization` | Vector of discretization angles (in radians) for `CRYGate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `CRZ_discretization` | Vector of discretization angles (in radians) for `CRZGate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `CU3_θ_discretization` | Vector of discretization angles (in radians) for θ parameter in `CU3Gate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `CU3_ϕ_discretization` | Vector of discretization angles (in radians) for ϕ parameter in `CU3Gate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `CU3_λ_discretization` | Vector of discretization angles (in radians) for λ parameter in `CU3Gate`, if this gate is part of the above-mentioned `elementary_gates`.|
| `objective` | Choose one of the following: (a) `"minimize_depth"`, which minimizes the total depth of decomposition. For this option, include `"Identity"` matrix in the above-mentioned `elementary_gates`, (b) `"minimize_cnot"`, which minimizes the number of CNOT gates in the decomposition. |
| `decomposition_type` | Choose one of the following: (a) `"exact"`, which finds an exact decomposition if it exists, (b) `"approximate"`, which finds an approximate decomposition if an exact one does not exist; otherwise it will return an exact solution. |
| `optimizer` | Mixed-integer programming (MIP) optimizer. For various MIP solver options, check [solver.jl](https://github.com/harshangrjn/QuantumCircuitOpt.jl/blob/master/examples/solver.jl). |

| Optional Inputs  | Description |
| -----------: | :----------- |
| `initial_gate` | Intitial-condition gate to the decomposition (gate at 0th depth) (default: `"Identity"`).  | 
| `identify_real_gates` | This option identifies if all the elementary and target gates have only real entries and formulates a compact MIP formulation accordingly (default: `false`).  | 
| `input_circuit` | Input circuit representing an ensemble of elementary gates which decomposes the given target gate. This input circuit, which serves as a warm-start, can accelerate the MIP solver's search for the incumbent solution. (default: empty circuit).  | 
| `relax_integrality` | This option transforms integer variables into continuous variables (default: `false`).  |
| `optimizer_presolve` | This option enables or disables the presolve option in the chosen `optimizer` (default: `true`). Turning it off can lead to slower run times.|
| `optimizer_log` | This option enables or disables console logging for the `optimizer` (default: `true`).|
| `slack_penalty` | This option allows the user to set the penalty for minimizing the slack term in the objective, when `decomposition_type` is set to `"approximate"` (default: `1E3`).  |
| `time_limit` | This option allows the user to set the maximum time limit for the optimizer in seconds (default: `10,800`).  |



# Sample circuit synthesis
Using some of the above-described user input options, here is a sample optimization model to minimize the total depth of the decomposition for a 2-qubit controlled-Z gate. With entangling CNOT gate and the universal rotation gate with three discretized Euler angles, (θ,ϕ,λ), here is the sample code:

```julia
import QuantumCircuitOpt as QCO
using JuMP
using CPLEX

# Target: CZGate
function target_gate()
    return Array{Complex{Float64},2}([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1]) 
end

params = Dict{String, Any}(
"num_qubits" => 2, 
"depth" => 4,    
"elementary_gates" => ["U3_2", "CNot_1_2", "Identity"], 
"target_gate" => target_gate(),
       
"U3_θ_discretization" => [-π/2, 0, π/2],
"U3_ϕ_discretization" => [-π/2, 0, π/2],
"U3_λ_discretization" => [-π/2, 0, π/2],

"objective" => "minimize_depth", 
"decomposition_type" => "exact",
"optimizer" => "cplex"
)

qcm_optimizer = JuMP.optimizer_with_attributes(CPLEX.Optimizer) 
QCO.run_QCModel(params, qcm_optimizer)
```
If you prefer to decompose a target gate of your choice, update the `target_gate()` function and the 
set of `elementary_gates` accordingly in the above sample code. For more such 2-qubit and 3-qubit gate decompositions, with and without the universal unitary in the elementary gates, refer to "[examples](https://github.com/harshangrjn/QuantumCircuitOpt.jl/tree/master/examples)" folder. 

!!! warning
    Note that [QuantumCircuitOpt.jl](https://github.com/harshangrjn/QuantumCircuitOpt.jl) tries to find the global minima of a specified objective function for a given set of input one- and two-qubit gates, target gate and the total depth of the decomposition. This combinatiorial optimization problem is known to be NP-hard to compute. Hence, unlike local optimization methods, such as machine learning, in the literature, the run times for larger number of qubits and depths can be prohibitively slow.

# Extracting results
The run commands (for example, `run_QCModel`) in QuantumCircuitOpt return detailed results in the form of a dictionary. This dictionary can be saved for further processing as follows,

```julia
results = QCO.run_QCModel(params, qcm_optimizer)
```
For example, for decomposing the above controlled-Z gate, the QuantumCircuitOpt's runtime and the optimal objective value (minimum depth) can be accessed using,
```julia
results["solve_time"]
results["objective"]
```
Also, `results["solution"]` contains detailed information about the solution produced by the optimization model, which can be utilized for further analysis. 

# Visualizing results
QuantumCircuitOpt currently supports the visualization of optimal circuit decompositions obtained from the results dictionary (from above), which can be executed using,
```julia
data = QCO.get_data(params)
QCO.visualize_solution(results, data)
```
For example, for the above controlled-Z gate decomposition, the processed output of QuantumCircuitOpt is as follows: 
```
=============================================================================
Quantum Circuit Model Data

  Number of qubits: 2
  Total number of elementary gates (including discretization): 19
  Maximum depth of decomposition: 4
  Input elementary gates: ["U3_1", "U3_2", "CNot_1_2", "Identity"]
    U3 - θ discretization: [-90.0, 0.0, 90.0]
    U3 - ϕ discretization: [-90.0, 0.0, 90.0]
    U3 - λ discretization: [-90.0, 0.0, 90.0]
  Type of decomposition: exact

Optimal Circuit Decomposition

  U3_2(-90.0,0.0,0.0) * CNot_12 * U3_2(90.0,0.0,0.0) = Target gate
  Minimum optimal depth: 3
  Optimizer run time: 1.74 sec.
=============================================================================
```

