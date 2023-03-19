# Quick Start Guide


## Framework
Building on the recent success of [Julia](https://julialang.org), [JuMP](https://github.com/jump-dev/JuMP.jl) and mixed-integer programming (MIP) solvers, [QuantumCircuitOpt](https://github.com/harshangrjn/QuantumCircuitOpt.jl) (or QCOpt), is an open-source toolkit for optimal quantum circuit design. As illustrated in the figure below, QCOpt is written in Julia, a relatively new and fast dynamic programming language used for technical computing with support for extensible type system and meta-programming. At a high level, QCOpt provides an abstraction layer to achieve two primary goals:
1. To capture user-specified inputs, such as a desired quantum computation and the available hardware gates, and build a JuMP model of an MIP formulation, and 
2. To extract, analyze and post-process the solution from the JuMP model to provide exact and approximate circuit decompositions, up to a global phase and machine precision.

```@raw html
<align="center"/>
<img width="550px" class="display-light-only" src="../assets/QCOpt_framework.png" alt="../assets/QCOpt_framework.png"/>
<img width="550px" class="display-dark-only" src="../assets/QCOpt_framework_dark.png" alt="../assets/QCOpt_framework.png"/>
```

## Video tutorials
- November 2022: Presentation [link](https://vimeo.com/771366943/01810daa4e) from the [Third Quantum Computing Software Workshop](https://sc22.supercomputing.org/session/?sess=sess423), held in conjunction with the International Conference on Super Computing ([SC22](https://sc22.supercomputing.org)). This video will introduce the importance of nonlinear programming formulations for optimal quantum circuit design. More technical details can be found in this [paper](https://doi.org/10.1109/QCS56647.2022.00009). 

- July 2022: Presentation at the [JuliaCon 2022](https://pretalx.com/juliacon-2022/talk/KJTGC3/) introduces the package in greater depth and how to use it's various features. 

```@raw html 
<align="center"/>
<a href="https://www.youtube.com/watch?v=OeONXwD4JJY">
    <img alt="Youtube-link" src="../assets/video_img_2.png"
    width=500" height="350">
</a>
```

- November 2021: Presentation from the [2nd Quantum Computing Software Workshop](https://sc21.supercomputing.org/session/?sess=sess345), held in conjunction with the International Conference on Super Computing ([SC21](https://sc21.supercomputing.org)), will introduce the technicalities underlying the package, which can also be found in this [paper](https://doi.org/10.1109/QCS54837.2021.00010).

```@raw html 
<align="center"/>
<a href="https://www.youtube.com/watch?v=sf1HJW5Vmio">
    <img alt="Youtube-link" src="../assets/video_img_1.png"
    width=500" height="350">
</a>
```

## Getting started
To get started, install [QCOpt](https://github.com/harshangrjn/QuantumCircuitOpt.jl) and [JuMP](https://github.com/jump-dev/JuMP.jl), a modeling language layer for optimization. QCOpt also needs a MIP solver such as [Gurobi](https://github.com/jump-dev/Gurobi.jl) or IBM's [CPLEX](https://github.com/jump-dev/CPLEX.jl). If you prefer an open-source MIP solver, install [HiGHS](https://github.com/jump-dev/HiGHS.jl) from the Julia package manager, though be warned that the run times of QCOpt can be substantially slower using any of the open-source MIP solvers. 

# User inputs
QCOpt takes two types of user-defined input specifications. The first type of input contains all the necessary circuit specifications. This is given by a dictionary in Julia, which is a collection of key-value pairs, where every key is of the type `String`, which admits values of various types. Below is the list of allowable keys for the dictionary, given in column 1, and it's respective values with descriptions, given in column 2. This input dictionary is represented as `params` in all the [example](https://github.com/harshangrjn/QuantumCircuitOpt.jl/tree/master/examples) circuit decompositions. 

| Mandatory circuit specifications  | Description |
| -----------: | :----------- |
| `num_qubits`      | Number of qubits of the circuit (≥ 2).  |
| `maximum_depth`   | Maximum allowable depth for decomposition of the circuit (≥ 2).   |
| `elementary_gates` | Vector of all one and two qubit elementary gates. The menagerie of quantum gates currently supported in QCOpt can be found in [gates.jl](https://github.com/harshangrjn/QuantumCircuitOpt.jl/blob/master/src/gates.jl). |
| `target_gate` | Target unitary gate which you wish to decompose using the above-mentioned `elementary_gates`.|
| `objective` | Choose one of the following: (a) `minimize_depth`, which minimizes the total number of one- and two-qubit gates. For this option, include `Identity` matrix in the above-mentioned `elementary_gates`, (b) `minimize_cnot`, which minimizes the number of CNOT gates in the decomposition. (default: `minimize_depth`) |
| `decomposition_type` | Choose one of the following: (a) `exact_optimal`, which finds an exact, provably optimal, decomposition if it exists, (b) `exact_feasible`, which finds any feasible exact decomposition, but not necessarily an optimal one if it exists, (c) `optimal_global_phase`, which finds an optimal (best) circuit decomposition if it exists, up to a global phase (d) `approximate`, which finds an approximate decomposition if an exact one does not exist; otherwise it will return an exact decomposition (default: `exact_optimal`)|

If the above-specified `elementary_gates` contain gates with continuous angle parameters, then the following mandarotry input angle discretizations have to be specified in addition to the above inputs: 

| Mandatory angle discretizations  | Description |
| -----------: | :----------- |
| `RX_discretization` | Vector of discretization angles (in radians) for `RXGate`. Input this only if this gate is part of the above-mentioned `elementary_gates`.|
| `RY_discretization` | Vector of discretization angles (in radians) for `RYGate`. Input this only if this gate is part of the above-mentioned `elementary_gates`.|
| `RZ_discretization` | Vector of discretization angles (in radians) for `RZGate`. Input this only if this gate is part of the above-mentioned `elementary_gates`.|
| `Phase_discretization` | Vector of discretization angles (in radians) for `PhaseGate`. Input this only if this gate is part of the above-mentioned `elementary_gates`.|
| `U3_θ_discretization` | Vector of discretization angles (in radians) for θ parameter in `U3Gate`. Input this only if this gate is part of the above-mentioned `elementary_gates`.|
| `U3_ϕ_discretization` | Vector of discretization angles (in radians) for ϕ parameter in `U3Gate`. Input this only if this gate is part of the above-mentioned `elementary_gates`.|
| `U3_λ_discretization` | Vector of discretization angles (in radians) for λ parameter in `U3Gate`. Input this only if this gate is part of the above-mentioned `elementary_gates`.|
| `CRX_discretization` | Vector of discretization angles (in radians) for `CRXGate`. Input this only if this gate is part of the above-mentioned `elementary_gates`.|
| `CRY_discretization` | Vector of discretization angles (in radians) for `CRYGate`. Input this only if this gate is part of the above-mentioned `elementary_gates`.|
| `CRZ_discretization` | Vector of discretization angles (in radians) for `CRZGate`. Input this only if this gate is part of the above-mentioned `elementary_gates`.|
| `CU3_θ_discretization` | Vector of discretization angles (in radians) for θ parameter in `CU3Gate`. Input this only if this gate is part of the above-mentioned `elementary_gates`.|
| `CU3_ϕ_discretization` | Vector of discretization angles (in radians) for ϕ parameter in `CU3Gate`. Input this only if this gate is part of the above-mentioned `elementary_gates`.|
| `CU3_λ_discretization` | Vector of discretization angles (in radians) for λ parameter in `CU3Gate`. Input this only if this gate is part of the above-mentioned `elementary_gates`.|

In addition, here is a list of *optional* circuit specifications, which can be added to the above set of inputs, to accelerate the performance of the QCOpt package:

| Optional circuit specifications  | Description |
| -----------: | :----------- |
| `initial_gate` | Intitial-condition gate to the decomposition (gate at 0th depth) (default: `Identity`).  | 
| `set_cnot_lower_bound` | This option sets a lower bound on the total number of CNot or CX gates which an optimal decomposition can admit.  |
| `set_cnot_upper_bound` | This option sets an upper bound on the total number of CNot/CX gates which an optimal decomposition can admit. Note that both `set_cnot_lower_bound` and `set_cnot_upper_bound` can also be set to an identitcal value to fix the number of CNot/CX gates in the optimal decomposition.|
| `input_circuit` | Input circuit representing an ensemble of elementary gates which decomposes the given target gate. This input circuit, which serves as a warm-start, can accelerate the MIP solver's search for the incumbent solution. (default: empty circuit).  | 


# Optimization model inputs
The second set of inputs for QCOpt contains all the optional specifications for the underlying optimization models. This is given by a dictionary in Julia, which is a collection of key-value pairs, where every key is of the type `Symbol`, which admits values of various types. Below is the list of allowable keys for this dictionary, given in column 1, and it's respective values with descriptions, given in column 2. This input dictionary is an optional one, as it's default values are already set in [`types.jl`](https://github.com/harshangrjn/QuantumCircuitOpt.jl/blob/master/src/types.jl) correspnding to an optimal performance of the QCOpt package. Further, this dictionary is an optional argument while executing functions, `build_QCModel` and `run_QCModel` only.

| Optional model inputs  | Description |
| -----------: | :----------- |
|`model_type`| The type of implemented MIP model to optimize in QCOpt (default: `compact_formulation`). | 
|`commute_gate_constraints`| This option activates the valid constraints to eliminate pairs of commuting gates in the elementary (native) gates set (default: `true`)| 
|`involutory_gate_constraints`| This option activates the valid constraints to eliminate pairs of involutory gates in the elementary (native) gates set (default: `true`)| 
|`redundant_gate_pair_constraints`| This option activates the valid constraints to eliminate redundant pairs of gates in the elementary (native) gates set (default: `true`)| 
|`identity_gate_symmetry_constraints`| This option activates the valid constraints to eliminate symmetry in the Identity gate in the decomposition (default: `true`)| 
|`idempotent_gate_constraints`| This option activates the valid constraints to eliminate idempotent gates in the elementary (native) gates set (default: `false`)| 
|`convex_hull_gate_constraints`| This option activates the valid constraints to apply convex hull of complex entries in the elementary (native) gates set (default: `false`)| 
|`fix_unitary_variables`| This option evaluates all the fixed-valued indices of unitary matrix varaibles (`U_var`) at every depth, and appropriately builds the optimization model (default: `true`)|
|`visualize_solution`| This option activates the visualization of the optimal circuit decomposition (default: `true`)| 
| `relax_integrality` | This option transforms integer variables into continuous variables (default: `false`).  |
| `optimizer_log` | This option enables or disables console logging for the `optimizer` (default: `true`).|
| `objective_slack_penalty` | This option sets the penalty for minimizing the slack term in the objective, when `decomposition_type` is set to `approximate` (default: `1E3`).  |
| `time_limit` | This option allows sets the maximum time limit for the optimizer in seconds (default: `10,800`).  |

# Sample circuit synthesis
Using the above-described mandatory and optional user inputs, here is a sample circuit decomposition to minimize the total depth for implementing a 2-qubit controlled-Z gate ([CZGate](https://harshangrjn.github.io/QuantumCircuitOpt.jl/dev/2_qubit_gates/#CZGate)), with entangling [CNOT](https://harshangrjn.github.io/QuantumCircuitOpt.jl/dev/2_qubit_gates/#CNotGate) gate and the one-qubit, universal rotation gate ([U3Gate](https://harshangrjn.github.io/QuantumCircuitOpt.jl/dev/1_qubit_gates/#U3Gate)) with three discretized Euler angles (θ,ϕ,λ):

```julia
import QuantumCircuitOpt as QCOpt
using JuMP
using Gurobi

# Target: CZGate
function target_gate()
    return Array{Complex{Float64},2}([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1]) 
end

# Circuit specifications (mandatory)
params = Dict{String, Any}(
"num_qubits" => 2, 
"maximum_depth" => 4,    
"elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "Identity"],
"target_gate" => target_gate(),
"objective" => "minimize_depth",
"decomposition_type" => "exact_optimal",
       
"U3_θ_discretization" => -π:π/2:π,
"U3_ϕ_discretization" => -π:π/2:π,
"U3_λ_discretization" => -π:π/2:π,
)

# Optimization model inputs (optional)
model_options = Dict{Symbol, Any}(:model_type => "compact_formulation",
                                  :visualize_solution => true)

qcm_optimizer = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "presolve" => 1)
QCOpt.run_QCModel(params, qcm_optimizer; options = model_options)
```

If you prefer to decompose a target gate of your choice, update the `target_gate()` function and the 
set of `elementary_gates` accordingly in the above sample code. For more such circuit decompositions, with various types of elementary gates, refer to [examples](https://github.com/harshangrjn/QuantumCircuitOpt.jl/tree/master/examples) folder. 

!!! warning
    Note that [QCOpt](https://github.com/harshangrjn/QuantumCircuitOpt.jl) tries to find the global minima of a specified objective function for a given set of input one- and two-qubit gates, target gate and the total depth of the decomposition. This combinatiorial optimization problem is known to be NP-hard to compute in the size of `num_qubits`, `maximum_depth` and `elementary_gates`. 

!!! tip
    Run times of [QCOpt](https://github.com/harshangrjn/QuantumCircuitOpt.jl)'s mathematical optimization models are significantly faster using [Gurobi](https://www.gurobi.com) as the underlying mixed-integer programming (MIP) solver. Note that this solver's individual-usage license is available [free](https://www.gurobi.com/academia/academic-program-and-licenses/) for academic purposes. 

# Extracting results
The run commands (for example, `run_QCModel`) in QCOpt return detailed results in the form of a dictionary. This dictionary can be saved for further processing as follows,

```julia
results = QCOpt.run_QCModel(params, qcm_optimizer)
```
For example, for decomposing the above controlled-Z gate, the QCOpt's runtime and the optimal objective value (minimum depth) can be accessed using,
```julia
results["solve_time"]
results["objective"]
```
Also, `results["solution"]` contains detailed information about the solution produced by the optimization model, which can be utilized for further analysis. 

# Visualizing results
QCOpt currently supports the visualization of optimal circuit decompositions obtained from the results dictionary (from above), which can be executed using,
```julia
data = QCOpt.get_data(params)
QCOpt.visualize_solution(results, data)
```
For example, for the above controlled-Z gate decomposition, the processed output of QCOpt is as follows: 
```
=============================================================================
QuantumCircuitOpt version: v0.5.1

Quantum Circuit Model Data
  Number of qubits: 2
  Total number of elementary gates (after presolve): 72
  Maximum depth of decomposition: 4
  Elementary gates: ["U3_1", "U3_2", "CNot_1_2", "Identity"]
    U3_θ discretization: [-180.0, -90.0, 0.0, 90.0, 180.0]
    U3_ϕ discretization: [-180.0, -90.0, 0.0, 90.0, 180.0]
    U3_λ discretization: [-180.0, -90.0, 0.0, 90.0, 180.0]
  Type of decomposition: exact_optimal
  MIP optimizer: Gurobi

Optimal Circuit Decomposition
  U3_2(-90.0,0.0,0.0) * CNot_1_2 * U3_2(90.0,0.0,0.0) = Target gate
  Minimum optimal depth: 3
  Optimizer run time: 3.01 sec.
=============================================================================
```

