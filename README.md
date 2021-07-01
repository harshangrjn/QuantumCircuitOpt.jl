<p align="center">
<img width="775px" src="https://github.com/harshangrjn/QuantumCircuitOpt.jl/blob/master/logo.png"/>
</p>


Status: 
[![CI](https://github.com/harshangrjn/QuantumCircuitOpt.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/harshangrjn/QuantumCircuitOpt.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/harshangrjn/QuantumCircuitOpt.jl/branch/master/graph/badge.svg?token=KGJWIV6QF4)](https://codecov.io/gh/harshangrjn/QuantumCircuitOpt.jl)
[![Documentation](https://github.com/harshangrjn/QuantumCircuitOpt.jl/actions/workflows/documentation.yml/badge.svg)](https://harshangrjn.github.io/QuantumCircuitOpt.jl/stable/)


<!-- # QuantumCircuitOpt.jl -->
**QuantumCircuitOpt.jl** is a Julia package which implements discrete optimization-based methods for optimal synthesis of the architecture for Quantum circuits. While programming Quantum Computers, a primary goal is to build useful and less-noisy quantum circuits from the basic building blocks, also termed as elementary gates which arise due to hardware constraints. Thus, given a desired target quantum gate and a set of elemental one- and two-qubit gates, this package provides a _provably optimal, exact_ (or, if preferable, an approximate) decomposition with minimum number of elemental gates and CNOT gates. _Note that QuantumCircuitOpt currently supports only decompositions up to three-qubit gates_.

## Installation
QuantumCircuitOpt.jl is a registered package and can be installed by entering the following in the Julia REPL-mode:
```julia
import Pkg
Pkg.add("QuantumCircuitOpt")
```

## Usage
- Clone the repository.
- Open a terminal in the repo folder and run `julia --project=.`.
- Hit `]` to open the project environment and run `test` to run unit tests. If
  you see an error because of missing packages, run `resolve`.

On how to use this package, check the [quick start guide](https://harshangrjn.github.io/QuantumCircuitOpt.jl/stable/quickguide/#Sample-circuit-decomposition) and the "[examples](https://github.com/harshangrjn/QuantumCircuitOpt.jl/tree/master/examples)" folder for more such decompositions.

## Sample Circuit Synthesis
Here is the a sample usage of QuantumCircuitOpt to optimally decompose a 2-qubit controlled-Z gate using the entangling CNOT gate the universal rotation gate with three discretized Euler angles, (θ,ϕ,λ) :

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
"elementary_gates" => ["U3", "cnot_12", "Identity"], 
"target_gate" => target_gate(),
       
"U_θ_discretization" => [-π/2, 0, π/2],
"U_ϕ_discretization" => [0, π/2],
"U_λ_discretization" => [0, π/4],

"objective" => "minimize_depth", 
"decomposition_type" => "exact",
"optimizer" => "cplex"
)

qcm_optimizer = JuMP.optimizer_with_attributes(CPLEX.Optimizer) 
QCO.run_QCModel(params, qcm_optimizer)
```
If you prefer to decompose a target gate of your choice, update the `target_gate()` function and the 
set of `elementary_gates` accordingly in the above sample code. 

## Bug reports and Contributing
Please report any issues via the Github **[issue tracker](https://github.com/harshangrjn/QuantumCircuitOpt.jl/issues)**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc. 

QuantumCircuitOpt is being actively developed and suggestions or other forms of contributions are encouraged. 

## Acknowledgement
This work was supported by Los Alamos National Laboratory's LDRD Early Career Research Award, *"20190590ECR: Discrete Optimization Algorithms for Provable Optimal Quantum Circuit Design"*. The primary developer of this package is [Harsha Nagarajan](http://harshanagarajan.com) ([@harshangrjn](https://github.com/harshangrjn)). 

<!-- ## Citation
If you find QuantumCircuitOpt.jl useful in your work, we request you to cite the following paper: 
```bibtex
@inproceedings{NagarajanHijaziCoffrin2021,
  title={Mixed-Integer Programming for Optimal Quantum Circuit Design},
  author={Nagarajan, Harsha and Hijazi, Hassan and Coffrin, Carleton},
  booktitle={Under preparation},
  year={2021}
}
``` -->