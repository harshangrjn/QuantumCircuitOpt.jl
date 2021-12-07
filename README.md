<p align="center">
<img width="790px" src="https://github.com/harshangrjn/QuantumCircuitOpt.jl/blob/master/logo.png" alt="https://github.com/harshangrjn/QuantumCircuitOpt.jl/tree/master/docs/src/assets/docs_header_dark.png"/>
</p>

Status: 
[![CI](https://github.com/harshangrjn/QuantumCircuitOpt.jl/actions/workflows/ci.yml/badge.svg)](https://github.com/harshangrjn/QuantumCircuitOpt.jl/actions/workflows/ci.yml)
[![codecov](https://codecov.io/gh/harshangrjn/QuantumCircuitOpt.jl/branch/master/graph/badge.svg?token=KGJWIV6QF4)](https://codecov.io/gh/harshangrjn/QuantumCircuitOpt.jl)

<!-- Stable version: [![Documentation](https://github.com/harshangrjn/QuantumCircuitOpt.jl/actions/workflows/documentation.yml/badge.svg)](https://harshangrjn.github.io/QuantumCircuitOpt.jl/stable/) -->
Dev version: [![Documentation](https://github.com/harshangrjn/QuantumCircuitOpt.jl/actions/workflows/documentation.yml/badge.svg)](https://harshangrjn.github.io/QuantumCircuitOpt.jl/dev/)


<!-- # QuantumCircuitOpt.jl -->
**QuantumCircuitOpt** is a Julia package which implements discrete optimization-based methods for provably optimal synthesis of the architecture for Quantum circuits. While programming Quantum Computers, a primary goal is to build useful and less-noisy quantum circuits from the basic building blocks, also termed as elementary gates which arise due to hardware constraints. Thus, given a desired quantum computation, as a target gate, and a set of elemental one- and two-qubit gates, this package provides a _provably optimal, exact_ (up to global phase and machine precision) or an approximate decomposition with minimum number of elemental gates and CNOT gates. _Note that QuantumCircuitOpt currently supports only decompositions of circuits up to ten qubits_.

## Installation
QuantumCircuitOpt is a registered package and can be installed by entering the following in the Julia REPL-mode:
```julia
import Pkg
Pkg.add("QuantumCircuitOpt")
```

## Usage
- Clone the repository.
- Open a terminal in the repo folder and run `julia --project=.`.
- Hit `]` to open the project environment and run `test` to run unit tests. If
  you see an error because of missing packages, run `resolve`.

On how to use this package, check the Documentation's [quick start guide](https://harshangrjn.github.io/QuantumCircuitOpt.jl/dev/quickguide/#Sample-circuit-decomposition) and the [examples](https://github.com/harshangrjn/QuantumCircuitOpt.jl/tree/master/examples) folder for several important circuit decompositions.

## Video Link
For more technical details about the package, check out this [video](https://www.youtube.com/watch?v=sf1HJW5Vmio) from the [2nd Quantum Computing Software Workshop](https://sc21.supercomputing.org/session/?sess=sess345), held in conjunction with the International Conference on Super Computing ([SC21](https://sc21.supercomputing.org)). 

## Sample Circuit Synthesis
Here is a sample usage of QuantumCircuitOpt to optimally decompose a 2-qubit controlled-Z gate ([CZGate](https://harshangrjn.github.io/QuantumCircuitOpt.jl/dev/2_qubit_gates/#CZGate)) using the entangling [CNOT](https://harshangrjn.github.io/QuantumCircuitOpt.jl/dev/2_qubit_gates/#CNotGate) gate and an one-qubit universal rotation gate ([U3Gate](https://harshangrjn.github.io/QuantumCircuitOpt.jl/dev/1_qubit_gates/#U3Gate)) with three discretized Euler angles (θ,ϕ,λ):

```julia
import QuantumCircuitOpt as QCOpt
using JuMP
using Gurobi

# Target: CZGate
function target_gate()
    return Array{Complex{Float64},2}([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1]) 
end

params = Dict{String, Any}(
"num_qubits" => 2, 
"maximum_depth" => 4,    
"elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "Identity"],
"target_gate" => target_gate(),
"objective" => "minimize_depth",
"decomposition_type" => "exact",
       
"U3_θ_discretization" => -π:π/2:π,
"U3_ϕ_discretization" => -π:π/2:π,
"U3_λ_discretization" => -π:π/2:π,
)

qcm_optimizer = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "presolve" => 1) 
QCOpt.run_QCModel(params, qcm_optimizer)
```
If you prefer to decompose a target gate of your choice, update the `target_gate()` function and the 
set of `elementary_gates` accordingly in the above sample code. 

## Bug reports and Contributing
Please report any issues via the Github **[issue tracker](https://github.com/harshangrjn/QuantumCircuitOpt.jl/issues)**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc. 

QuantumCircuitOpt is being actively developed and suggestions or other forms of contributions are encouraged. 

## Acknowledgement
This work was supported by Los Alamos National Laboratory's LDRD Early Career Research award. The primary developer of this package is [Harsha Nagarajan](http://harshanagarajan.com) ([@harshangrjn](https://github.com/harshangrjn)).

## Citing QuantumCircuitOpt
If you find QuantumCircuitOpt useful in your work, we request you to cite the following paper ([arXiv link](https://arxiv.org/abs/2111.11674)): 
```bibtex
@inproceedings{QCOpt_SC2021,
  title={{QuantumCircuitOpt}: An Open-source Framework for Provably Optimal Quantum Circuit Design},
  author={Nagarajan, Harsha and Lockwood, Owen and Coffrin, Carleton},
  booktitle={SC21: The International Conference for High Performance Computing, Networking, Storage, and Analysis},
  series={Workshop on Quantum Computing Software},
  pages={1-7},
  year={2021},
  organization={IEEE Computer Society}
}
```