# QuantumCircuitOpt.jl Documentation

```@meta
CurrentModule = QuantumCircuitOpt
```
## Overview
 **QuantumCircuitOpt** is a Julia package which implements discrete optimization algorithms for optimization of Quantum circuits architecture. Given a desired target quantum gate and a set of elemental one and two qubit gates, this package provides an exact (and approximate) decomposition with minimum number of elemental gates and CNOT gates.

## Installation 
To use QuantumCircuitOpt, first [download and install](https://julialang.org/downloads/) Julia. Note that the current version of QuantumCircuitOpt is compatible with Julia 1.0 and later. 

The latest stable release of QuantumCircuitOpt can be installed using the Julia package manager with

```julia
import Pkg
Pkg.add("QuantumCircuitOpt")
```

At least one mixed-integer programming solver is required for running QuantumCircuitOpt. The well-known [CPLEX](https://github.com/jump-dev/CPLEX.jl) or the [Gurobi](https://github.com/jump-dev/Gurobi.jl) solver is highly recommended, as it is fast, scaleable and can be used to solve on fairly large-scale graphs. However, open-source solvers such as [Cbc](https://github.com/jump-dev/Cbc.jl) or [GLPK](https://github.com/jump-dev/GLPK.jl) is also compatible with QuantumCircuitOpt which can be installed via the package manager with

```julia
import Pkg
Pkg.add("Cbc")
```

## Unit Tests
To run the tests in the package, run the following command after installing the QuantumCircuitOpt package.

```julia
import Pkg
Pkg.test("QuantumCircuitOpt")
```

## Citing QuantumCircuitOpt
If you find QuantumCircuitOpt useful in your work, we request you to cite the following paper: 
```bibtex
@inproceedings{NagarajanHijaziCoffrinVuffrayMisra2021,
  title={Optimal Quantum Circuit Decompositions using Discrete Optimization Approach},
  author={Nagarajan, Harsha and Hijazi, Hassan and Coffrin, Carleton and Vuffray, Marc and Misra, Sidhant},
  booktitle={Mathematical Programming Computation},
  year={2021}
}
```