```@raw html
<align="center"/>
<img width="790px" class="display-light-only" src="assets/docs_header.png" alt="assets/docs_header.png"/>
<img width="790px" class="display-dark-only" src="assets/docs_header_dark.png" alt="assets/docs_header.png"/>
```

# Documentation

```@meta
CurrentModule = QuantumCircuitOpt
```
## Overview
**[QuantumCircuitOpt.jl](https://github.com/harshangrjn/QuantumCircuitOpt.jl)** is a Julia package which implements discrete optimization-based methods for optimal synthesis of the architecture for Quantum circuits. While programming Quantum Computers, a primary goal is to build useful and less-noisy quantum circuits from the basic building blocks, also termed as elementary gates which arise due to hardware constraints. Thus, given a desired target quantum gate and a set of elemental one- and two-qubit gates, this package provides a _provably optimal, exact_ (or, if preferable, an approximate) decomposition with minimum number of elemental gates and CNOT gates. _Note that QuantumCircuitOpt currently supports only decompositions up to five qubit gates_.

## Installation 
To use QuantumCircuitOpt, first [download and install](https://julialang.org/downloads/) Julia. Note that the current version of QuantumCircuitOpt is compatible with Julia 1.0 and later. 

The latest stable release of QuantumCircuitOpt can be installed by entering the following in the Julia REPL-mode:

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

## Acknowledgement
This work was supported by Los Alamos National Laboratory's LDRD Early Career Research Award, *"20190590ECR: Discrete Optimization Algorithms for Provable Optimal Quantum Circuit Design"*. The primary developer of this package is [Harsha Nagarajan](http://harshanagarajan.com) ([@harshangrjn](https://github.com/harshangrjn)). 

## Citing QuantumCircuitOpt
If you find QuantumCircuitOpt useful in your work, we request you to cite the following paper: 
```bibtex
@inproceedings{NagarajanHijaziCoffrin2021,
  title={Mixed-Integer Programming for Optimal Quantum Circuit Design},
  author={Nagarajan, Harsha and Hijazi, Hassan and Coffrin, Carleton},
  booktitle={Under preparation},
  year={2021}
}
```