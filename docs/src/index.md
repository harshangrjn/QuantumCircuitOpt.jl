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
**[QuantumCircuitOpt](https://github.com/harshangrjn/QuantumCircuitOpt.jl)** is a Julia package which implements discrete optimization-based methods for provably optimal synthesis of an architecture for Quantum circuits. While programming Quantum Computers, a primary goal is to build useful and less-noisy quantum circuits from the basic building blocks, also termed as elementary gates which arise due to hardware constraints. Thus, given a desired quantum computation, as a target gate, and a set of elemental one- and two-qubit gates, this package provides a _provably optimal, exact_ (up to global phase and machine precision) or an approximate decomposition with minimum number of elemental gates and CNOT gates. Now, this package also supports multi-qubit gates in the elementary gates set, such as the [global rotation](https://harshangrjn.github.io/QuantumCircuitOpt.jl/dev/multi_qubit_gates/#GRGate) gate which is native to trapped ion quantum computers. _Note that QuantumCircuitOpt currently supports only decompositions of circuits up to ten qubits_.

## Installation 
To use QuantumCircuitOpt, first [download and install](https://julialang.org/downloads/) Julia. Note that the current version of QuantumCircuitOpt is compatible with Julia 1.0 and later. 

The latest stable release of QuantumCircuitOpt can be installed by entering the following in the Julia REPL-mode:

```julia
import Pkg
Pkg.add("QuantumCircuitOpt")
```

At least one mixed-integer programming (MIP) solver is required for running QuantumCircuitOpt. The well-known [Gurobi](https://github.com/jump-dev/Gurobi.jl) or IBM's [CPLEX](https://github.com/jump-dev/CPLEX.jl) solver is highly recommended, as it is fast, scaleable and can be used to solve on fairly large-scale circuits. However, the open-source MIP solver [HiGHS](https://github.com/jump-dev/HiGHS.jl) is also compatible with QuantumCircuitOpt. Gurobi (or any other MIP solver) can be installed via the package manager with

```julia
import Pkg
Pkg.add("Gurobi")
```

## Unit Tests
To run the tests in the package, run the following command after installing the QuantumCircuitOpt package.

```julia
import Pkg
Pkg.test("QuantumCircuitOpt")
```

## Acknowledgement
This work was supported by Los Alamos National Laboratory's LDRD Early Career Research award. The primary developer of this package is [Harsha Nagarajan](http://harshanagarajan.com) ([@harshangrjn](https://github.com/harshangrjn)). 

## Citing QuantumCircuitOpt
If you find QuantumCircuitOpt useful in your work, we request you to cite the following publication ([IEEE link](https://doi.org/10.1109/QCS54837.2021.00010), [arXiv link](https://arxiv.org/abs/2111.11674)):  
```bibtex
@inproceedings{NagarajanLockwoodCoffrin2021,
  title={{QuantumCircuitOpt}: An Open-source Framework for Provably Optimal Quantum Circuit Design},
  author={Nagarajan, Harsha and Lockwood, Owen and Coffrin, Carleton},
  booktitle={SC21: The International Conference for High Performance Computing, Networking, Storage, and Analysis},
  series={Second Workshop on Quantum Computing Software},
  pages={55--63},
  year={2021},
  doi={10.1109/QCS54837.2021.00010},
  organization={IEEE Computer Society}
}
```
Another publication which explores the potential of non-linear programming formulations in the QuantumCircuitOpt package is the following ([IEEE link](https://doi.org/10.1109/QCS56647.2022.00009): 
```bibtex
@inproceedings{HendersonNagarajanCoffrin2022,
  title={Exploring Non-linear Programming Formulations in {QuantumCircuitOpt} for Optimal Circuit Design},
  author={Henderson. R, Elena and Nagarajan, Harsha and Coffrin, Carleton},
  booktitle={SC22: The International Conference for High Performance Computing, Networking, Storage, and Analysis},
  series={Third Workshop on Quantum Computing Software},
  pages={36--42},
  year={2022},
  doi={10.1109/QCS56647.2022.00009},
  organization={IEEE Computer Society}
}
```