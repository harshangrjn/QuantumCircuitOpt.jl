```@raw html
<img width="790" class="display-light-only" src="assets/logo_header_light.png" alt="QuantumCircuitOpt.jl"/>
<img width="790" class="display-dark-only" src="assets/logo_header_dark.png" alt="QuantumCircuitOpt.jl"/>
```

# Documentation

```@meta
CurrentModule = QuantumCircuitOpt
```

## Overview

**[QuantumCircuitOpt](https://github.com/harshangrjn/QuantumCircuitOpt.jl)** is
a Julia package that implements discrete optimization methods for provably
optimal quantum circuit synthesis. Quantum processors expose restricted,
architecture-dependent native gate sets, so a desired unitary operation must
be compiled into a sequence of elementary gates that is both functionally
correct and economical in scarce quantum resources.

Given a target unitary matrix and a set of available elementary gates,
QuantumCircuitOpt formulates the synthesis problem as a mathematical
optimization model. It can compute provably optimal exact decompositions (up
to machine precision), optimal decompositions up to a global phase, or
approximate decompositions. Supported objectives include minimizing circuit
depth, CNOT count, and T-gate count.

The elementary gate library supports one- and two-qubit gates as well as
multi-qubit operations such as the [global rotation gate](@ref GRGate), which
is native to trapped-ion architectures. QuantumCircuitOpt currently supports
circuit decompositions involving up to ten qubits. By providing optimality
certificates, it is useful for benchmarking heuristic quantum compilers,
studying hardware-aware gate sets, and identifying minimum-resource
implementations of small but important quantum operations.

## Installation

QuantumCircuitOpt requires Julia 1.10 or later. Install the latest registered
release from the Julia REPL:

```julia
import Pkg
Pkg.add("QuantumCircuitOpt")
```

At least one JuMP-compatible mixed-integer programming solver is also required.
[Gurobi](https://github.com/jump-dev/Gurobi.jl) is recommended for difficult
circuit-synthesis instances; IBM
[CPLEX](https://github.com/jump-dev/CPLEX.jl) is also supported. The
open-source [HiGHS](https://github.com/jump-dev/HiGHS.jl) solver is compatible
with QuantumCircuitOpt, although commercial solvers may be faster on large
instances. Install Gurobi from the Julia package manager with:

```julia
import Pkg
Pkg.add("Gurobi")
```

## Unit Tests

Run the registered package's tests with:

```julia
import Pkg
Pkg.test("QuantumCircuitOpt")
```

## Acknowledgement

This work was supported by Los Alamos National Laboratory's LDRD Early Career
Research award. The primary developer is
[Harsha Nagarajan](https://harshanagarajan.com)
([@harshangrjn](https://github.com/harshangrjn)).

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
A detailed version of the mixed-integer programming framework is
available in this [paper](https://arxiv.org/abs/2510.00649):

```bibtex
@misc{nagarajan2025provablyoptimalquantumcircuits,
  title={Provably Optimal Quantum Circuits with Mixed-Integer Programming},
  author={Harsha Nagarajan and Zsolt Szabó},
  year={2025},
  eprint={2510.00649},
  archivePrefix={arXiv},
  primaryClass={quant-ph},
  url={https://arxiv.org/abs/2510.00649},
}
```

Another publication which explores the potential of non-linear programming formulations in the QuantumCircuitOpt package is the following ([IEEE link](https://doi.org/10.1109/QCS56647.2022.00009), [arXiv link](https://arxiv.org/abs/2310.18281)): 
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