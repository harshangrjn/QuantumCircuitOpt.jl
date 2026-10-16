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

If you find QuantumCircuitOpt useful in your work, please cite the publication(s)
most relevant to your use of the package. The [Research & Publications](@ref)
page presents the original software paper, subsequent nonlinear-programming
research, and the detailed mixed-integer programming framework in chronological
order, with paper links, and copyable BibTeX citations.