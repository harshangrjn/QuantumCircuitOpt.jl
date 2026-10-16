<h1 align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/src/assets/logo_header_dark.png">
    <img alt="QuantumCircuitOpt.jl" src="docs/src/assets/logo_header_light.png" width="900">
  </picture>
  <br>
  A Julia Package for Optimal Quantum Circuit Design
</h1>

[![CI](https://github.com/harshangrjn/QuantumCircuitOpt.jl/actions/workflows/ci.yml/badge.svg?branch=master)](https://github.com/harshangrjn/QuantumCircuitOpt.jl/actions/workflows/ci.yml)
[![Documentation](https://github.com/harshangrjn/QuantumCircuitOpt.jl/actions/workflows/documentation.yml/badge.svg?branch=master)](https://harshangrjn.github.io/QuantumCircuitOpt.jl/dev/)
[![Codecov](https://codecov.io/gh/harshangrjn/QuantumCircuitOpt.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/harshangrjn/QuantumCircuitOpt.jl/branch/master)
[![Latest release](https://img.shields.io/github/v/release/harshangrjn/QuantumCircuitOpt.jl)](https://github.com/harshangrjn/QuantumCircuitOpt.jl/releases/latest)

**QuantumCircuitOpt** is a Julia package which implements discrete optimization-based methods for provably optimal synthesis of an architecture for quantum circuits. While programming quantum computers, a primary goal is to build useful and less-noisy quantum circuits from the basic building blocks, also termed as elementary gates which arise due to hardware constraints. Thus, given a desired quantum computation, as a target gate, and a set of elemental one- and two-qubit gates, this package provides a _provably optimal, exact_ (up to global phase and machine precision) or an approximate decomposition with minimum number of elemental gates and CNOT gates. Now, this package also supports multi-qubit gates in the elementary gates set, such as the [global rotation](https://harshangrjn.github.io/QuantumCircuitOpt.jl/dev/multi_qubit_gates/#GRGate) gate. _Note that QuantumCircuitOpt currently supports only decompositions of circuits up to ten qubits_.

Overall, QuantumCircuitOpt can be a useful tool for researchers and developers working on quantum algorithms or quantum computing applications, as it can help to reduce the resource requirements of quantum computations, making them more practical and efficient.

## Installation

QuantumCircuitOpt requires Julia 1.10 or later. It is a registered package and
can be installed by entering the following in the Julia REPL:

```julia
import Pkg
Pkg.add("QuantumCircuitOpt")
```

## Quick start

The following example decomposes a two-qubit controlled-Z gate using CNOT and
discretized U3 gates:

```julia
import QuantumCircuitOpt as QCOpt
using Gurobi
using JuMP

params = Dict{String,Any}(
    "num_qubits" => 2,
    "maximum_depth" => 4,
    "elementary_gates" => ["U3_1", "U3_2", "CNot_1_2", "Identity"],
    "target_gate" => ComplexF64[
        1 0 0 0
        0 1 0 0
        0 0 1 0
        0 0 0 -1
    ],
    "objective" => "minimize_depth",
    "decomposition_type" => "exact_optimal",
    "U3_θ_discretization" => -π:(π / 2):π,
    "U3_ϕ_discretization" => -π:(π / 2):π,
    "U3_λ_discretization" => -π:(π / 2):π,
)

optimizer = optimizer_with_attributes(Gurobi.Optimizer, "Presolve" => 1)
results = QCOpt.run_QCModel(params, optimizer)
```

See the [quick-start guide](https://harshangrjn.github.io/QuantumCircuitOpt.jl/dev/quickguide/)
for all circuit and model options, and the
[`examples`](https://github.com/harshangrjn/QuantumCircuitOpt.jl/tree/master/examples)
directory for additional decompositions.

## Testing

To test a development checkout:

```sh
git clone https://github.com/harshangrjn/QuantumCircuitOpt.jl.git
cd QuantumCircuitOpt.jl
julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.test()'
```

## Video Links

For more technical details about the package, check out these video links:

- November 2022: Presentation [link](https://vimeo.com/771366943/01810daa4e) from the [Third Quantum Computing Software Workshop](https://sc22.supercomputing.org/session/?sess=sess423), held in conjunction with the International Conference on Super Computing ([SC22](https://sc22.supercomputing.org)).
- July 2022: Presentation [link](https://www.youtube.com/watch?v=OeONXwD4JJY) from the [JuliaCon 2022](https://pretalx.com/juliacon-2022/talk/KJTGC3/) conference.
- November 2021: Presentation [link](https://www.youtube.com/watch?v=sf1HJW5Vmio) from the [Second Quantum Computing Software Workshop](https://sc21.supercomputing.org/session/?sess=sess345), held in conjunction with the International Conference on Super Computing ([SC21](https://sc21.supercomputing.org)).

## Contributing

Bug reports, documentation corrections, and feature requests are welcome in the
[issue tracker](https://github.com/harshangrjn/QuantumCircuitOpt.jl/issues).

## Acknowledgement

This work was supported by Los Alamos National Laboratory's LDRD Early Career
Research award. The primary developer is
[Harsha Nagarajan](https://harshanagarajan.com)
([@harshangrjn](https://github.com/harshangrjn)).

## Citing QuantumCircuitOpt

If you find QuantumCircuitOpt useful in your work, we request you to cite the following paper ([IEEE link](https://doi.org/10.1109/QCS54837.2021.00010), [arXiv link](https://arxiv.org/abs/2111.11674)):

```bibtex
@inproceedings{QCOpt_SC2021,
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
@misc{nagarajan2025provably,
  title={Provably Optimal Quantum Circuits with Mixed-Integer Programming},
  author={Nagarajan, Harsha and Szabó, Zsolt},
  year={2025},
  eprint={2510.00649},
  archivePrefix={arXiv},
  url={https://arxiv.org/abs/2510.00649},
}
```

Another publication which explores the potential of non-linear programming formulations in QuantumCircuitOpt is the following: [IEEE link](https://doi.org/10.1109/QCS56647.2022.00009), [arXiv link](https://arxiv.org/abs/2310.18281).