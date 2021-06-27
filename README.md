<p align="center">
<img width="330px" src="https://github.com/harshangrjn/QuantumCircuitOpt.jl/blob/main/logo.png"/>
</p>

<!-- # QuantumCircuitOpt.jl -->
**QuantumCircuitOpt.jl** is a Julia package which implements discrete optimization-based algorithms for optimization of Quantum circuits architecture. Given a desired target quantum gate and a set of elemental one and two qubit gates, this package provides an exact (and approximate) decomposition with minimum number of elemental gates and CNOT gates.   

## Installation
QuantumCircuitOpt.jl is a registered package and can be installed by entering the following in the package manager:
```julia
using Pkg
Pkg.add("QuantumCircuitOpt")
```

## Usage
- Clone the repository.
- Open a terminal in the repo folder and run `julia --project=.`.
- Hit `]` to open the project environment and run `test` to run unit tests. If
  you see an error because of missing packages, run `resolve`.

On how to use this package, check the [quick start guide](https://harshangrjn.github.io/QuantumCircuitOpt.jl/stable/quickguide/#Sample-circuit-decomposition) and the "examples" folder for more such decompositions.

## Bug reports and Contributing
Please report any issues via the Github **[issue tracker](https://github.com/harshangrjn/QuantumCircuitOpt.jl/issues)**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc. 

QuantumCircuitOpt is being actively developed and suggestions or other forms of contributions are encouraged. 

## Acknowledgement
This work was supported by Los Alamos National Laboratory's LDRD Early Career Research Award, *"20190590ECR: Discrete Optimization Algorithms for Provable Optimal Quantum Circuit Design"*. The primary developer of this package is [Harsha Nagarajan](http://harshanagarajan.com) ([@harshangrjn](https://github.com/harshangrjn)). 

## Citation
If you find QuantumCircuitOpt.jl useful in your work, we request you to cite the following paper: 
```bibtex
@inproceedings{NagarajanHijaziCoffrin2021,
  title={Optimal Quantum Circuit Decompositions using Discrete Optimization Approach},
  author={Nagarajan, Harsha and Hijazi, Hassan and Coffrin, Carleton},
  booktitle={arXiv},
  year={2021}
}
```