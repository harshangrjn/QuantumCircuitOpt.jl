<p align="center">
<img width="400px" src="https://github.com/harshangrjn/QuantumCircuitOpt.jl/blob/main/logo.png"/>
</p>

# QuantumCircuitOpt.jl
This is a Julia package which implements discrete optimization algorithms for optimization of Quantum circuits architecture. Given a desired target quantum gate and a set of elemental one and two qubit gates, this package provides an exact (and approximate) decomposition with minimum number of elemental gates and  CNOT gates.   

## Usage
- Clone the repository.
- Open a terminal in the repo folder and run `julia --project=.`.
- Hit `]` to open the project environment and run `test` to run unit tests. If
  you see an error because of missing packages, run `resolve`.

Check the "examples" folder on how to use this package.

## Bug reports and support
Please report any issues via the Github **[issue tracker]**. All types of issues are welcome and encouraged; this includes bug reports, documentation typos, feature requests, etc.

## Acknowledgement
This work was supported by Los Alamos National Laboratory's LDRD Early Career Research Award, *"20190590ECR: Discrete Optimization Algorithms for Provable Optimal Quantum Circuit Design"*. The primary developer of this package is [Harsha Nagarajan](http://harshanagarajan.com) (@harshangrjn). 

## Citation
If you find QuantumCircuitOpt.jl useful in your work, we request you to cite the following paper: 
