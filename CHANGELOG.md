QuantumCircuitOpt.jl Change Log
===============================

### v0.1.4
- Optimizer time limit option added for early termination of solvers
- `log.jl` updates to support optimizer hitting time limits 
- Updated commuting matrix valid inequalites 
- Added support for eliminating pairs of M*M^-1 = M^-1*M = I
- Updated unit tests for the above changes

### v0.1.3
- Update to minimize_depth from maximizing IGate to minimizing non-IGates, which also improved run times. 
- Clean-up and update to unit tests to reflect the above change
- Fixed Gurobi solver bug in `examples/solver.jl`
- Involutory matrix constraints added and corresponding unit tests

### v0.1.2
- Update to warm start with an initial circuit (without U3 and R gates)
- Update to unit tests to reflect the above changes
- Update to docs to reflect the above changes
- Bug fixes in error/warn messages 

### v0.1.1
- docs/make.jl updates for sidebar name
- Updated docs logo to handle dark theme

### v0.1.0
- Initial function-based working implementation 
- Support for 2- and 3-qubit decompositions
- Added preliminary documentation 
- Added preliminary unit tests
- Github Actions set-up
