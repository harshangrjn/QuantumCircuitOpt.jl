QuantumCircuitOpt.jl Change Log
===============================

### v0.1.6
- Added support for obtaining idempotent matrices and adding corresponding valid constraints 
- Updated `all_valid_constraints` from `Bool` to values in `[-1,0,1]`
- `qc_model.jl` clean-up 
- Updated unit tests to reflect above changes

### v0.1.5
- Added support for pair-wise commuting products with Identity, (AB=BA=I)
- Bug fix to remove redundant pairs in `get_commutative_gate_pairs` when corresponding product is an Identity gate
- Added support for eliminating redundant pairs whose product matches an input elementary gate (`get_redundant_gate_product_pairs` in utility) - improved run times
- Added support to turn on/off all valid constraints (`all_valid_constraints`)
- Updated docs to reflect the improved run times

### v0.1.4
- Optimizer time limit option added for early termination of solvers
- `log.jl` updates to support optimizer hitting time limits 
- Revamping of commuting gate valid constraints (dropped support for triplets)
- Revamping data.jl to fix default values of user inputs
- Updated unit tests for the above changes

### v0.1.3
- Update to minimize_depth from maximizing IGate to minimizing non-IGates, which also improved run times. 
- Clean-up and update to unit tests to reflect the above change
- Fixed Gurobi solver bug in `examples/solver.jl`
- Involutory gate constraints added and corresponding unit tests

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