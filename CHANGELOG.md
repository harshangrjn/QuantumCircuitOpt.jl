QuantumCircuitOpt.jl Change Log
===============================

### v0.2.0
- MAJOR framework updates which support gates upto *five* qubit gates. Framework is now flexible to generalize it to even larger qubit circuits, by updating functions `kron_single_qubit_gate` and `kron_double_qubit_gate`
- Bug fixes and major re-factoring of `src/data.jl` to support any permutation of gates on qubits 
- Enabling Controlled-R (CRX, CRY, CRZ) and Controlled-U3 (CU3) gates
- Adding functions to data.jl that are the equivalent of the single-qubit ones (e.g. get_all_CR_gates)
- Enabling this functionality in the log.jl
- Kron symbol updated from `⊗` to `x`
- Enabling functionality to take elementary gates with kronecker symbols over all qubits. For example, `"H_1xCNot_23xI_4xI_5"`, `H_1xH_2xT_3`, `I_1xSX_2xTdagger_3`. This functionality only supports one and 2 qubit gates without angle parameters (excluded controlled R and U3 gates)
- Adding various new 2 qubit gates (including controlled R and U3) in `src/gates.jl`
- Expanding the elementary gates set with cleaner and flexible framework for any `num_qubits` in functions, `get_full_sized_gate` and `get_full_sized_kron_symbol_gate`
- Adding new functions: `_parse_gates_with_kron_symbol`, `kron_double_qubit_gate`
- Updating Docs and adding new unit tests to reflect above changes

### v0.1.9
- Updated toffoli input circuit to make it feasible for the commuting gate constraints 
- Added support for 3-qubit Toffoli using Rotation gates, Fredkin gate and `CNot_13` (in `src/examples`)
- Infeasibility bug fix in commuting gate constraints (for `T_2` expressed as `U3_2(0,0,π/4)`)

### v0.1.8
- Added support for specificity of qubits in U3 and Rotation (RX, RY, RZ) gates (Major change). Example: U3_1, RX_2, RZ_3, etc
- Added support for kronecker product gates (`X_1⊗X_2`, `S_1⊗S_2`, etc) and bug fix to handle this change in `get_compressed_decomposition` in `src/log.jl`
- Added `is_multi_qubit_gate` function in `src/log.jl`
- Docs and tests updated to reflect above changes

### v0.1.7
- Added support for a few controlled versions of V gates in 2 & 3 qubits
- Added support for Grover's diffusion operator
- Naming convention updated from NameQubit to Name_Qubit for 1 qubit gates (consistent with 2 qubit gates)
- `cnot_ij` updated to `CNot_ij` to be consistent with naming convention in `src/gates.jl`
- Docs and tests updated to reflect above changes

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
