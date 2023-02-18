QuantumCircuitOpt.jl Change Log
===============================

### v0.5.5
- `optimal_global_phase` now recognizes commuting elementary gate pairs that commute, are redundant and idempotent up to a global phase
- Docs updated with the SC22 (IEEE) paper link
- Added `isapprox_global_phase` utility to evaluate equivalences of any two complex gates up to a global phase
- Updated docs and unit tests to reflect above changes
- Removed compatability on Pkg

### v0.5.4
- Added a generalized function to obtain controlled gates (`controlled_gate()`) in any number of qubits
- Cleaned-up controlled gates in `src/gates.jl`
- Updated docs to reflect above changes

### v0.5.3
- Minor update: SC22 publication added in docs
- README banner update

### v0.5.2
- Added dependency on `Pkg` package for logging purposes
- Renamed gates in `src/examples` folder

### v0.5.1
- New feature: Added support for compiling optimal circuits upto a global phase using linear constraints (`optimal_global_phase`)
- Deactivated determinant feasiblity check for global phase model
- Bug fix in evaluation of kron-symboled elementary gates (#63)
- `get_full_sized_gate` function now supports gates with kronecker operations
- Updated docs and unit tests to reflect above changes

### v0.5.0
- Minor: Fixed result `primal_status` issue in `log.jl`
- Added helper functions for obtaining fixed indices of `U_var` unitary variables to zeros or constant values. Dropped linearization constraints for fixed `U_var` variables
- Dropped support for `unit_magnitude_constraints`
- Reformulated quadratic objective function of `slack_var` variables into a linear outer-approximation - increases code coverage due to MILP reformulation 
- Fixed `slack_var` based on fixed `U_var` variables  
- Minor updates in error messages for kron symboled gates
- Changed Cbc test dependency to HiGHS (MIP) solver 
- Updated docs and unit tests to reflect above changes

### v0.4.1
- Added support for returning exact feasible decomposition (without optimality certificate) using `exact_feasible`
- Added support for unitary constraints using binary linearizations (speeds up feasiblity problems)
- Added `W_using_HSTCnot` in examples
- Added `CNOT_using_GR` in examples (for Rydberg atom array-based simulators)
- Added support for `CSGate`, `CSdaggerGate`, `CTGate`, `CTdaggerGate` and `SSwapGate` (sqrt(`SwapGate`))
- Dropped support for `CZRevGate` as it is invariant to qubit flip
- Updated README and docs with the Youtube link for the JuliaCon's presentation
- Updated `decomposition_type` to `exact_optimal` from `exact`
- Updated docs and unit tests to reflect above changes

### v0.4.0
- Added support for two angle param gates througout data and log
- Added support for Rotation gate (`RGate`) with two angle params
- Added support for multi-qubit gates with angle params
- Added support for Global R gate (`GRGate`) with two angle params
- Clean-up in `data.jl` and `utility.jl` to handle multi-qubit gates
- Decomposition of Pauli-X using global rotation added in 2-qubit examples
- Updated planar-hull cuts evaluation to account for repeated sets of extreme points (improved run times)
- Clean-up of `constraint_convex_hull_complex_gates` 
- Updated docs and unit tests to reflect above changes

### v0.3.6
- DOI link for publication added
- Added support for JuMP v1.0

### v0.3.5
- Dropped support for redundant `constraint_complex_to_real_symmetry_compact` function
- Added unit magnitude (outer-approximation) constraints for unitary variables
- Updated docs and unit tests to reflect above changes
- Updated README.md and docs with the Youtube link for the SC21 presentation
- Added support for JuMP v0.22+
- Dropped explicit support for MathOptInterface

### v0.3.4
- Updated error messages in `src/data.jl` to better represent one and two qubit gates
- More determinant-based tests (to handle complex det values) to detect infeasibility
- Updated unit tests to cover CNOT gate bound constraints

### v0.3.3 
- Added determinant-based test to detect infeasibility in the pre-processing step
- Added support for 3-qubit RCCXGate (relative Toffoli)
- Added support for 3-qubit MargolusGate (simplified Toffoli) - #46
- Added support for 3-qubit CiSwapGate (controlled iSwap)
- Updated docs and unit tests to reflect above changes

### v0.3.2
- Moved all the model options (including valid constraints) from `qc_model.jl` to `types.jl` into `QCModelOptions`, a struct form 
- Streamlined default options for `QCModelOptions`
- Moved `relax_integrality`, `time_limit`, `objective_slack_penalty` and `optimizer_log` options from `data.jl` to `QCModelOptions`
- Removed `identify_real_gates`, `optimizer` and `optimizer_presolve` options from `params` and updated `solvers.jl`
- Default recognition of real elementary and target gates, and implements a compact MIP
- Updated docs and unit tests to reflect above changes

### v0.3.1
- Added decomposition of the QFT2 gate using H, T and CNOT gates to `2qubits_gates.jl`
- Added CITATION.bib 
- Added QCOpt's framework to quick start guide in docs
- Minor updates in README

### v0.3.0
- Enabling support for convex hull evaluation (Graham's algorithm) within QCOpt (`QCO.convex_hull()`)
- Dropping dependency on QHull.jl
- Added decomposition of quantum full-adder, double-toffoli and double-peres gates to `4qubit_gates.jl`
- Added decomposition of toffoli gate using controlled gates and miller gate to `3qubit_gates.jl`
- Added decomposition of QFT2 gates using R gates to `2qubits_gates.jl`
- Issue #37 fixed and more meaningful messages for input gate parsing
- Updated unit tests to cover `convex_hull` function
- Bug fix in `_get_cnot_idx` to account for gates with kronecker symbols
- Inculded `CXGate` for evaluation of `_get_cnot_idx` in `data.jl`
- Separate files added in the `examples` folder for SC21 paper's results
- Options to add lower and upper bounds on number of CNOT gates in user-defined `params`
- User-specified "depth" changed to "maximum_depth"
- Updated docs and unit tests to reflect above changes

### v0.2.5
- Major updates in input gate convention for two qubit gates. For example, `CNot_12` becomes `CNot_1_2`. This update now makes the package flexible to be able to compute on any number of qubit circuits. `CNot_2_10` is a valid input for a 10 qubit circuit
- Major clean up in `data.jl` including generalized parsing in `get_full_sized_gate` and `get_full_sized_kron_gate` functions
- Clean up of `log.jl` to make it generic and consistent with above changes 
- Enabling support for convex hull-based cuts on complex variables at every depth (improved run times in CPLEX, default = turned off)
- Added `constraint_convex_hull_complex_gates` function in support with depndency of QHull.jl 
- Added `round_real_value` function and cleaned up `round_complex_values` in `utility.jl`
- Added toffoli-left block and toffoli-right block decompositions to `3qubit_gates.jl`
- Update `is_multi_qubit_gate` function to make it generic to any gate type

### v0.2.4 
- Enabling support for 2-qubit Sycamore gate, native to Google's quantum processor
- Citation for the SC21 paper (accepted)
- Clean-up in real to complex matrix conversion and vice-versa in `src/utility.jl` 
- Better rounding in `relaxation_bilinear` function
- Updated docs with more function references to reflect the above changes

### v0.2.3
- Enabling support for Identity-gate symmetry-breaking constraints - improved run times 
- Enabling support for `identify_real_gates` which applies a compact MIP formulation (turned-off default)
- Updated unit tests and docs to reflect the above changes
- Minor bug fix in `is_multi_qubit_gate` function to handle gates with kron symbols

### v0.2.2
- Bug fix in `is_multi_qubit_gate` function in `src/log.jl` to handle any 2 qubit gates
- Decomposition for 4-qubit CNot_41 added in `examples/4qubit_gates.jl`
- Moved involutory gate evaluation from `data.jl` to  a `get_involutory_gates` function in `utility.jl` (to be consistent with other similar functions). `constraint_involutory_gates` function is updated accordingly to reflect these changes
- Meaningful error messages in `data.jl` for input gates with missing continuous angle parameters
- Updated unit tests on involutory gates

### v0.2.1
- Framework updates which support upto *nine* qubit circuit decompositions
- Generalizing and condensing `kron_single_qubit_gate` function to handle any number of integer-valued qubits in `src/utility.jl`
- Generalizing and condensing `kron_two_qubit_gate` function to handle any number of integer-valued qubits without explicit enumeration in `src/utility.jl`. Also, computes slightly faster with lesser memory on larger qubits than the previous version
- Docs updates for missing inputs options
- Enabling support for PhaseGate with discretizations in `src/data.jl`
- Enabling support for DCXGate in `src/data.jl`
- Minor clean ups and removal of redundant functions in `src/data.jl`

### v0.2.0
- MAJOR framework updates which support gates upto *five* qubit gates. Framework is now flexible to generalize it to even larger qubit circuits, by updating functions `kron_single_qubit_gate` and `kron_two_qubit_gate`
- Bug fixes and major re-factoring of `src/data.jl` to support any permutation of gates on qubits 
- Enabling Controlled-R (CRX, CRY, CRZ) and Controlled-U3 (CU3) gates
- Adding functions to data.jl that are the equivalent of the single-qubit ones (e.g. get_all_CR_gates)
- Enabling this functionality in the log.jl
- Kron symbol updated from `⊗` to `x`
- Enabling functionality to take elementary gates with kronecker symbols over all qubits. For example, `"H_1xCNot_23xI_4xI_5"`, `H_1xH_2xT_3`, `I_1xSX_2xTdagger_3`. This functionality only supports one and 2 qubit gates without angle parameters (excluded controlled R and U3 gates)
- Adding various new 2 qubit gates (including controlled R and U3) in `src/gates.jl`
- Expanding the elementary gates set with cleaner and flexible framework for any `num_qubits` in functions, `get_full_sized_gate` and `get_full_sized_kron_gate`
- Adding new functions: `_parse_gates_with_kron_symbol`, `kron_two_qubit_gate`
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
