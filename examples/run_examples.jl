"""

Input data format for `QuantumCircuitOpt.jl`

`n_qubits`
    Number of qubits (>= 2)

`depth`
    Maximum depth allowed for decomposition (>= 2)

`elementary_gates`
    Input all the one and two qubit elementary gates here. Comprehensive list of gates currently supported in `QuantumCircuitOpt` 
can be found in `src/data.jl`.

`target_gate`
    Input the target gate which you wish to decompose using the above-mentioned elementary gates. 

`R_x_discretization`, `R_y_discretization`,`R_z_discretization`
    Input discretization angles for each of the `R_x`, `R_y` and `R_z` gates, respectively, which are part of the elementary gates. 

`U_θ_discretization`, `U_ϕ_discretization`, `U_λ_discretization`
    Input discretization angles for every angle `θ`, `ϕ` and `λ`, respectively, which are part of the U3 elementary gate. 

`initial_gate`
    This gate will form as an intitial condition to the decomposition. `Identity` will be the default gate for all the tests. 

`objective`
    `minimize_depth`: Minimizes the total depth of decomposition. For this option, include `Identity` matrix in the above-mentioned `elementary_gates`. 
    `minimize_cnot` : Minimizes the number of CNOT gates in the decomposition. 

`decomposition_type`
    `exact`      : `QuantumCircuitOpt` finds an exact decomposition if it is a feasible solution. 
    `approximate`: `QuantumCircuitOpt` finds an approximate decomposition if an exact one does not exist; otherwise it will return an exact solution.

`optimizer`
    Choose the mixed-integer optimizer here. For multiple solver options, check `examples/solver.jl`. 
    Note that CPLEX (`cplex`) or Gurobi (`gurobi`) will perform numerically best.

`presolve`
    This option enables and disables the presolve option in the above-mentioned `optimizer`. Turning it off can lead to slower run times. 

`optimizer_log`
    This option enables and disables console logging for the `optimizer`.

`relax_integrality`
    This option transforms integer variables into continuous variables. 

"""

using QuantumCircuitOpt
using JuMP
using CPLEX
using LinearAlgebra
using Cbc

include("solver.jl")
include("2qubit_gates.jl")


gate_params = [test_hadamard(), 
               test_S(), 
               test_magic_M(), 
               test_controlled_R2(), 
               test_cnot_21(), 
               test_cnot_21_with_U(), 
               test_controlled_Z(),
               test_magic_M_using_SHCnot(),
               test_swap(),
               test_controlled_V(),
               test_W(),
               test_W_using_HCnot()]

#------------------------------#
#      Optimization model      #
#------------------------------#
function run_QuantumCircuitOpt(params)

    qcm_optimizer = get_solver(params)
    data = QuantumCircuitOpt.get_data(params)

    model_qc  = QuantumCircuitOpt.build_QCModel(data, model_type = "compact_formulation", commute_matrix_cuts = false)
    result_qc = QuantumCircuitOpt.optimize_QCModel!(model_qc, optimizer = qcm_optimizer)
    QuantumCircuitOpt.visualize_QCModel_solution(result_qc, data)

end

for params in gate_params 
    run_QuantumCircuitOpt(params)
end

