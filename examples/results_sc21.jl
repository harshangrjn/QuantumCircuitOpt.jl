#= 
Results for the paper "QuantumCircuitOpt: An Open-source Frameworkfor Provably Optimal 
Quantum Circuit Design" can be reproduced by executing this file. 

Version of QCOpt: v0.3.0

Last updated: Oct 13, 2021
=#

import QuantumCircuitOpt as QCOpt
using JuMP
using Gurobi

include("gates_sc21.jl")

# Choosing a different MIP solver other than Gurobi can slow down the run times
qcm_optimizer = JuMP.optimizer_with_attributes(Gurobi.Optimizer, "Presolve" => 1)

#-----------------------------#
#      Results of Table-I     #
#-----------------------------#
table_I_gates = ["controlled_Z",
                 "controlled_V",
                 "controlled_H",
                 "magic_using_U3_CNot_2_1",
                 "iSwapGate",
                 "GroverDiffusion_using_U3"]

result_with_VC = Dict{String,Any}()
result_without_VC = Dict{String,Any}()
run_times_tab1 = zeros(length(table_I_gates),2)

k = 1
for gates in table_I_gates 
    params = getfield(Main, Symbol(gates))()
    
    model_options = Dict{Symbol, Any}(:all_valid_constraints => 0)

    global result_with_VC = QCOpt.run_QCModel(params, 
                                              qcm_optimizer;
                                              options = model_options)
    run_times_tab1[k,1] = result_with_VC["solve_time"]

    model_options = Dict{Symbol, Any}(:all_valid_constraints => -1)
    
    global result_without_VC = QCOpt.run_QCModel(params, 
                                                 qcm_optimizer;
                                                 options = model_options)
    run_times_tab1[k,2] = result_without_VC["solve_time"]
    global k += 1 
end

#--------------------------------------------#
#      Results of Figure 3 (Magic basis)     #
#--------------------------------------------#
result = Dict{String,Any}()
figure_3_gates = ["magic_using_SH_CNot_2_1",
                  "magic_using_U3_CNot_2_1",
                  "magic_using_SH_CNot_1_2",
                  "magic_using_U3_CNot_1_2",
                  ]

run_times_fig3 = zeros(length(figure_3_gates),1)

k = 1
for gates in figure_3_gates 
    params = getfield(Main, Symbol(gates))()

    global result = QCOpt.run_QCModel(params, qcm_optimizer)

    run_times_fig3[k,1] = result["solve_time"]
    global k += 1 
end

#------------------------------------------------#
#      Results of Figure 4 (Grover operator)     #
#------------------------------------------------#
result = Dict{String,Any}()
figure_4_gates = ["GroverDiffusion_using_Clifford",
                  "GroverDiffusion_using_U3"]

run_times_fig4 = zeros(length(figure_4_gates),1)

k = 1
for gates in figure_4_gates 
    params = getfield(Main, Symbol(gates))()

    global result = QCOpt.run_QCModel(params, qcm_optimizer)

    run_times_fig4[k,1] = result["solve_time"]
    global k += 1 
end

#------------------------------------------------#
#      Results of Figure 5 (larger circuits)     #
#------------------------------------------------#
result = Dict{String,Any}()
figure_5_gates = ["toffoli_with_controlled_gates",
                  "Fredkin",
                  "double_toffoli",
                  "quantum_fulladder",
                  ]

run_times_fig5 = zeros(length(figure_5_gates),1)

k = 1
for gates in figure_5_gates 
    params = getfield(Main, Symbol(gates))()

    global result = QCOpt.run_QCModel(params, qcm_optimizer)

    run_times_fig5[k,1] = result["solve_time"]
    global k += 1 
end