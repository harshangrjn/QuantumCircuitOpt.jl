using JuMP
using CPLEX
using QuantumCircuitOptimization
using Ipopt
using Juniper

# MIP solvers
const cplex = optimizer_with_attributes(CPLEX.Optimizer, 
                                        MOI.Silent() => true, 
                                        "CPX_PARAM_PREIND" => false) 

# Local solvers
const ipopt = optimizer_with_attributes(Ipopt.Optimizer, 
                                        MOI.Silent() => true, 
                                        "sb" => "yes", 
                                        "max_iter" => 9999)

const juniper = optimizer_with_attributes(Juniper.Optimizer, 
                                          MOI.Silent() => true, 
                                          "mip_solver" => cplex, 
                                          "nl_solver" => ipopt)

const qco = optimizer_with_attributes(QuantumCircuitOptimization.QCOoptimizer, 
                                        "nlp_solver" => ipopt,
                                        "minlp_solver" => juniper,  
                                        "mip_solver" => cplex)
                                        
m = Model(qco)
JuMP.optimize!(m)
@show objective_value(m)                                        
            
