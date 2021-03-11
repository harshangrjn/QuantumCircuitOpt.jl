#-------------------------------------------------------#
# Build the objective function of the formulation here  #
#-------------------------------------------------------#

function objective_QCOModel(m::QCOoptimizer)
    JuMP.@objective(m.Qmodel, Min, x[1] + x[2])
    return
end

