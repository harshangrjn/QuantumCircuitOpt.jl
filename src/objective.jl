#-------------------------------------------------------#
# Build the objective function of the formulation here  #
#-------------------------------------------------------#

function build_QCO_mip_objective(m::QCOoptimizer)
    JuMP.@objective(m.Qmodel, Min, x[1] + x[2])
    return
end

