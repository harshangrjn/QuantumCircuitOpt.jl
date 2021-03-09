#--------------------------------------------------#
# Add all the constraints of the formulation here  #
#--------------------------------------------------#

function build_QCO_mip_constraints(m::QCOoptimizer)
    JuMP.@constraint(m.Qmodel, x[1] + x[2] <= 5)
    return
end
