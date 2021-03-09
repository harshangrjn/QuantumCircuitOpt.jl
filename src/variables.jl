#-------------------------------------------------------#
# Initialize all the variables of the formulation here  #
#-------------------------------------------------------#

function build_QCO_mip_variables(m::QCOoptimizer)
    JuMP.@variable(m.Qmodel, 1 <= x[1:2] <= 10)
    return
end


