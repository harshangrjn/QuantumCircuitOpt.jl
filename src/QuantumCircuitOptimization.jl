module QuantumCircuitOptimization

using JuMP
import LinearAlgebra
import Random
import Memento

using MathOptInterface
const MOI = MathOptInterface

# Create our module level logger (this will get precompiled)
const _LOGGER = Memento.getlogger(@__MODULE__)

# Register the module level logger at runtime so that folks can access the logger via `getlogger(QuantumCircuitOptimization)`
# NOTE: If this line is not included then the precompiled `QuantumCircuitOptimization._LOGGER` won't be registered at runtime.
__init__() = Memento.register(_LOGGER)

"Suppresses information and warning messages output by PowerModels, for fine grained control use the Memento package"
function silence()
    Memento.info(_LOGGER, "Suppressing information and warning messages for the rest of this session.  Use the Memento package for more fine-grained control of logging.")
    Memento.setlevel!(Memento.getlogger(QuantumCircuitOptimization), "error")
end

"alows the user to set the logging level without the need to add Memento"
function logger_config!(level)
    Memento.config!(Memento.getlogger("PowerModels"), level)
end

include("types.jl")
include("utility.jl")
include("relaxations.jl")
include("log.jl")

end # module
