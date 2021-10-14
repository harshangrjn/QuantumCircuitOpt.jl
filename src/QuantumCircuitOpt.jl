module QuantumCircuitOpt

import JuMP
import LinearAlgebra
import Memento
import MathOptInterface

const MOI = MathOptInterface
const LA = LinearAlgebra
const QCO = QuantumCircuitOpt
const kron_symbol = 'x'
const qubit_separator = '_'

# Create our module level logger (this will get precompiled)
const _LOGGER = Memento.getlogger(@__MODULE__)

# Register the module level logger at runtime so that folks can access the logger via `getlogger(QuantumCircuitOpt)`
# NOTE: If this line is not included then the precompiled `QuantumCircuitOpt._LOGGER` won't be registered at runtime.
__init__() = Memento.register(_LOGGER)

"Suppresses information and warning messages output by QuantumCircuitOpt, for fine grained control use of the Memento package"
function silence()
    Memento.info(_LOGGER, "Suppressing information and warning messages for the rest of this session.  Use the Memento package for more fine-grained control of logging.")
    Memento.setlevel!(Memento.getlogger(QuantumCircuitOpt), "error")
end

"allows the user to set the logging level without the need to add Memento"
function logger_config!(level)
    Memento.config!(Memento.getlogger("QuantumCircuitOpt"), level)
end

include("data.jl")
include("gates.jl")
include("types.jl")
include("utility.jl")
include("chull.jl")
include("relaxations.jl")
include("qc_model.jl")
include("variables.jl")
include("constraints.jl")
include("objective.jl")
include("log.jl")
include("solution.jl")

end # module
