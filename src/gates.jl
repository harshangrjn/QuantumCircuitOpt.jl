#=
IMPORTANT: Update the lists below whenever a 1 or 2 qubit gate is added to this file. 
           Update QCO._get_angle_gates_idx in src/data.jl if gates with angle parameters are added to this file.
=#

const ONE_QUBIT_GATES_ANGLE_PARAMETERS     = ["U3", "U2", "U1", "RX", "RY", "RZ", "Phase"]

const ONE_QUBIT_GATES_CONSTANTS            = ["Identity", "I", "H", "X", "Y", "Z", "S", 
                                              "Sdagger", "T", "Tdagger", "SX", "SXdagger"]

const TWO_QUBIT_GATES_ANGLE_PARAMETERS     = ["CRX", "CRXRev", "CRY", "CRYRev", "CRZ", "CRZRev", 
                                              "CU3", "CU3Rev"]

# Gates invariant to qubit flip
const TWO_QUBIT_GATES_CONSTANTS_SYMMETRIC  = ["Swap", "iSwap", "Sycamore", "DCX", "W", "M", 
                                              "QFT2", "HCoin", "GroverDiffusion"]

# Gates non-invariant to qubit flip
const TWO_QUBIT_GATES_CONSTANTS_ASYMMETRIC = ["CNot", "CNotRev", "CX", "CXRev", "CY", "CYRev", "CZ", 
                                              "CZRev", "CH", "CHRev", "CV", "CVRev", "CVdagger", 
                                              "CVdaggerRev", "CSX", "CSXRev"]

const ONE_QUBIT_GATES           = union(QCO.ONE_QUBIT_GATES_CONSTANTS, 
                                        QCO.ONE_QUBIT_GATES_ANGLE_PARAMETERS)

const TWO_QUBIT_GATES_CONSTANTS = union(QCO.TWO_QUBIT_GATES_CONSTANTS_SYMMETRIC, 
                                        QCO.TWO_QUBIT_GATES_CONSTANTS_ASYMMETRIC)

const TWO_QUBIT_GATES           = union(QCO.TWO_QUBIT_GATES_CONSTANTS, 
                                        QCO.TWO_QUBIT_GATES_ANGLE_PARAMETERS)

#----------------------------------------#
#            Single-qubit gates          #
#----------------------------------------#

@doc raw"""
    IGate(num_qubits::Int64)

Identity matrix for an input number of qubits.

**Matrix Representation (num_qubits = 1)**

```math
I = \begin{pmatrix}
        1 & 0 \\
        0 & 1
    \end{pmatrix}
```
"""
function IGate(num_qubits::Int64)

    return Array{Complex{Float64},2}(Matrix(LA.I, 2^num_qubits, 2^num_qubits))

end

@doc raw"""
    U3Gate(θ::Number, ϕ::Number, λ::Number)

Universal single-qubit rotation gate with three Euler angles, ``\theta``, ``\phi`` and ``\lambda``.

**Matrix Representation**

```math
\newcommand{\th}{\frac{\theta}{2}}

U3(\theta, \phi, \lambda) =
    \begin{pmatrix}
        \cos(\th)          & -e^{i\lambda}\sin(\th) \\
        e^{i\phi}\sin(\th) & e^{i(\phi+\lambda)}\cos(\th)
    \end{pmatrix}
```
"""
function U3Gate(θ::Number, ϕ::Number, λ::Number)

    QCO._verify_θ_bounds(θ)
    QCO._verify_ϕ_bounds(ϕ)
    QCO._verify_λ_bounds(λ)

    U3 = Array{Complex{Float64},2}([           cos(θ/2)               -(cos(λ) + (sin(λ))im)*sin(θ/2) 
                                    (cos(ϕ) + (sin(ϕ))im)*sin(θ/2)  (cos(λ+ϕ) + (sin(λ+ϕ))im)*cos(θ/2)])

    return QCO.round_complex_values(U3)
end

@doc raw"""
    U2Gate(ϕ::Number, λ::Number)

Universal single-qubit rotation gate with two Euler angles, ``\phi`` and ``\lambda``. U2Gate is the special case of 
[U3Gate](@ref). 

**Matrix Representation**

```math
U2(\phi, \lambda) = \frac{1}{\sqrt{2}}
\begin{pmatrix}
    1          & -e^{i\lambda} \\
    e^{i\phi} & e^{i(\phi+\lambda)}
\end{pmatrix}
```
"""
function U2Gate(ϕ::Number, λ::Number)
    
    θ = π/2

    U2 = QCO.U3Gate(θ, ϕ, λ)    
    return U2
end

@doc raw"""
    U1Gate(λ::Number)

Universal single-qubit rotation gate with one Euler angle, ``\lambda``. U1Gate represents rotation about the Z axis and 
is the special case of [U3Gate](@ref), which also known as the [PhaseGate](@ref). Also note that ``U1(\pi) = ``[ZGate](@ref), ``U1(\pi/2) = ``[SGate](@ref) and 
``U1(\pi/4) = ``[TGate](@ref).

**Matrix Representation**

```math
U1(\lambda) =
\begin{pmatrix}
    1 & 0 \\
    0 & e^{i\lambda}
\end{pmatrix}
```
"""
function U1Gate(λ::Number)

    θ = 0
    ϕ = 0

    U1 = QCO.U3Gate(θ, ϕ, λ)
    return U1
end

@doc raw"""
    RXGate(θ::Number)

A single-qubit Pauli gate which represents rotation about the X axis.

**Matrix Representation**
```math
\newcommand{\th}{\frac{\theta}{2}}

RX(\theta) = exp(-i \th X) =
    \begin{pmatrix}
        \cos{\th}   & -i\sin{\th} \\
        -i\sin{\th} & \cos{\th}
    \end{pmatrix}
```
"""
function RXGate(θ::Number)

    QCO._verify_θ_bounds(θ)

    RX = Array{Complex{Float64},2}([cos(θ/2) -(sin(θ/2))im; -(sin(θ/2))im cos(θ/2)])
    
    return QCO.round_complex_values(RX)
end

@doc raw"""
    RYGate(θ::Number)

A single-qubit Pauli gate which represents rotation about the Y axis.

**Matrix Representation**
```math
\newcommand{\th}{\frac{\theta}{2}}

RY(\theta) = exp(-i \th Y) =
    \begin{pmatrix}
        \cos{\th} & -\sin{\th} \\
        \sin{\th} & \cos{\th}
    \end{pmatrix}
```
"""
function RYGate(θ::Number)

    QCO._verify_θ_bounds(θ)

    RY = Array{Complex{Float64},2}([cos(θ/2) -(sin(θ/2)); (sin(θ/2)) cos(θ/2)])
    
    return QCO.round_complex_values(RY)
end

@doc raw"""
    RZGate(θ::Number)

A single-qubit Pauli gate which represents rotation about the Z axis. This gate is also equivalent to [U1Gate](@ref) up to a phase factor, 
that is, ``RZ(\theta) = e^{-i{\theta}/2}U1(\theta)``.

**Matrix Representation**
```math
\newcommand{\th}{\frac{\theta}{2}}

RZ(\theta) = exp(-i\th Z) =
\begin{pmatrix}
    e^{-i\th} & 0 \\
    0 & e^{i\th}
\end{pmatrix}
```
"""
function RZGate(θ::Number)

    QCO._verify_θ_bounds(θ)

    RZ = Array{Complex{Float64},2}([(cos(θ/2) - (sin(θ/2))im) 0; 0 (cos(θ/2) + (sin(θ/2))im)])
    
    return QCO.round_complex_values(RZ)
end

@doc raw"""
    HGate(num_qubits::Int64)

Single-qubit Hadamard gate, which is a ``\pi`` rotation about the X+Z axis, thus equivalent to [U3Gate](@ref)(``\frac{\pi}{2},0,\pi``)

**Matrix Representation**

```math
H = \frac{1}{\sqrt{2}}
        \begin{pmatrix}
            1 & 1 \\
            1 & -1
        \end{pmatrix}
```
"""
function HGate()

    return Array{Complex{Float64},2}(1/sqrt(2)*[1 1; 1 -1])

end

@doc raw"""
    XGate()

Single-qubit Pauli-X gate (``\sigma_x``), equivalent to [U3Gate](@ref)(``\pi,0,\pi``)

**Matrix Representation**

```math
X = \begin{pmatrix}
0 & 1 \\
1 & 0
\end{pmatrix}
```
"""
function XGate()

    return Array{Complex{Float64},2}([0 1; 1 0])

end

@doc raw"""
    YGate()

Single-qubit Pauli-Y gate (``\sigma_y``), equivalent to [U3Gate](@ref)(``\pi,\frac{\pi}{2},\frac{\pi}{2}``)

**Matrix Representation**

```math
Y = \begin{pmatrix}
0 & -i \\
i & 0
\end{pmatrix}
```
"""
function YGate()

    return Array{Complex{Float64},2}([0 -im; im 0])

end

@doc raw"""
    ZGate()

Single-qubit Pauli-Z gate (``\sigma_z``), equivalent to [U3Gate](@ref)(``0,0,\pi``)

**Matrix Representation**

```math
Z = \begin{pmatrix}
1 & 0 \\
0 & -1
\end{pmatrix}
```
"""
function ZGate()

    return Array{Complex{Float64},2}([1 0; 0 -1])

end

@doc raw"""
    SGate()

Single-qubit S gate, equivalent to [U3Gate](@ref)(``0,0,\frac{\pi}{2}``). This 
gate is also referred to as a Clifford gate, P gate or a square-root of Pauli-[ZGate](@ref). Historically, this is also 
called as the phase gate (denoted by P), since it shifts the phase of the one state relative to the zero state.

**Matrix Representation**

```math
S = \begin{pmatrix}
1 & 0 \\
0 & i
\end{pmatrix}
```
"""
function SGate()

    return Array{Complex{Float64},2}([1 0; 0 im])

end

@doc raw"""
    SdaggerGate()

Single-qubit, hermitian conjugate of the [SGate](@ref). This is also an alternative square root of 
the [ZGate](@ref). 

**Matrix Representation**

```math
S = \begin{pmatrix}
1 & 0 \\
0 & -i
\end{pmatrix}
```
"""
function SdaggerGate()

    return Array{Complex{Float64},2}([1 0; 0 -im])

end

@doc raw"""
    TGate()

Single-qubit T gate, equivalent to [U3Gate](@ref)(``0,0,\frac{\pi}{4}``). This 
gate is also referred to as a ``\frac{\pi}{8}`` gate or as a fourth-root of Pauli-[ZGate](@ref). 

**Matrix Representation**

```math
T = \begin{pmatrix}
1 & 0 \\
0 & e^{i\pi/4}
\end{pmatrix}
```
"""
function TGate()

    return Array{Complex{Float64},2}([1 0; 0 (1/sqrt(2)) + (1/sqrt(2))im])

end

@doc raw"""
    TdaggerGate()

Single-qubit, hermitian conjugate of the [TGate](@ref). This gate is equivalent to [U3Gate](@ref)(``0,0,-\frac{\pi}{4}``). This 
gate is also referred to as the fourth-root of Pauli-[ZGate](@ref). 

**Matrix Representation**

```math
T^{\dagger} = \begin{pmatrix}
1 & 0 \\
0 & e^{-i\pi/4}
\end{pmatrix}
```
"""
function TdaggerGate()

    return Array{Complex{Float64},2}([1 0; 0 (1/sqrt(2)) - (1/sqrt(2))im])

end

@doc raw"""
    SXGate()

Single-qubit square root of pauli-[XGate](@ref).

**Matrix Representation**

```math
\sqrt{X} = \frac{1}{2} \begin{pmatrix}
1 + i & 1 - i \\
1 - i & 1 + i
\end{pmatrix}
```
"""
function SXGate()

    return Array{Complex{Float64},2}(1/2*[1+im 1-im; 1-im 1+im])

end

@doc raw"""
    SXdaggerGate()

Single-qubit hermitian conjugate of the square root of pauli-[XGate](@ref), or the [SXGate](@ref).

**Matrix Representation**

```math
\sqrt{X}^{\dagger} = \frac{1}{2} \begin{pmatrix}
1 - i & 1 + i \\
1 + i & 1 - i
\end{pmatrix}
```
"""
function SXdaggerGate()

    return Array{Complex{Float64},2}(1/2*[1-im 1+im; 1+im 1-im])

end

@doc raw"""
    PhaseGate()

Single-qubit rotation gate about the Z axis. This is also equivalent to [U3Gate](@ref)(``0,0,\lambda``). This 
gate is also referred to as the [U1Gate](@ref). 

**Matrix Representation**

```math
P(\lambda) = \begin{pmatrix}
    1 & 0 \\
    0 & e^{i\lambda}
\end{pmatrix}
```
"""
function PhaseGate(λ::Number)
    return QCO.U1Gate(λ)
end

#-------------------------------------#
#            Two-qubit gates          #
#-------------------------------------#
@doc raw"""
    CNotGate()

Two-qubit controlled NOT gate with control and target on first and second qubits, respectively. This is also 
called the controlled X gate ([CXGate](@ref)). 

**Circuit Representation**
```
q_0: ──■──
     ┌─┴─┐
q_1: ┤ X ├
     └───┘
```

**Matrix Representation**

```math
CNot = \begin{pmatrix}
    1 & 0 & 0 & 0 \\
    0 & 1 & 0 & 0 \\
    0 & 0 & 0 & 1 \\
    0 & 0 & 1 & 0
    \end{pmatrix}
```
"""
function CNotGate()

    return Array{Complex{Float64},2}([1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]) 

end

@doc raw"""
    CNotRevGate()

Two-qubit reverse controlled NOT gate, with target and control on first and second qubits, respectively. 

**Circuit Representation**
```
     ┌───┐
q_0: ┤ X ├
     └─┬─┘
q_1: ──■──
```

**Matrix Representation**

```math
CNotRev = \begin{pmatrix}
            1 & 0 & 0 & 0 \\
            0 & 0 & 0 & 1 \\
            0 & 0 & 1 & 0 \\
            0 & 1 & 0 & 0
            \end{pmatrix}
```
"""
function CNotRevGate()

    return Array{Complex{Float64},2}([1 0 0 0; 0 0 0 1; 0 0 1 0; 0 1 0 0]) 

end

@doc raw"""
    DCXGate()

Two-qubit double controlled NOT gate consisting of two back-to-back [CNotGate](@ref)s with alternate controls. 

**Circuit Representation**
```
          ┌───┐
q_0: ──■──┤ X ├
     ┌─┴─┐└─┬─┘
q_1: ┤ X ├──■──
     └───┘
```

**Matrix Representation**

```math
DCX = \begin{pmatrix}
    1 & 0 & 0 & 0 \\
    0 & 0 & 0 & 1 \\
    0 & 1 & 0 & 0 \\
    0 & 0 & 1 & 0
    \end{pmatrix}
```
"""
function DCXGate()

    return Array{Complex{Float64},2}([1 0 0 0; 0 0 0 1; 0 1 0 0; 0 0 1 0]) 

end

@doc raw"""
    CXGate()

Two-qubit controlled [XGate](@ref), which is also the same as [CNotGate](@ref). 

**Circuit Representation**
```
q_0: ──■──
     ┌─┴─┐
q_1: ┤ X ├
     └───┘
```

**Matrix Representation**

```math
CX = \begin{pmatrix}
    1 & 0 & 0 & 0 \\
    0 & 1 & 0 & 0 \\
    0 & 0 & 0 & 1 \\
    0 & 0 & 1 & 0
    \end{pmatrix}
```
"""
function CXGate()

    return Array{Complex{Float64},2}([1 0 0 0; 0 1 0 0; 0 0 0 1; 0 0 1 0]) 

end

@doc raw"""
    CXRevGate()

Two-qubit reverse controlled-X gate, with target and control on first and second qubits, respectively. 
This is also the same as [CNotRevGate](@ref). 

**Circuit Representation**
```
     ┌───┐
q_0: ┤ X ├
     └─┬─┘
q_1: ──■──
```

**Matrix Representation**

```math
CXRev = I \otimes |0 \rangle\langle 0| + X \otimes |1 \rangle\langle 1| = \begin{pmatrix}
        1 & 0 & 0 & 0 \\
        0 & 0 & 0 & 1 \\
        0 & 0 & 1 & 0 \\
        0 & 1 & 0 & 0
    \end{pmatrix}
```
"""
function CXRevGate()
    
    # I ⊗ |0⟩⟨0|
    control_0 = kron(QCO.IGate(1), Array{Complex{Float64},2}([1 0; 0 0]))
    # X ⊗ |1⟩⟨1| 
    control_1 = kron(QCO.XGate(), Array{Complex{Float64},2}([0 0; 0 1]))
    
    return control_0 + control_1
end

@doc raw"""
    CYGate()

Two-qubit controlled [YGate](@ref). 

**Circuit Representation**
```
q_0: ──■──
     ┌─┴─┐
q_1: ┤ Y ├
     └───┘
```

**Matrix Representation**

```math
CY = |0 \rangle\langle 0| \otimes I + |1 \rangle\langle 1| \otimes Y = \begin{pmatrix}
    1 & 0 & 0 & 0 \\
    0 & 1 & 0 & 0 \\
    0 & 0 & 0 & -i \\
    0 & 0 & i & 0
    \end{pmatrix}
```
"""
function CYGate()

    # |0⟩⟨0| ⊗ I
    control_0 = kron(Array{Complex{Float64},2}([1 0; 0 0]) , QCO.IGate(1))
    # |1⟩⟨1| ⊗ Y 
    control_1 = kron(Array{Complex{Float64},2}([0 0; 0 1]) , QCO.YGate())
    
    return control_0 + control_1
end

@doc raw"""
    CYRevGate()

Two-qubit reverse controlled-Y gate, with target and control on first and second qubits, respectively. 

**Circuit Representation**
```
     ┌───┐
q_0: ┤ Y ├
     └─┬─┘
q_1: ──■──
```

**Matrix Representation**

```math
CYRev = I \otimes |0 \rangle\langle 0| + Y \otimes |1 \rangle\langle 1| = \begin{pmatrix}
        1 & 0 & 0 & 0 \\
        0 & 0 & 0 & -i \\
        0 & 0 & 1 & 0 \\
        0 & i & 0 & 0
    \end{pmatrix}
```
"""
function CYRevGate()
    
    # I ⊗ |0⟩⟨0|
    control_0 = kron(QCO.IGate(1), Array{Complex{Float64},2}([1 0; 0 0]))
    # Y ⊗ |1⟩⟨1| 
    control_1 = kron(QCO.YGate(), Array{Complex{Float64},2}([0 0; 0 1]))
    
    return control_0 + control_1
end

@doc raw"""
    CZGate()

Two-qubit, symmetric, controlled [ZGate](@ref). 

**Circuit Representation**
```
q_0: ──■──     ─■─
     ┌─┴─┐  ≡   │
q_1: ┤ Z ├     ─■─
     └───┘
```

**Matrix Representation**

```math
CZ = |0 \rangle\langle 0| \otimes I + |1 \rangle\langle 1| \otimes Z = \begin{pmatrix}
    1 & 0 & 0 & 0 \\
    0 & 1 & 0 & 0 \\
    0 & 0 & 1 & 0 \\
    0 & 0 & 0 & -1
    \end{pmatrix}
```
"""
function CZGate()
    
    # |0⟩⟨0| ⊗ I
    control_0 = kron(Array{Complex{Float64},2}([1 0; 0 0]) , QCO.IGate(1))
    # |1⟩⟨1| ⊗ Z 
    control_1 = kron(Array{Complex{Float64},2}([0 0; 0 1]) , QCO.ZGate())
    
    return control_0 + control_1 
end

@doc raw"""
    CZRevGate()

Two-qubit reverse controlled-Z gate, with target and control on first and second qubits, respectively. 

**Circuit Representation**
```
     ┌───┐
q_0: ┤ Z ├
     └─┬─┘
q_1: ──■──
```

**Matrix Representation**

```math
CZRev = I \otimes |0\rangle\langle0| + Z \otimes |1\rangle\langle1| = \begin{pmatrix}
    1 & 0 & 0 & 0 \\
    0 & 1 & 0 & 0 \\
    0 & 0 & 1 & 0 \\
    0 & 0 & 0 & -1
    \end{pmatrix}
```
"""
function CZRevGate()
    
    # I ⊗ |0⟩⟨0|
    control_0 = kron(QCO.IGate(1), Array{Complex{Float64},2}([1 0; 0 0]))
    # Z ⊗ |1⟩⟨1| 
    control_1 = kron(QCO.ZGate(), Array{Complex{Float64},2}([0 0; 0 1]))
    
    return control_0 + control_1
end

@doc raw"""
    CHGate()

Two-qubit, symmetric, controlled Hadamard gate ([HGate](@ref)). 

**Circuit Representation**
```
q_0: ──■──
     ┌─┴─┐  
q_1: ┤ H ├    
     └───┘
```

**Matrix Representation**

```math
CH = |0\rangle\langle 0| \otimes I + |1\rangle\langle 1| \otimes H = \begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \\
0 & 0 & \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{2}} \\
0 & 0 & \frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{2}}
\end{pmatrix}
```
"""
function CHGate()

    # |0⟩⟨0| ⊗ I
    control_0 = kron(Array{Complex{Float64},2}([1 0; 0 0]) , QCO.IGate(1))
    # |1⟩⟨1| ⊗ H 
    control_1 = kron(Array{Complex{Float64},2}([0 0; 0 1]) , QCO.HGate())
    
    return control_0 + control_1 
end

@doc raw"""
    CHRevGate()

Two-qubit reverse controlled-H gate, with target and control on first and second qubits, respectively. 

**Circuit Representation**
```
     ┌───┐
q_0: ┤ H ├
     └─┬─┘
q_1: ──■──
```

**Matrix Representation**

```math
CHRev = I \otimes |0\rangle\langle 0| + H \otimes |1\rangle\langle 1| = \begin{pmatrix}
        1 & 0 & 0 & 0 \\
        0 & \frac{1}{\sqrt{2}} & 0 & \frac{1}{\sqrt{2}} \\
        0 & 0 & 1 & 0 \\
        0 & \frac{1}{\sqrt{2}} & 0 & -\frac{1}{\sqrt{2}}
    \end{pmatrix}
```
"""
function CHRevGate()
    
    # I ⊗ |0⟩⟨0|
    control_0 = kron(QCO.IGate(1), Array{Complex{Float64},2}([1 0; 0 0]))
    # H ⊗ |1⟩⟨1| 
    control_1 = kron(QCO.HGate(), Array{Complex{Float64},2}([0 0; 0 1]))
    
    return control_0 + control_1
end

@doc raw"""
    CVGate()

Two-qubit, controlled-V gate, which is also the same as Controlled square-root of X gate ([CSXGate](@ref)).  

**Circuit Representation**
```
q_0: ──■──     
     ┌─┴─┐    
q_1: ┤ V ├     
     └───┘
```

**Matrix Representation**

```math
CV = \begin{pmatrix}
        1 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 \\
        0 & 0 & 0.5+0.5i & 0.5-0.5i \\
        0 & 0 & 0.5-0.5i & 0.5+0.5i
    \end{pmatrix}
```
"""
function CVGate()

    return QCO.CSXGate() 

end


@doc raw"""
    CVRevGate()

Two-qubit reverse controlled-V gate, with target and control on first and second qubits, respectively. 

**Circuit Representation**
```
     ┌───┐
q_0: ┤ V ├
     └─┬─┘
q_1: ──■──
```

**Matrix Representation**

```math
CVRev = \begin{pmatrix}
        1 & 0 & 0 & 0 \\
        0 & 0.5+0.5i & 0 & 0.5-0.5i \\
        0 & 0 & 1 & 0 \\
        0 & 0.5-0.5i & 0 & 0.5+0.5i
    \end{pmatrix}
```
"""
function CVRevGate()

    return Array{Complex{Float64},2}([1 0 0 0; 0 0.5+0.5im 0 0.5-0.5im; 0 0 1 0; 0 0.5-0.5im 0 0.5+0.5im])

end

@doc raw"""
    CVdaggerGate()

Two-qubit hermitian conjugate of controlled-V gate, which is also the same as hermitian conjugate Controlled square-root of X gate ([CSXGate](@ref)).  

**Circuit Representation**
```
q_0: ──■──     
     ┌─┴─┐    
q_1: ┤ V'├     
     └───┘
```

**Matrix Representation**

```math
CVdagger = \begin{pmatrix}
        1 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 \\
        0 & 0 & 0.5-0.5i & 0.5+0.5i \\
        0 & 0 & 0.5+0.5i & 0.5-0.5i
    \end{pmatrix}
```
"""
function CVdaggerGate()

    return Array{Complex{Float64},2}([1 0 0 0; 0 1 0 0; 0 0 0.5-0.5im 0.5+0.5im; 0 0 0.5+0.5im 0.5-0.5im])

end


@doc raw"""
    CVdaggerRevGate()

Two-qubit hermitian conjugate of reverse controlled-V gate, with target and control on first and second qubits, respectively. 

**Circuit Representation**
```
     ┌───┐
q_0: ┤ V'├
     └─┬─┘
q_1: ──■──
```

**Matrix Representation**

```math
CVdaggerRev = \begin{pmatrix}
        1 & 0 & 0 & 0 \\
        0 & 0.5-0.5i & 0 & 0.5+0.5i \\
        0 & 0 & 1 & 0 \\
        0 & 0.5+0.5i & 0 & 0.5-0.5i
    \end{pmatrix}
```
"""
function CVdaggerRevGate()

    return Array{Complex{Float64},2}([1 0 0 0; 0 0.5-0.5im 0 0.5+0.5im; 0 0 1 0; 0 0.5+0.5im 0 0.5-0.5im])

end

@doc raw"""
    WGate()

Two-qubit, W hermitian gate, typically useful to diagonlize the ([SwapGate](@ref)).  

**Matrix Representation**

```math
W = \begin{pmatrix}
        1 & 0 & 0 & 0 \\
        0 & \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{2}} & 0 \\
        0 & \frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{2}} & 0 \\
        0 & 0 & 0 & 1
    \end{pmatrix}
```
"""
function WGate()

    return Array{Complex{Float64},2}([1 0 0 0; 0 1/sqrt(2) 1/sqrt(2) 0; 0 1/sqrt(2) -1/sqrt(2) 0; 0 0 0 1])

end

@doc raw"""
    CRXGate(θ::Number)

Two-qubit controlled [RXGate](@ref). 

**Circuit Representation**
```
q_0: ────■────
     ┌───┴───┐
q_1: ┤ RX(ϴ) ├
     └───────┘
```

**Matrix Representation**

```math
\newcommand{\th}{\frac{\theta}{2}}

CRX(\theta)\ q_1, q_0 =
|0\rangle\langle0| \otimes I + |1\rangle\langle1| \otimes RX(\theta) =
    \begin{pmatrix}
        1 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 \\
        0 & 0 & \cos{\th}   & -i\sin{\th} \\
        0 & 0 & -i\sin{\th} & \cos{\th}
    \end{pmatrix}
```
"""
function CRXGate(θ::Number)

    CRX = Array{Complex{Float64},2}([ 1 0 0 0            
                                      0 1 0 0       
                                      0 0  cos(θ/2)  -(sin(θ/2))im 
                                      0 0  -(sin(θ/2))im  cos(θ/2)])

    return QCO.round_complex_values(CRX)
end

@doc raw"""
    CRXRevGate(θ::Number)

Two-qubit controlled reverse [RXGate](@ref). 

**Circuit Representation**
```
     ┌───────┐
q_1: ┤ RX(ϴ) ├
     └───┬───┘
q_0: ────■────
```

**Matrix Representation**

```math
\newcommand{\th}{\frac{\theta}{2}}

CRXRev(\theta)\ q_1, q_0 =
|0\rangle\langle0| \otimes I + |1\rangle\langle1| \otimes RX(\theta) =
    \begin{pmatrix}
        1 & 0 & 0 & 0 \\
        0 & \cos{\th} & 0 & -i\sin{\th} \\
        0 & 0 & 1 & 0\\
        0 & -i\sin{\th} & 0 & \cos{\th}
    \end{pmatrix}
```
"""
function CRXRevGate(θ::Number)

    CRXRev = Array{Complex{Float64},2}([ 1 0 0 0            
                                      0 cos(θ/2) 0 -(sin(θ/2))im      
                                      0 0 1 0
                                      0 -(sin(θ/2))im 0 cos(θ/2)])

    return QCO.round_complex_values(CRXRev)
end

@doc raw"""
    CRYGate(θ::Number)

Two-qubit controlled [RYGate](@ref). 

**Circuit Representation**
```
q_0: ────■────
     ┌───┴───┐
q_1: ┤ RY(ϴ) ├
     └───────┘
```

**Matrix Representation**

```math
\newcommand{\th}{\frac{\theta}{2}}

CRY(\theta)\ q_1, q_0 =
|0\rangle\langle0| \otimes I + |1\rangle\langle1| \otimes RY(\theta) =
    \begin{pmatrix}
        1 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 \\
        0 & 0 & \cos{\th}   & -\sin{\th} \\
        0 & 0 & \sin{\th} & \cos{\th}
    \end{pmatrix}
```
"""
function CRYGate(θ::Number)
    
    QCO._verify_θ_bounds(θ)

    CRY = Array{Complex{Float64},2}([ 1 0 0 0            
                                      0 1 0 0       
                                      0 0  cos(θ/2)  -(sin(θ/2)) 
                                      0 0  (sin(θ/2))  cos(θ/2)])

    return QCO.round_complex_values(CRY)
end

@doc raw"""
    CRYRevGate(θ::Number)

Two-qubit controlled reverse [RYGate](@ref). 

**Circuit Representation**
```
     ┌───────┐
q_1: ┤ RY(ϴ) ├
     └───┬───┘
q_0: ────■────
```

**Matrix Representation**

```math
\newcommand{\th}{\frac{\theta}{2}}

CRYRev(\theta)\ q_1, q_0 =
|0\rangle\langle0| \otimes I + |1\rangle\langle1| \otimes RY(\theta) =
    \begin{pmatrix}
        1 & 0 & 0 & 0 \\
        0 & \cos{\th} & 0 & -\sin{\th} \\
        0 & 0 & 1 & 0 \\
        0 & \sin{\th} & 0 & \cos{\th}
    \end{pmatrix}
```
"""
function CRYRevGate(θ::Number)
    
    QCO._verify_θ_bounds(θ)

    CRYRev = Array{Complex{Float64},2}([ 1 0 0 0            
                                      0 cos(θ/2) 0 -(sin(θ/2))      
                                      0 0 1 0
                                      0 (sin(θ/2)) 0 cos(θ/2)])

    return QCO.round_complex_values(CRYRev)
end

@doc raw"""
    CRZGate(θ::Number)

Two-qubit controlled [RZGate](@ref). 

**Circuit Representation**
```
q_0: ────■────
     ┌───┴───┐
q_1: ┤ RZ(ϴ) ├
     └───────┘
```

**Matrix Representation**

```math
\newcommand{\th}{\frac{\theta}{2}}

CRZ(\theta)\ q_1, q_0 =
|0\rangle\langle0| \otimes I + |1\rangle\langle1| \otimes RZ(\theta) =
    \begin{pmatrix}
        1 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 \\
        0 & 0 & e^{-i\th}  &  0 \\
        0 & 0 & 0 & e^{i\th}
    \end{pmatrix}
```
"""
function CRZGate(θ::Number)
    
    QCO._verify_θ_bounds(θ)

    CRZ = Array{Complex{Float64},2}([ 1 0 0 0            
                                      0 1 0 0       
                                      0 0  (cos(θ/2) - (sin(θ/2))im)  0 
                                      0 0  0  (cos(θ/2) + (sin(θ/2))im)])

    return QCO.round_complex_values(CRZ)
end

@doc raw"""
    CRZRevGate(θ::Number)

Two-qubit controlled reverse [RZGate](@ref). 

**Circuit Representation**
```
     ┌───────┐
q_1: ┤ RZ(ϴ) ├
     └───┬───┘
q_0: ────■────
```

**Matrix Representation**

```math
\newcommand{\th}{\frac{\theta}{2}}

CRZRev(\theta)\ q_1, q_0 =
|0\rangle\langle0| \otimes I + |1\rangle\langle1| \otimes RZ(\theta) =
    \begin{pmatrix}
        1 & 0 & 0 & 0 \\
        0 & e^{-i\th} & 0 & 0 \\
        0 & 0 & 1 &  0 \\
        0 & 0 & 0 & e^{i\th}
    \end{pmatrix}
```
"""
function CRZRevGate(θ::Number)
    
    QCO._verify_θ_bounds(θ)

    CRZRev = Array{Complex{Float64},2}([ 1 0 0 0            
                                      0 (cos(θ/2) - (sin(θ/2))im) 0 0       
                                      0 0 1 0 
                                      0 0  0  (cos(θ/2) + (sin(θ/2))im)])

    return QCO.round_complex_values(CRZRev)
end

@doc raw"""
    CU3Gate(θ::Number, ϕ::Number, λ::Number)

Two-qubit, controlled version of the universal rotation gate with three Euler angles ([U3Gate](@ref)). 

**Circuit Representation**
```
q_0: ──────■──────
     ┌─────┴─────┐
q_1: ┤ U3(ϴ,φ,λ) ├
     └───────────┘
```

**Matrix Representation**

```math
\newcommand{\th}{\frac{\theta}{2}}

CU3(\theta, \phi, \lambda)\ q_1, q_0 =
                |0\rangle\langle 0| \otimes I +
                |1\rangle\langle 1| \otimes U3(\theta,\phi,\lambda) =
                \begin{pmatrix}
                    1 & 0   & 0                  & 0 \\
                    0 & 1   & 0                  & 0 \\
                    0 & 0   & \cos(\th)          & -e^{i\lambda}\sin(\th) \\
                    0 & 0   & e^{i\phi}\sin(\th) & e^{i(\phi+\lambda)}\cos(\th)
                \end{pmatrix}
```
"""
function CU3Gate(θ::Number, ϕ::Number, λ::Number)

    QCO._verify_θ_bounds(θ)
    QCO._verify_ϕ_bounds(ϕ)
    QCO._verify_λ_bounds(λ)

    CU3 = Array{Complex{Float64},2}([ 1 0 0 0            
                                      0 1 0 0       
                                      0 0  cos(θ/2)  -(cos(λ)+(sin(λ))im)*sin(θ/2) 
                                      0 0 (cos(ϕ)+(sin(ϕ))im)*sin(θ/2)  (cos(λ+ϕ)+(sin(λ+ϕ))im)*cos(θ/2)])

    return QCO.round_complex_values(CU3)
end

@doc raw"""
    CU3RevGate(θ::Number, ϕ::Number, λ::Number)

Two-qubit, reverse controlled version of the universal rotation gate with three Euler angles ([U3Gate](@ref)). 

**Circuit Representation**
```
     ┌────────────┐
q_1: ┤  U3(ϴ,φ,λ) ├
     └──────┬─────┘
q_0: ───────■──────
```

**Matrix Representation**

```math
\newcommand{\th}{\frac{\theta}{2}}

CU3(\theta, \phi, \lambda)\ q_1, q_0 =
                |0\rangle\langle 0| \otimes I +
                |1\rangle\langle 1| \otimes U3(\theta,\phi,\lambda) =
                \begin{pmatrix}
                    1 & 0   & 0  & 0 \\
                    0 & \cos(\th)   & 0 & -e^{i\lambda}\sin(\th) \\
                    0 & 0   &  1 & 0 \\
                    0 & e^{i\phi}\sin(\th)  & 0  & e^{i(\phi+\lambda)}\cos(\th)
                \end{pmatrix}
```
"""
function CU3RevGate(θ::Number, ϕ::Number, λ::Number)

    QCO._verify_θ_bounds(θ)
    QCO._verify_ϕ_bounds(ϕ)
    QCO._verify_λ_bounds(λ)

    CU3Rev = Array{Complex{Float64},2}([ 1 0 0 0            
                                      0 cos(θ/2) 0 -(cos(λ)+(sin(λ))im)*sin(θ/2)    
                                      0 0 1 0 
                                      0 (cos(ϕ)+(sin(ϕ))im)*sin(θ/2) 0 (cos(λ+ϕ)+(sin(λ+ϕ))im)*cos(θ/2)])

    return QCO.round_complex_values(CU3Rev)
end

@doc raw"""
    SwapGate()

Two-qubit, symmetric, SWAP gate. 

**Circuit Representation**
```
q_0: ─X─
      │
q_1: ─X─
```

**Matrix Representation**

```math
SWAP = \begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 0 & 1 & 0 \\
0 & 1 & 0 & 0 \\
0 & 0 & 0 & 1
\end{pmatrix}
```
"""
function SwapGate()

    return Array{Complex{Float64},2}([1 0 0 0; 0 0 1 0; 0 1 0 0; 0 0 0 1])

end

@doc raw"""
    iSwapGate()

Two-qubit, symmetric and clifford, iSWAP gate. This is an entangling swapping gate where the qubits 
obtain a phase of ``i`` if the state of the qubits is swapped.

**Circuit Representation**
```
q_0: ─⨂─
      │     
q_1: ─⨂─    
```
Minimum depth representation
```
      ┌───┐     ┌───┐ ┌───┐
q_0: ─┤ X ├──■──┤ S ├─┤ X ├─
      └─┬─┘┌─┴─┐└───┘ └─┬─┘
q_1: ───■──┤ X ├────────■──
           └───┘
```

**Matrix Representation**

```math
iSWAP = \begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 0 & i & 0 \\
0 & i & 0 & 0 \\
0 & 0 & 0 & 1
\end{pmatrix}
```
"""
function iSwapGate()

    return Array{Complex{Float64},2}([1 0 0 0; 0 0 im 0; 0 im 0 0; 0 0 0 1])

end

@doc raw"""
    CSXGate()

Two-qubit controlled version of ([SXGate](@ref)). 

**Circuit Representation**
```
q_0: ─────■─────
     ┌────┴────┐
q_1: ┤ sqrt(X) ├
     └─────────┘
```

**Matrix Representation**

```math
CSXGate = |0 \rangle\langle 0| \otimes I + |1 \rangle\langle 1| \otimes SX = \begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \\
0 & 0 & 0.5+0.5i & 0.5-0.5i \\
0 & 0 & 0.5-0.5i & 0.5+0.5i
    \end{pmatrix}
```
"""
function CSXGate()

    # |0⟩⟨0| ⊗ I
    control_0 = kron(Array{Complex{Float64},2}([1 0; 0 0]) , QCO.IGate(1))
    # |1⟩⟨1| ⊗ SX 
    control_1 = kron(Array{Complex{Float64},2}([0 0; 0 1]) , QCO.SXGate())
        
    return control_0 + control_1 
end

@doc raw"""
    CSXRevGate()

Two-qubit controlled version of the reverse ([SXGate](@ref)). 

**Circuit Representation**
```
     ┌─────────┐
q_1: ┤ sqrt(X) ├
     └────┬────┘
q_0: ─────■────
```

**Matrix Representation**

```math
CSXRevGate = I \otimes |0\rangle\langle 0| + SX \otimes |1\rangle\langle 1| = \begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 0.5+0.5i & 0 & 0.5-0.5i \\
0 & 0 & 1 & 0 \\
0 & 0.5-0.5i & 0 & 0.5+0.5i
    \end{pmatrix}
```
"""
function CSXRevGate()

    # I ⊗ |0⟩⟨0|
    control_0 = kron(QCO.IGate(1), Array{Complex{Float64},2}([1 0; 0 0]))
    # SX ⊗ |1⟩⟨1| 
    control_1 = kron(QCO.SXGate(), Array{Complex{Float64},2}([0 0; 0 1]))
    
    return control_0 + control_1
end


@doc raw"""
    MGate()

Two-qubit Magic gate, also known as the Ising coupling or the XX gate.

Reference: [https://doi.org/10.1103/PhysRevA.69.032315](https://doi.org/10.1103/PhysRevA.69.032315)

**Circuit Representation**
```
      ┌───┐        ┌───┐
q_0: ─┤ X ├────────┤ S ├
      └─┬─┘        └─┬─┘        
        │   ┌───┐  ┌─┴─┐
q_1: ───■───┤ H ├──┤ S ├
            └───┘  └───┘
```

**Matrix Representation**

```math
M = \frac{1}{\sqrt{2}} \begin{pmatrix}
1 & i & 0 & 0 \\
0 & 0 & i & 1 \\
0 & 0 & i & -1 \\
1 & -i & 0 & 0
\end{pmatrix}
```
"""
function MGate()

    return Array{Complex{Float64},2}(1/sqrt(2)*[1 im 0 0; 0 0 im 1; 0 0 im -1; 1 -im 0 0])

end

@doc raw"""
    QFT2Gate()

Two-qubit Quantum Fourier Transform (QFT) gate, where the QFT operation on n-qubits is given by: 
```math
|j\rangle \mapsto \frac{1}{2^{n/2}} \sum_{k=0}^{2^n - 1} e^{2\pi ijk / 2^n} |k\rangle
```

**Circuit Representation**
```
     ┌──────┐
q_0: ┤      ├
     │ QFT2 │   
q_1: ┤      ├ 
     └──────┘ 
```

**Matrix Representation**

```math
M = \frac{1}{2} \begin{pmatrix}
1 & 1 & 1 & 1 \\
1 & i & -1 & -i \\
1 & -1 & 1 & -1 \\
1 & -i & -1 & i
\end{pmatrix}
```
"""
function QFT2Gate()

    return Array{Complex{Float64},2}(0.5*[1 1 1 1; 1 im -1 -im; 1 -1 1 -1; 1 -im -1 im])

end

@doc raw"""
    HCoinGate()

Two-qubit, Hadamard Coin gate when implemented in tune with the quantum cellular automata. 
Reference: [https://doi.org/10.1007/s11128-018-1983-x](https://doi.org/10.1007/s11128-018-1983-x), [https://arxiv.org/pdf/2106.03115.pdf](https://arxiv.org/pdf/2106.03115.pdf)

**Matrix Representation**

```math
HCoinGate = \begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{2}} & 0 \\
0 & \frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{2}} & 0 \\
0 & 0 & 0 & 1
\end{pmatrix}
```
"""
function HCoinGate()

    return Array{Complex{Float64},2}([1 0 0 0; 0 1/sqrt(2) 1/sqrt(2) 0; 0 1/sqrt(2) -1/sqrt(2) 0; 0 0 0 1])

end

@doc raw"""
    GroverDiffusionGate()

Two-qubit, Grover's diffusion operator, a key building block of the Glover's algorithm used to find a specific
item (with probability > 0.5) within a randomly ordered database of N items in O(sqrt(N)) operations. 
Reference: [https://arxiv.org/pdf/1804.03719.pdf](https://arxiv.org/pdf/1804.03719.pdf)

**Matrix Representation**

```math
GroverDiffusionGate = \frac{1}{2}\begin{pmatrix}
1 & -1 & -1 & -1 \\
-1 & 1 & -1 & -1 \\
-1 & -1 & 1 & -1 \\ 
-1 & -1 & -1 & 1
\end{pmatrix}
```
"""
function GroverDiffusionGate()

    return Array{Complex{Float64},2}(0.5*[1 -1 -1 -1; -1 1 -1 -1; -1 -1 1 -1; -1 -1 -1 1]) 

end

@doc raw"""
    SycamoreGate()

Two-qubit Sycamore Gate, native to Google's universal quantum processor.
Reference: [quantumai.google/cirq/google/devices](https://quantumai.google/cirq/google/devices)

**Circuit Representation**
```
     ┌──────┐
q_0: ┤      ├
     │ SYC  │   
q_1: ┤      ├ 
     └──────┘ 
```

**Matrix Representation**

```math

SycamoreGate() = \begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 0 & -i & 0 \\
0 & -i & 0 & 0 \\ 
0 & 0 & 0 & e^{-i \frac{\pi}{6}}
\end{pmatrix}

```
"""
function SycamoreGate()

    return  Array{Complex{Float64},2}([ 1 0 0 0            
                                      0 0 -im 0    
                                      0 -im 0 0
                                      0 0 0 cos(pi/6)-(sin(pi/6))im])

end

#---------------------------------------#
#            Three-qubit gates          #
#---------------------------------------#

@doc raw"""
    ToffoliGate()

Three-qubit Toffoli gate, also known as the CCX (controlled-controlled-NOT) gate. 

**Circuit Representation**
```
q_0: ──■──
       │
q_1: ──■──
     ┌─┴─┐
q_2: ┤ X ├
     └───┘
```

**Matrix Representation**

```math
Toffoli     =
            |0 \rangle \langle 0| \otimes I \otimes I + |1 \rangle \langle 1| \otimes CXGate =
            \begin{pmatrix}
                1 & 0 & 0 & 0 & 0 & 0 & 0 & 0\\
                0 & 1 & 0 & 0 & 0 & 0 & 0 & 0\\
                0 & 0 & 1 & 0 & 0 & 0 & 0 & 0\\
                0 & 0 & 0 & 1 & 0 & 0 & 0 & 0\\
                0 & 0 & 0 & 0 & 1 & 0 & 0 & 0\\
                0 & 0 & 0 & 0 & 0 & 1 & 0 & 0\\
                0 & 0 & 0 & 0 & 0 & 0 & 0 & 1\\
                0 & 0 & 0 & 0 & 0 & 0 & 1 & 0
            \end{pmatrix}
```
"""
function ToffoliGate()

    return Array{Complex{Float64},2}([1  0  0  0  0  0  0  0
                                      0  1  0  0  0  0  0  0
                                      0  0  1  0  0  0  0  0
                                      0  0  0  1  0  0  0  0
                                      0  0  0  0  1  0  0  0
                                      0  0  0  0  0  1  0  0
                                      0  0  0  0  0  0  0  1
                                      0  0  0  0  0  0  1  0])

end

@doc raw"""
    CSwapGate()

Three-qubit, controlled [SwapGate](@ref), also known as the [Fredkin gate](https://en.wikipedia.org/wiki/Fredkin_gate).

**Circuit Representation**
```
q_0: ─■─
      │
q_1: ─X─
      │
q_2: ─X─
```

**Matrix Representation**

```math
CSwapGate =
            \begin{pmatrix}
            1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
            0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
            0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
            0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
            0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
            0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
            0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\
            0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \\
        \end{pmatrix}
```
"""
function CSwapGate()

    return Array{Complex{Float64},2}([1  0  0  0  0  0  0  0
                                      0  1  0  0  0  0  0  0
                                      0  0  1  0  0  0  0  0
                                      0  0  0  1  0  0  0  0
                                      0  0  0  0  1  0  0  0
                                      0  0  0  0  0  0  1  0
                                      0  0  0  0  0  1  0  0
                                      0  0  0  0  0  0  0  1])

end

@doc raw"""
    CCZGate()

Three-qubit controlled-controlled Z gate. 

**Circuit Representation**
```
q_0: ─■─
      │
q_1: ─■─
      │
q_2: ─■─
```

**Matrix Representation**

```math
CCZGate =
            \begin{pmatrix}
            1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
            0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
            0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
            0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
            0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
            0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\
            0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
            0 & 0 & 0 & 0 & 0 & 0 & 0 & -1 \\
        \end{pmatrix}
```
"""
function CCZGate()

    return Array{Complex{Float64},2}([1  0  0  0  0  0  0  0
                                      0  1  0  0  0  0  0  0
                                      0  0  1  0  0  0  0  0
                                      0  0  0  1  0  0  0  0
                                      0  0  0  0  1  0  0  0
                                      0  0  0  0  0  1  0  0
                                      0  0  0  0  0  0  1  0
                                      0  0  0  0  0  0  0  -1])

end

@doc raw"""
    PeresGate()

Three-qubit Peres gate. This gate is equivalent to [ToffoliGate](@ref) followed by the [CNotGate](@ref) in 3 qubits. 
Reference: [https://doi.org/10.1103/PhysRevA.32.3266](https://doi.org/10.1103/PhysRevA.32.3266)

**Circuit Representation**
```
q_0: ──■─────■──          
       │   ┌─┴─┐
q_1: ──■───┤ X ├
     ┌─┴─┐ └───┘
q_2: ┤ X ├──────
     └───┘
```

**Matrix Representation**

```math
PeresGate =
            \begin{pmatrix}
            1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
            0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
            0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
            0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
            0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
            0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \\
            0 & 0 & 0 & 0 & 0 & 1 & 0 & 0 \\
            0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
        \end{pmatrix}
```
"""
function PeresGate()

    return Array{Complex{Float64},2}([1  0  0  0  0  0  0  0
                                      0  1  0  0  0  0  0  0
                                      0  0  1  0  0  0  0  0
                                      0  0  0  1  0  0  0  0
                                      0  0  0  0  0  0  1  0
                                      0  0  0  0  0  0  0  1
                                      0  0  0  0  0  1  0  0
                                      0  0  0  0  1  0  0  0])

end

@doc raw"""
    RCCXGate()

Three-qubit relative (or simplified) Toffoli gate, or the CCX gate. This gate is equivalent to [ToffoliGate](@ref) upto relative phases. 
The advantage of this gate is that it's implementation requires only three CNot (or CX) gates. 
Reference: [https://arxiv.org/pdf/1508.03273.pdf](https://arxiv.org/pdf/1508.03273.pdf)

**Matrix Representation**

```math
RCCXGate =
            \begin{pmatrix}
            1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
            0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
            0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
            0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
            0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
            0 & 0 & 0 & 0 & 0 & -1 & 0 & 0 \\
            0 & 0 & 0 & 0 & 0 & 0 & 0 & -i \\
            0 & 0 & 0 & 0 & 0 & 0 & i & 0 \\
        \end{pmatrix}
```
"""
function RCCXGate()

    return Array{Complex{Float64},2}([1  0  0  0  0  0  0  0
                                      0  1  0  0  0  0  0  0
                                      0  0  1  0  0  0  0  0
                                      0  0  0  1  0  0  0  0
                                      0  0  0  0  1  0  0  0
                                      0  0  0  0  0  -1  0  0
                                      0  0  0  0  0  0  0  -im
                                      0  0  0  0  0  0  im  0])
end

@doc raw"""
    MargolusGate()

Three-qubit Margolus gate, which is a simplified [ToffoliGate](@ref) and coincides with the Toffoli gate up 
to a single change of sign. The advantage of this gate is that its implementation requires only three CNot (or CX) gates. 
Reference: [https://arxiv.org/pdf/quant-ph/0312225.pdf](https://arxiv.org/pdf/quant-ph/0312225.pdf)

**Matrix Representation**

```math
MargolusGate =
            \begin{pmatrix}
            1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
            0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
            0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
            0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
            0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
            0 & 0 & 0 & 0 & 0 & -1 & 0 & 0 \\
            0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \\
            0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
        \end{pmatrix}
```
"""
function MargolusGate()

    return Array{Complex{Float64},2}([1  0  0  0  0  0  0  0
                                      0  1  0  0  0  0  0  0
                                      0  0  1  0  0  0  0  0
                                      0  0  0  1  0  0  0  0
                                      0  0  0  0  1  0  0  0
                                      0  0  0  0  0 -1  0  0
                                      0  0  0  0  0  0  0  1
                                      0  0  0  0  0  0  1  0])
end

@doc raw"""
    CiSwapGate()

Three-qubit controlled version of the [iSwapGate](@ref).
Reference: [https://doi.org/10.1103/PhysRevResearch.2.033097](https://doi.org/10.1103/PhysRevResearch.2.033097)

**Circuit Representation**
```
q_0: ─────■─────
          │
      ┌───────┐
q_1: ─┤       ├─
      │ iSwap │   
q_2: ─┤       ├─ 
      └───────┘ 
```

**Matrix Representation**

```math
CiSwapGate =
            \begin{pmatrix}
            1 & 0 & 0 & 0 & 0 & 0 & 0 & 0 \\
            0 & 1 & 0 & 0 & 0 & 0 & 0 & 0 \\
            0 & 0 & 1 & 0 & 0 & 0 & 0 & 0 \\
            0 & 0 & 0 & 1 & 0 & 0 & 0 & 0 \\
            0 & 0 & 0 & 0 & 1 & 0 & 0 & 0 \\
            0 & 0 & 0 & 0 & 0 & 0 & i & 0 \\
            0 & 0 & 0 & 0 & 0 & i & 0 & 0 \\
            0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \\
        \end{pmatrix}
```
"""
function CiSwapGate()

    return Array{Complex{Float64},2}([1  0  0  0  0  0  0  0
                                      0  1  0  0  0  0  0  0
                                      0  0  1  0  0  0  0  0
                                      0  0  0  1  0  0  0  0
                                      0  0  0  0  1  0  0  0
                                      0  0  0  0  0  0 im  0
                                      0  0  0  0  0 im  0  0
                                      0  0  0  0  0  0  0  1])
end