#----------------------------------------#
#            Single-qubit gates          #
#----------------------------------------#

@doc raw"""
    IGate(num_qubits::Int64)

Identity gate corresponds to the regular identity matrix for a given number of qubits.

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

    if !(-π <= θ <= π)
        Memento.error(_LOGGER, "θ angle in U3Gate is not within valid bounds")
    end
    if !(-2*π <= ϕ <= 2*π)
        Memento.error(_LOGGER, "ϕ angle in U3Gate is not within valid bounds")
    end
    if !(-2*π <= λ <= 2*π)
        Memento.error(_LOGGER, "λ angle in U3Gate is not within valid bounds")
    end

    U3 = Array{Complex{Float64},2}([           cos(θ/2)               -(cos(λ) + (sin(λ))im)*sin(θ/2) 
                                    (cos(ϕ) + (sin(ϕ))im)*sin(θ/2)  (cos(λ+ϕ) + (sin(λ+ϕ))im)*cos(θ/2)])

    return round_complex_values(U3)
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

    if !(-2*π <= ϕ <= 2*π)
        Memento.error(_LOGGER, "ϕ angle in U2Gate is not within valid bounds")
    end
    if !(-2*π <= λ <= 2*π)
        Memento.error(_LOGGER, "λ angle in U2Gate is not within valid bounds")
    end
    
    θ = π/2

    U2 = U3Gate(θ, ϕ, λ)
    
    return U2
end

@doc raw"""
    U1Gate(λ::Number)

Universal single-qubit rotation gate with one Euler angle, ``\lambda``. U1Gate represents rotation about the Z axis and 
is the special case of [U3Gate](@ref). Also note that ``U1(\pi) = ``[ZGate](@ref), ``U1(\pi/2) = ``[SGate](@ref) and 
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

    if !(-2*π <= λ <= 2*π)
        Memento.error(_LOGGER, "λ angle in U1Gate is not within valid bounds")
    end
    
    θ = 0
    ϕ = 0

    U1 = U3Gate(θ, ϕ, λ)
    
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

    if !(-2*π <= θ <= 2*π)
        Memento.error(_LOGGER, "θ angle in RXGate is not within valid bounds")
    end

    RX = Array{Complex{Float64},2}([cos(θ/2) -(sin(θ/2))im; -(sin(θ/2))im cos(θ/2)])
    
    return round_complex_values(RX)
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

    if !(-2*π <= θ <= 2*π)
        Memento.error(_LOGGER, "θ angle in RYGate is not within valid bounds")
    end

    RY = Array{Complex{Float64},2}([cos(θ/2) -(sin(θ/2)); (sin(θ/2)) cos(θ/2)])
    
    return round_complex_values(RY)
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

    if !(-2*π <= θ <= 2*π)
        Memento.error(_LOGGER, "θ angle in RZGate is not within valid bounds")
    end

    RZ = Array{Complex{Float64},2}([(cos(θ/2) - (sin(θ/2))im) 0; 0 (cos(θ/2) + (sin(θ/2))im)])
    
    return round_complex_values(RZ)
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
    #input angles in radians
    θ = 0
    ϕ = 0

    if !(-2*π <= λ <= 2*π)
        Memento.error(_LOGGER, "λ angle in Phase gate is not within valid bounds")
    end

    return U3Gate(θ, ϕ, λ)

end

#-------------------------------------#
#            Two-qubit gates          #
#-------------------------------------#
@doc raw"""
    NegIGate()

Two-qubit negative Identity gate. 

**Circuit Representation**
```  ┌───────┐
q_1: ┤ RZ(ϴ) ├
     └───────┘
q_0: ─────────
```

**Matrix Representation**

```math
NegI = \begin{pmatrix}
    -1 & 0 & 0 & 0 \\
    0 & -1 & 0 & 0 \\
    0 & 0 & -1 & 0 \\
    0 & 0 & 0 & -1
    \end{pmatrix}
```
"""
function NegIGate()

    return Array{Complex{Float64},2}([-1 0 0 0; 0 -1 0 0; 0 0 -1 0; 0 0 0 -1]) 

end

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
DCNot = \begin{pmatrix}
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
CY = \begin{pmatrix}
    1 & 0 & 0 & 0 \\
    0 & 1 & 0 & 0 \\
    0 & 0 & 0 & -i \\
    0 & 0 & i & 0
    \end{pmatrix}
```
"""
function CYGate()

    return Array{Complex{Float64},2}([1 0 0 0; 0 1 0 0; 0 0 0 -im; 0 0 im 0]) 

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
CZ = \begin{pmatrix}
    1 & 0 & 0 & 0 \\
    0 & 1 & 0 & 0 \\
    0 & 0 & 1 & 0 \\
    0 & 0 & 0 & -1
    \end{pmatrix}
```
"""
function CZGate()

    return Array{Complex{Float64},2}([1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 -1]) 

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
CH = \begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \\
0 & 0 & \frac{1}{\sqrt{2}} & \frac{1}{\sqrt{2}} \\
0 & 0 & \frac{1}{\sqrt{2}} & -\frac{1}{\sqrt{2}}
\end{pmatrix}
```
"""
function CHGate()

    return Array{Complex{Float64},2}([1 0 0 0; 0 1 0 0; 0 0 1/sqrt(2) 1/sqrt(2); 0  0 1/sqrt(2) -1/sqrt(2)])

end

@doc raw"""
    CVGate()

Two-qubit, controlled-V gate, which is also the same as Controlled square-root of X gate ([C2SXGate](@ref)).  

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

    return QCO.C2SXGate() 

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
CV = \begin{pmatrix}
        1 & 0 & 0 & 0 \\
        0 & 0.5+0.5i & 0 & 0.5-0.5i \\
        0 & 0 & 1 & 0 \\
        0 & 0.5-0.5i & 0 & 0.5+0.5i
    \end{pmatrix}
```
"""
function CVRevGate()

    return return Array{Complex{Float64},2}([1 0 0 0; 0 0.5+0.5im 0 0.5-0.5im; 0 0 1 0; 0 0.5-0.5im 0 0.5+0.5im])

end

@doc raw"""
    CVdaggerGate()

Two-qubit hermitian conjugate of controlled-V gate, which is also the same as hermitian conjugate Controlled square-root of X gate ([C2SXGate](@ref)).  

**Circuit Representation**
```
q_0: ──■──     
     ┌─┴─┐    
q_1: ┤ V'├     
     └───┘
```

**Matrix Representation**

```math
CV = \begin{pmatrix}
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
    CVRevdaggerGate()

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
CV = \begin{pmatrix}
        1 & 0 & 0 & 0 \\
        0 & 0.5-0.5i & 0 & 0.5+0.5i \\
        0 & 0 & 1 & 0 \\
        0 & 0.5+0.5i & 0 & 0.5-0.5i
    \end{pmatrix}
```
"""
function CVRevdaggerGate()

    return return Array{Complex{Float64},2}([1 0 0 0; 0 0.5-0.5im 0 0.5+0.5im; 0 0 1 0; 0 0.5+0.5im 0 0.5-0.5im])

end

@doc raw"""
    WGate()

Two-qubit, W hermitian gate, typically useful to diagonlize the ([SwapGate](@ref)).  

**Matrix Representation**

```math
CV = \begin{pmatrix}
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

    return round_complex_values(CRX)
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
    
    if !(-2*π <= θ <= 2*π)
        Memento.error(_LOGGER, "θ angle in CRYGate is not within valid bounds")
    end

    CRY = Array{Complex{Float64},2}([ 1 0 0 0            
                                      0 1 0 0       
                                      0 0  cos(θ/2)  -(sin(θ/2)) 
                                      0 0  (sin(θ/2))  cos(θ/2)])

    return round_complex_values(CRY)
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
    
    if !(-2*π <= θ <= 2*π)
        Memento.error(_LOGGER, "θ angle in CRZGate is not within valid bounds")
    end

    CRZ = Array{Complex{Float64},2}([ 1 0 0 0            
                                      0 1 0 0       
                                      0 0  (cos(θ/2) - (sin(θ/2))im)  0 
                                      0 0  0  (cos(θ/2) + (sin(θ/2))im)])

    return round_complex_values(CRZ)
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

    if !(-π <= θ <= π)
        Memento.error(_LOGGER, "θ angle in CU3Gate is not within valid bounds")
    end
    if !(-2*π <= ϕ <= 2*π)
        Memento.error(_LOGGER, "ϕ angle in CU3Gate is not within valid bounds")
    end
    if !(-2*π <= λ <= 2*π)
        Memento.error(_LOGGER, "λ angle in CU3Gate is not within valid bounds")
    end

    CU3 = Array{Complex{Float64},2}([ 1 0 0 0            
                                      0 1 0 0       
                                      0 0  cos(θ/2)  -(cos(λ)+(sin(λ))im)*sin(θ/2) 
                                      0 0 (cos(ϕ)+(sin(ϕ))im)*sin(θ/2)  (cos(λ+ϕ)+(sin(λ+ϕ))im)*cos(θ/2)])

    return round_complex_values(CU3)
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

Two-qubit, symmetric and clifford, iSWAP gate.

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
    C2SXGate()

Two-qubit controlled ([SXGate](@ref)), which is also the square root of the Pauli-[XGate](@ref). 

**Circuit Representation**
```
q_0: ─────■─────
     ┌────┴────┐
q_1: ┤ sqrt(X) ├
     └─────────┘
```

**Matrix Representation**

```math
C2SXGate = \begin{pmatrix}
1 & 0 & 0 & 0 \\
0 & 1 & 0 & 0 \\
0 & 0 & 0.5+0.5i & 0.5-0.5i \\
0 & 0 & 0.5-0.5i & 0.5+0.5i
    \end{pmatrix}
```
"""
function C2SXGate()

    return Array{Complex{Float64},2}([1 0 0 0; 0 1 0 0; 0 0 0.5+0.5im 0.5-0.5im; 0 0 0.5-0.5im 0.5+0.5im]) 

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
            0 & 0 & 0 & 0 & 0 & 0 & 0 & 1 \\
            0 & 0 & 0 & 0 & 0 & 0 & 1 & 0 \\
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
                                      0  0  0  0  0  0  0  1
                                      0  0  0  0  0  0  1  0
                                      0  0  0  0  0  1  0  0
                                      0  0  0  0  1  0  0  0])

end
