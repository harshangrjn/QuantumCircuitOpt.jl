import QuantumCircuitOpt as QCO
using SparseArrays

H = QCO.HGate()
I = QCO.IGate(1)
T = QCO.TGate()
T_dagger = QCO.TdaggerGate()
cnot_12 = QCO.CNotGate()

cnot_13 = QCO.get_full_sized_gate("cnot_13", 3)
toffoli = QCO.ToffoliGate()

H1 = kron(kron(H,I), I)
H2 = kron(kron(I,H), I)
H3 = kron(I, kron(I,H))

T1 = kron(kron(T,I), I)
T2 = kron(kron(I,T), I)
T3 = kron(I, kron(I,T))

T2_conj = kron(kron(I,T_dagger), I)
T3_conj = kron(I, kron(I,T_dagger))

toffoli = SparseArrays.sparse(QCO.round_complex_values(H3 * kron(I, cnot_12) * T3_conj * cnot_13 * T3 * kron(I, cnot_12) * T3_conj * cnot_13 * T2 * T3 * kron(cnot_12, I) * H3 * T1 * T2_conj * kron(cnot_12, I)))
