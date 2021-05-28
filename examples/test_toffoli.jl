using QuantumCircuitOpt
using SparseArrays

const QCO = QuantumCircuitOpt

H = QCO.get_elementary_gates(2)["hadamard_H"]
I = QCO.get_elementary_gates(2)["I_2"]
T = QCO.get_elementary_gates(2)["ph_shift_T"]
T_conjugate = QCO.get_elementary_gates(2)["ph_shift_T_conj"]
cnot_12 = QCO.get_elementary_gates(2)["cnot_12"]

cnot_13 = QCO.get_elementary_gates(3)["cnot_13"]
toffoli = QCO.get_elementary_gates(3)["toffoli"]

H1 = kron(kron(H,I), I)
H2 = kron(kron(I,H), I)
H3 = kron(I, kron(I,H))

T1 = kron(kron(T,I), I)
T2 = kron(kron(I,T), I)
T3 = kron(I, kron(I,T))

T2_conj = kron(kron(I,T_conjugate), I)
T3_conj = kron(I, kron(I,T_conjugate))

toffoli = SparseArrays.sparse(QCO.round_complex_values(H3 * kron(I, cnot_12) * T3_conj * cnot_13 * T3 * kron(I, cnot_12) * T3_conj * cnot_13 * T2 * T3 * kron(cnot_12, I) * H3 * T1 * T2_conj * kron(cnot_12, I)))
