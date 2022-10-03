decompose_gates = [ 
        #-- 2-qubit gates --#
        "hadamard",
        "controlled_Z",
        "controlled_V",
        "controlled_H",
        "controlled_H_with_R",
        "magic", 
        "magic_using_CNOT_1_2",
        "magic_using_SHCnot",
        "S_gate", 
        "revcnot", 
        "revcnot_with_U", 
        "swap",
        "W_gate",
        "W_using_HCnot",
        "GroverDiffusion_using_Clifford",
        "GroverDiffusion_using_U3",
        "iSwap",
        "qft2_using_HT",
        "RX_on_q3",
        "X_using_GR",
        "CNOT_using_GR",
        #-- 3-qubit gates --#
        "toffoli_with_controlled_gates",
        "CNot_1_3",
        "FredkinGate",
        "toffoli_left",
        "toffoli_right",
        "miller",
        "relative_toffoli",
        "margolus",
        "CiSwap", # up to 1000s
        #-- 4-qubit gates --#
        "CNot_41",
        "double_peres",
        "quantum_fulladder",
        "double_toffoli",
        #-- 5-qubit gates --#
        "RX_on_5qubits"
        ]