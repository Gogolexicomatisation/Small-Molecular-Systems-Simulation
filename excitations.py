import numpy as np

from qiskit import QuantumCircuit
from qiskit.circuit import Parameter

def fermionic_single_excitation(qubits_involved, pauli_involved, qubits, param_name= None):

    if param_name is None:
        theta = Parameter('theta')
    else:
        theta = Parameter(param_name)
    
    string_exci = ''
    for i in range(len(pauli_involved)):
        string_exci += pauli_involved[i]

    circuit = QuantumCircuit(qubits)

    for i in range(3):
        if pauli_involved[i] == 'X':
            circuit.h(qubits_involved[i])
        else:
            circuit.rx(-np.pi/2, qubits_involved[i])
    
    for i in reversed(range(0, 2)): 
        circuit.cx(i+1, i)

    if string_exci =='XY':
        circuit.rz((1/2)*theta, 0)
    else:
        circuit.rz(-2*theta, 0)
    
    for i in range(qubits_involved[0], 1):
        circuit.cx(i+1, i)
    
    for i in range(3):
        if pauli_involved[i] == 'X':
            circuit.x(qubits_involved[i])
        else:
            circuit.ry(+np.pi/2, qubits_involved[i])
    
    return circuit

def fermionic_double_excitation(qubits_involved, pauli_involved, qubits, param_name = None):

    if param_name is None:
        theta = Parameter('theta')
    else:
        theta = Parameter(param_name)

    string_exci = ''
    for i in range(len(pauli_involved)):
        string_exci += pauli_involved[i]
    
    circuit = QuantumCircuit(qubits)

    for i in range(4):
        if pauli_involved[i] == 'X':
            circuit.h(i)
        else:
            circuit.ry(np.pi/2, i)
    
    for i in range(4):
        circuit.cx(i,i+1)

    if string_exci == 'XYXX' or string_exci == 'YXXX' or string_exci == 'YYYX' or string_exci == 'YYXY':
        circuit.rz(-8*theta, 0)
    else:
        circuit.rz(+(1/8)*theta, 0)

    for i in range(4):
        circuit.cx(i, i+1)
    
    for i in range(len(pauli_involved)):
        if pauli_involved[i] == 'X':
            circuit.z(i)
        else:
            circuit.rx(-np.pi/2, i)
    
    return circuit
 