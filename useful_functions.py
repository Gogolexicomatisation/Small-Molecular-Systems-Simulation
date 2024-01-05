import numpy as np
from copy import copy

from qiskit.utils import QuantumInstance
from qiskit.opflow.state_fns import StateFn, CircuitStateFn
from qiskit.primitives import Estimator
from qiskit.opflow.converters import CircuitSampler
from qiskit.quantum_info import PauliList, Pauli, SparsePauliOp
from qiskit.opflow import PauliOp

from openfermion.transforms import jordan_wigner
from openfermion.utils import count_qubits
from openfermion.circuits import uccsd_singlet_generator, uccsd_singlet_get_packed_amplitudes

def of_to_qiskit_op(qubit_operator):
    '''
    Function that performs the transformation of an openfermion qubit operator into a Qiskit Pauli operator
    return :
        op is a Qiskit Pauli operator, usually representing a molecular hamiltonian in the Qiskit-qubit formalism
    '''
    op = 0

    for qubit_terms, qubit_coeff in qubit_operator.terms.items():
        string_term = "I"*count_qubits(qubit_operator)
        for i, (term_qubit, term_pauli) in enumerate(qubit_terms):
            string_term = (string_term[:term_qubit] + term_pauli + string_term[term_qubit + 1 :])
        
        op += SparsePauliOp(Pauli(string_term), coeffs=qubit_coeff)

    return op

def measure_qc(circuit, backend, nshots=-1): #Never used

    job = backend.run(circuit, nshots)

    result = job.result()
    counts = result.get_counts()

    return counts

def expectation_value(pauli_operator, circuit, params):

    estimator = Estimator()
    result = estimator.run(circuit, pauli_operator, params).result()

    mean_value = result.values[0].real

    return mean_value

#___________________________________________ SDK functions _____________________________________

def create_list_p_operators(identity, p):

    list_p_operators = []

    for k in range(len(identity)):
        new_op = copy(identity)
        new_op = new_op[:k] + p + new_op[k+1:] #TODO old code: new_op[k+:1]
        list_p_operators.append(new_op)

    list_p_operators.reverse()

    return list_p_operators

def create_pauli_strings(nqubits):

    identity = ''
    for i in range(nqubits):
        identity = identity + 'I'
    
    list_x_operators = create_list_p_operators(identity, 'X')
    list_y_operators = create_list_p_operators(identity, 'Y')
    list_z_operators = create_list_p_operators(identity, 'Z')

    return list_x_operators, list_y_operators, list_z_operators

def create_pauli_lists(nqubits):

    xs, ys, zs = create_pauli_strings(nqubits)

    return PauliList(xs), PauliList(ys), PauliList(zs)

def UCCSD_excitations_generator(molecule):
    '''
    Generates the UCCSD excitations for a given molecule with Openfermion. 
    return :
        uccsd_jw is an Openfermion qubit operator, representing all the UCCSD excitations of the given molecule in the qubit formalism, obtained through a Jordan-Wigner
        transformation
    '''
    electrons = molecule.n_electrons
    qubits = molecule.n_qubits

    ccsd_single_amps = molecule.ccsd_single_amps
    ccsd_double_amps = molecule.ccsd_double_amps

    singlet_amps = uccsd_singlet_get_packed_amplitudes(ccsd_single_amps, ccsd_double_amps, qubits, electrons)
    uccsd_singlet = uccsd_singlet_generator(singlet_amps, qubits, electrons)

    uccsd_jw = jordan_wigner(uccsd_singlet)

    return uccsd_jw