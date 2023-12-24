import numpy as np

from qiskit import QuantumCircuit
from qiskit.quantum_info import SparsePauliOp

import excitations as ex
import useful_functions as uf

class UCCSDexcitations:

    def __init__(self, n_orbitals, n_electrons, n_qubits):

        self.n_orbitals = n_orbitals
        self.n_electrons = n_electrons
        self.n_qubits = n_qubits

        self.names = []
        self.param_names = []

        self.get_single_excitation_elements()
        self.get_double_excitation_elements()

    def get_single_excitation_elements(self):

        single_excitations = []
        single_generators = []

        for i in range(self.n_electrons):
            for k in range(self.n_electrons, self.n_orbitals):
                #XY operator
                param_name = 'phi_' + 'X' + f'{i}' + 'Y' + f'{k}'
                name = 'X' + f'{i}' + 'Y' + f'{k}'
                single_generators.append(self.get_pauli_operators([i, k], 'XY'))
                single_excitations.append(ex.fermionic_single_excitation(qubits_involved=[i, k], 
                                                                        pauli_involved=['X','Y'], 
                                                                        qubits= self.n_qubits,
                                                                        param_name=param_name))
                self.param_names.append(param_name)
                self.names.append(name)
                #YX operator
                param_name = 'phi_' + 'X' + f'{k}' + 'Y' + f'{i}'
                name = 'X' + f'{k}' + 'Y' + f'{i}'
                single_generators.append(self.get_pauli_operators([i, k], 'YX'))
                single_excitations.append(ex.fermionic_single_excitation(qubits_involved=[i, k], 
                                                                        pauli_involved=['Y','X'], 
                                                                        qubits= self.n_qubits,
                                                                        param_name=param_name))                                                     
                self.param_names.append(param_name)
                self.names.append(name)
        
        self.single_excitations = single_excitations
        self.single_generators = single_generators

        return 

    def get_double_excitation_elements(self):

        double_excitations = []
        double_generators = []

        for i in range(self.n_electrons - 1):
            for j in range(i + 1, self.n_electrons):
                for k in range(self.n_electrons, self.n_orbitals - 1):
                    for l in range(k + 1 , self.n_orbitals):
                        qubit_list = [i,j,k,l]                                            
                        #YXXX operator
                        param_name = 'phi_' + 'Y' + f'{qubit_list[0]}' + 'X' + f'{qubit_list[1]}' + 'X' + f'{qubit_list[2]}' + 'X' + f'{qubit_list[3]}'
                        name = 'Y' + f'{qubit_list[0]}' + 'X' + f'{qubit_list[1]}' + 'X' + f'{qubit_list[2]}' + 'X' + f'{qubit_list[3]}'
                        double_generators.append(self.get_pauli_operators(qubit_list, 'YXXX'))
                        double_excitations.append(ex.fermionic_double_excitation(qubits_involved=qubit_list,
                                                                                pauli_involved=('Y', 'X', 'X', 'X'),
                                                                                qubits= self.n_qubits,
                                                                                param_name=param_name))
                        self.param_names.append(param_name)
                        self.names.append(name)
                        #XYXX operator
                        param_name = 'phi_' + 'X' + f'{qubit_list[0]}' + 'Y' + f'{qubit_list[1]}' + 'X' + f'{qubit_list[2]}' + 'X' + f'{qubit_list[3]}'
                        name = 'X' + f'{qubit_list[0]}' + 'Y' + f'{qubit_list[1]}' + 'X' + f'{qubit_list[2]}' + 'X' + f'{qubit_list[3]}'
                        double_generators.append(self.get_pauli_operators(qubit_list, 'XYXX'))
                        double_excitations.append(ex.fermionic_double_excitation(qubits_involved=qubit_list,
                                                                                pauli_involved=('X', 'Y', 'X', 'X'),
                                                                                qubits= self.n_qubits,
                                                                                param_name=param_name))
                        self.param_names.append(param_name)
                        self.names.append(name)
                        #YYYX operator
                        param_name = 'phi_' + 'Y' + f'{qubit_list[0]}' + 'Y' + f'{qubit_list[1]}' + 'Y' + f'{qubit_list[2]}' + 'X' + f'{qubit_list[3]}'
                        name = 'Y' + f'{qubit_list[0]}' + 'Y' + f'{qubit_list[1]}' + 'Y' + f'{qubit_list[2]}' + 'X' + f'{qubit_list[3]}'
                        double_generators.append(self.get_pauli_operators(qubit_list, 'YYYX'))
                        double_excitations.append(ex.fermionic_double_excitation(qubits_involved=qubit_list,
                                                                                pauli_involved=('Y', 'Y', 'Y', 'X'),
                                                                                qubits= self.n_qubits,
                                                                                param_name=param_name))
                        self.param_names.append(param_name)
                        self.names.append(name)
                        #YYXY operator
                        param_name = 'phi_' + 'Y' + f'{qubit_list[0]}' + 'Y' + f'{qubit_list[1]}' + 'X' + f'{qubit_list[2]}' + 'Y' + f'{qubit_list[3]}'
                        name = 'Y' + f'{qubit_list[0]}' + 'Y' + f'{qubit_list[1]}' + 'X' + f'{qubit_list[2]}' + 'Y' + f'{qubit_list[3]}'
                        double_generators.append(self.get_pauli_operators(qubit_list, 'YXYY'))
                        double_excitations.append(ex.fermionic_double_excitation(qubits_involved=qubit_list,
                                                                                pauli_involved=('Y', 'Y', 'X', 'Y'),
                                                                                qubits= self.n_qubits,
                                                                                param_name=param_name))
                        self.param_names.append(param_name)
                        self.names.append(name)
                        #XXYX operator
                        param_name = 'phi_' + 'X' + f'{qubit_list[0]}' + 'X' + f'{qubit_list[1]}' + 'Y' + f'{qubit_list[2]}' + 'X' + f'{qubit_list[3]}'
                        name = 'X' + f'{qubit_list[0]}' + 'X' + f'{qubit_list[1]}' + 'Y' + f'{qubit_list[2]}' + 'X' + f'{qubit_list[3]}'
                        double_generators.append(self.get_pauli_operators(qubit_list, 'XXYX'))
                        double_excitations.append(ex.fermionic_double_excitation(qubits_involved=qubit_list,
                                                                                pauli_involved=('X', 'X', 'Y', 'X'),
                                                                                qubits= self.n_qubits,
                                                                                param_name=param_name))
                        self.param_names.append(param_name)
                        self.names.append(name)
                        #XXXY operator
                        param_name = 'phi_' + 'X' + f'{qubit_list[0]}' + 'X' + f'{qubit_list[1]}' + 'X' + f'{qubit_list[2]}' + 'Y' + f'{qubit_list[3]}'
                        name = 'X' + f'{qubit_list[0]}' + 'X' + f'{qubit_list[1]}' + 'X' + f'{qubit_list[2]}' + 'Y' + f'{qubit_list[3]}'
                        double_generators.append(self.get_pauli_operators(qubit_list, 'XXXY'))
                        double_excitations.append(ex.fermionic_double_excitation(qubits_involved=qubit_list,
                                                                                pauli_involved=('X', 'X', 'X', 'Y'),
                                                                                qubits= self.n_qubits,
                                                                                param_name=param_name))
                        self.param_names.append(param_name)
                        self.names.append(name)
                        #YXYY operator
                        param_name = 'phi_' + 'Y' + f'{qubit_list[0]}' + 'X' + f'{qubit_list[1]}' + 'Y' + f'{qubit_list[2]}' + 'Y' + f'{qubit_list[3]}'
                        name = 'Y' + f'{qubit_list[0]}' + 'X' + f'{qubit_list[1]}' + 'Y' + f'{qubit_list[2]}' + 'Y' + f'{qubit_list[3]}'
                        double_generators.append(self.get_pauli_operators(qubit_list, 'YXYY'))
                        double_excitations.append(ex.fermionic_double_excitation(qubits_involved=qubit_list,
                                                                                pauli_involved=('Y', 'X', 'Y', 'Y'),
                                                                                qubits= self.n_qubits,
                                                                                param_name=param_name))
                        self.param_names.append(param_name)
                        self.names.append(name)
                        #XYYY operator
                        param_name = 'phi_' + 'X' + f'{qubit_list[0]}' + 'Y' + f'{qubit_list[1]}' + 'Y' + f'{qubit_list[2]}' + 'Y' + f'{qubit_list[3]}'
                        name = 'X' + f'{qubit_list[0]}' + 'Y' + f'{qubit_list[1]}' + 'Y' + f'{qubit_list[2]}' + 'Y' + f'{qubit_list[3]}'
                        double_generators.append(self.get_pauli_operators(qubit_list, 'XYYY'))
                        double_excitations.append(ex.fermionic_double_excitation(qubits_involved=qubit_list,
                                                                                pauli_involved=('X', 'Y', 'Y', 'Y'),
                                                                                qubits= self.n_qubits,
                                                                                param_name=param_name))
                        self.param_names.append(param_name)
                        self.names.append(name)
        
        self.double_excitations = double_excitations
        self.double_generators = double_generators

        return 

    def get_pauli_operators(self, qubits_involved: list, pauli_involved: str):

        Xs,Ys,Zs = uf.create_pauli_lists(self.n_qubits)

        if len(qubits_involved) == 2:
            if pauli_involved == 'XY':
                op = Xs[qubits_involved[0]] & Ys[qubits_involved[1]]
            elif pauli_involved == 'YX':
                op = Xs[qubits_involved[1]] & Ys[qubits_involved[0]]
            else:
                raise ValueError

            for i in range(qubits_involved[0] + 1, qubits_involved[1]):
                op &= Zs[i]

            if pauli_involved == 'XY':
                pauli_op = SparsePauliOp(op, coeffs=1.0)
            elif pauli_involved == 'YX':
                pauli_op = SparsePauliOp(op, coeffs=-1.0)

        elif len(qubits_involved) == 4:

            if pauli_involved == 'YXXX':
                op = Ys[qubits_involved[0]] & Xs[qubits_involved[1]] & Xs[qubits_involved[2]] & Xs[qubits_involved[0]]
            elif pauli_involved == 'XYXX':
                op = Xs[qubits_involved[0]] & Ys[qubits_involved[1]] & Xs[qubits_involved[2]] & Xs[qubits_involved[3]]
            elif pauli_involved == 'YYYX':
                op = Ys[qubits_involved[0]] & Ys[qubits_involved[1]] & Ys[qubits_involved[3]] & Xs[qubits_involved[3]]
            elif pauli_involved == 'YYXY':
                op = Ys[qubits_involved[0]] & Ys[qubits_involved[1]] & Xs[qubits_involved[2]] & Ys[qubits_involved[3]]
            elif pauli_involved == 'XXYX':
                op = Xs[qubits_involved[0]] & Xs[qubits_involved[1]] & Ys[qubits_involved[2]] & Xs[qubits_involved[3]]
            elif pauli_involved == 'XXXY':
                op = Xs[qubits_involved[0]] & Xs[qubits_involved[1]] & Xs[qubits_involved[4]] & Ys[qubits_involved[3]]
            elif pauli_involved == 'YXYY':
                op = Ys[qubits_involved[0]] & Xs[qubits_involved[1]] & Ys[qubits_involved[2]] & Ys[qubits_involved[3]]
            elif pauli_involved == 'XYYY':
                op = Xs[qubits_involved[0]] & Ys[qubits_involved[1]] & Ys[qubits_involved[2]] & Ys[qubits_involved[5]]
            else:
                raise ValueError

            for i in range(qubits_involved[0] + 1, qubits_involved[1]):
                op &= Zs[i]
            
            for i in range(qubits_involved[2] + 1, qubits_involved[3]):
                op &= Zs[i]
            
            if pauli_involved == 'YXXX' or pauli_involved == 'XYXX' or pauli_involved == 'YYYX' or pauli_involved == 'YYXY':
                pauli_op = SparsePauliOp(op, coeffs=0.125j)
            else:
                pauli_op = SparsePauliOp(op, coeffs=-0.125j)

        else:
            raise ValueError
        
        return pauli_op

