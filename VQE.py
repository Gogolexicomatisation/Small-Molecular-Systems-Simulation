import numpy as np
from copy import copy

from scipy.optimize import minimize

import molecule
import ansatze
import useful_functions as uf

class MoleculeSimulator:

    def __init__(self, mol:molecule.Molecule, nshots, backend):

        mol_attributes = [attr for attr in dir(mol) if '__' not in attr]
        for attr in mol_attributes:
            setattr(self, attr, mol.__dict__[attr])
        
        #___________________________________________ SDK-Options _____________________________________
        
        self.current_state = self.initial_state
        self.nshots = min([nshots,1])

        self.backend = backend

        self.qiskit_ham = uf.of_to_qiskit_op(self.hamiltonian)

        #___________________________________________ Initialize _____________________________________

        self.optimizer = 'cobyla'
        self.energy_threshold = 1e-10

        self.variational_energies = {}
        self.optimal_parameters = []

        self.list_of_operators = {}
        self.qubits_acted_on = []

        self.pool = ansatze.UCCSDexcitations(self.n_orbitals, self.n_electrons, self.n_qubits)
        self.generators = self.pool.single_generators + self.pool.double_generators
        self.excitations = self.pool.single_excitations + self.pool.double_excitations

        self.current_parameters = [0]

        self.max_steps = 200
        self.num_steps = 0

    def __repr__(self):

        output_message = '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&' + '\n'
        output_message += 'BACKEND PROPERTIES' + '\n'
        output_message += '- Number of qubits: ' + str(self.nqubits) + '\n'
        output_message += '- Backend: ' + str(self.backend) + '\n'
        output_message += '- Number of shots: ' + str(self.nshots) + '\n'
        output_message += '\n'
        output_message += 'SIMULATION COSTS' + '\n'
        output_message += '- Estimate for the cost of a single shot: ' + self.shot_cost + '\n'
        output_message += '- Estimate for the cost of launching a measurement: ' + str(self.launching_measurement_cost) + '\n'
        output_message += '- Estimate for the total simulation cost: ' + str(self.total_simulation_cost) + '\n'
        output_message += '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'        

        return output_message

    def beginning(self):

        output_message = '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&' + '\n'
        output_message += 'SIMULATION PROPERTIES' + '\n'
        output_message += 'Molecular system simulated: ' + self.name + '\n'
        output_message += 'Molecular orbital basis used: ' + 'sto-3g' + '\n'
        output_message += 'Number of orbital considered: ' + str(len(self.active_indices)) + '\n'
        output_message += 'Hartree-Fock energy of the system: ' + self.variational_energies[0] + '\n'
        output_message += '&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&'

    def verbose(self):

        print('--------------------------------------------------')
        print('Number of iterations accomplished so far: ', self.num_steps)
        print('Number of operators in the ansatz: ', len(self.generators))
        print('List of operators', self.list_of_operators)
        print('Variational energies: ', self.variational_energies)
        print('Current optimal parameters: ', self.current_parameters)
        print('--------------------------------------------------')
    
    def results(self):

        print('--------------------------------------------------')
        if self.num_steps != self.max_steps:
            print('CONVERGED!')
            print('Number of iterations until convergence: ', self.num_steps)
        else:
            print('DID NOT CONVERGED!')
            print(f'The maximal number {self.max_steps} has been reached! Need of algorithm review')
        print('Number of operators in the ansatz: ', len(self.generators))
        print('List of operators', self.list_of_operators)
        print('Final variational energy: ', self.variational_energies[-1])
        print('Final optimal parameters: ', self.optimal_parameters)
        print('--------------------------------------------------')

    def compose_ansatz(self):

        circuit = copy(self.current_state)
        
        for i in range(len(self.excitations)):
            circuit.append(self.excitations[i], self.n_qubits)
            self.list_of_operators.append(self.pool.names[i])
        
        return 

    def vqe(self):

        #Adding HF energy as the first variational energy
        self.variational_energies.append(uf.expectation_value(self.qiskit_ham, copy(self.current_state)))
        self.beginning()

        local = False

        self.compose_ansatz()
        while not(local) :
            print('Number of operators in the ansatz: ', len(self.generators))
            self.vqe_step()
            local = (np.abs(self.variational_energies[-1] - self.fci_energy) < self.energy_threshold)
            
            local = local and (self.num_steps >= self.max_steps)

            self.num_steps += 1
            self.verbose()
        
        self.optimal_parameters = self.current_parameters
        self.results()

        return
    
    def vqe_step(self):

        x0 = list(self.current_parameters)

        optimization = minimize(loss_function, x0=x0, method=self.optimizer)
        self.current_parameters = optimization.x
        self.variational_energies.append(optimization.fun)

        return

    def loss_function(self, params):
        
        return uf.expectation_value(self.qiskit_ham, self.current_state(params), self.nshots)
