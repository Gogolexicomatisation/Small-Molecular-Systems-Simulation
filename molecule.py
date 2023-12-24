class Molecule():
    """Molecular setup based on Openfermion"""

    def __init__(self, mol_name, bond_length) -> None:

        from qiskit.circuit import QuantumCircuit

        from openfermion import jordan_wigner
        from openfermion.chem import MolecularData
        from openfermionpyscf import run_pyscf

  #<<<<<<<<<<<<<<<<<<<<     Molecular Data     >>>>>>>>>>>>>>>>>>>>>>>>>>>

        self.name = mol_name

        if self.name == "H2" : atom_list = ['H', 'H']; coord_list = [(0.0, 0.0, 0.0), (0.0, 0.0, bond_length)]
        elif self.name == "H4" : atom_list = ['H', 'H', 'H', 'H']; coord_list = [(0.0, 0.0, 0.0), (0.0, 0.0, bond_length), (0.0, 0.0, 2*bond_length), (0.0, 0.0, 3*bond_length)]
        elif self.name == "LiH" : atom_list = ['Li', 'H']; coord_list = [(0.0, 0.0, 0.0), (0.0, 0.0, bond_length)]
        elif self.name == 'BeH2' : atom_list = ['Be', 'H', 'H']; coord_list = [(0.0, 0.0, 0.0), (0.0, 0.0, -bond_length/2), (0.0, 0.0, bond_length/2)]

        else: raise ValueError("Unkown molecule name. Please check the contents of the Molecule class")

        basis = 'sto-3g'
        multiplicity = 1
        charge = 0

        geometry = [[atom_list[i], coord_list[i]] for i in range(len(atom_list))]

        self.molecule = MolecularData(
            geometry=geometry,
            basis=basis,
            multiplicity=multiplicity,
            charge=charge,
            filename=f"{self.name}_openfermion.hdf5"
        )

        #run pyscf to get Hartree-Fock or Full CI wave function from it
        run_pyscf(self.molecule,
                  run_scf=True,
                  run_mp2=False,
                  run_cisd=False,
                  run_ccsd=True,
                  run_fci=True,
                  verbose=True
        )        

        print(self.molecule.n_orbitals)

        self.occupied_indices = range(0)
        self.active_indices = range(0,self.molecule.n_orbitals)

        #Creating the Qubit Hamiltonian
        hamiltonian = self.molecule.get_molecular_hamiltonian(occupied_indices=self.occupied_indices, active_indices=self.active_indices)
        self.hamiltonian = jordan_wigner(hamiltonian)

        self.hf_energy = self.molecule.hf_energy
        self.fci_energy = self.molecule.fci_energy

        self.n_electrons = self.molecule.n_electrons
        self.n_orbitals = 2*len(self.active_indices)
        self.n_qubits = 2*len(self.active_indices)

        #HF initial state 
        initial_state = QuantumCircuit(self.n_qubits)
        for k in range(self.n_qubits - self.n_electrons, self.n_qubits):
            initial_state.x(k)

        self.initial_state = initial_state

        return 
