
def run(vqe_config: dict) -> None:

    import os

    from VQE import MoleculeSimulator
    import molecule
    import useful_functions as uf

    from qiskit_aer import Aer

    #Setup molecule
    molecule_sys = molecule.Molecule(mol_name=vqe_config['molecule_name'],
                                    bond_length=vqe_config['bond_length'])
    assert molecule_sys.initial_state is not None
    assert molecule_sys.n_qubits > 1

    #Setup VQE
    backend = Aer.get_backend('aer_simulator')

    #Initialize simulator
    molSim = MoleculeSimulator(mol=molecule_sys, nshots=vqe_config['nshots'], backend=backend)

    #Run VQE

    molSim.vqe()

if __name__ == "__main__":

    import argparse

    #define the parser
    parser = argparse.ArgumentParser(
        epilog="VQE runner",
        usage="python main.py --help",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "-m",
        "--molecule_name",
        help="Name of the simulated molecule, e.g. H2",
        type=str,
        required=True,
    )
    parser.add_argument(
        "-b",
        "--bond_length",
        help="Bond-length in Angstrom (A)",
        type=float,
        required=True,
    )
    parser.add_argument(
        "-s",
        "--shots",
        help="Number of shots",
        type=int,
        required=False,
        default=1000,
    )
    parser.add_argument(
        "-e",
        "--excitations",
        help="Type of excitations used (UCCSD or QEB)",
        type=str,
        required=False,
        choices=['uccsd'],
        default='uccsd'
    )

    #parse the input values
    args = parser.parse_args()

    vqe_config = {
        'molecule_name'     : args.molecule_name,
        'bond_length'       : args.bond_length,
        'nshots'            : args.shots,
        'excitation_method' : args.excitations
    }

    #Run the VQE driver
    run(vqe_config)
    print("Exit main")
