import sys
import json
from pprint import pprint


def consolidate_basis(basis):
    new_basis = {}

    for shell in basis:
        angmom = shell["Shell"]
        exp_hash = str(shell["Exponents"])

        if angmom not in new_basis:
            new_basis[angmom] = {}

        if exp_hash not in new_basis[angmom]:
            new_basis[angmom][exp_hash] = (shell["Exponents"],
                                           [shell["Coefficients"]])
        else:
            new_basis[angmom][exp_hash][1].append(shell["Coefficients"])

    return {angmom: list(shell_basis.values())
            for angmom, shell_basis in new_basis.items()}


def organize_basis(orca_json):
    # Basis as key, atom index as value
    basis_sets = {}

    for atom in orca_json["Molecule"]["Atoms"]:
        basis_hash = str(atom["Basis"])
        if basis_hash not in basis_sets:
            basis_sets[basis_hash] = (atom["Basis"], [atom["Idx"]])
        else:
            basis_sets[basis_hash][1].append(atom["Idx"])

    atoms_perm = []

    for _, idxs in basis_sets.values():
        for i in idxs:
            atoms_perm.append(i)


if __name__ == "__main__":
    with open(sys.argv[1]) as orca_file:
        orca_json = json.load(orca_file)

    organize_basis(orca_json)

    # atoms = read_geometry(orca_file)
    # basis = read_basis(orca_file)

    # print_mrchem_bas_file(atoms, basis)

    # mo_perm = make_basis_permutation(atoms, basis)

    # mo = read_occ_mo(orca_file, mo_perm)

    # print_mrchem_mo_files(mo, len(mo_perm))
