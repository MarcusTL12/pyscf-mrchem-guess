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


def print_mrchem_bas_file(atoms, basis_sets):
    with open("mrchem.bas", "w") as f:
        f.write("Gaussian basis from ORCA gbw json file\n")
        f.write(f"{len(basis_sets):9}\n")

        for symbol, q, basis, atom_idxs in basis_sets:
            f.write(f"{q:9} {len(atom_idxs):4}")
            f.write(f" {len(basis):4}")
            for shells in basis.values():
                f.write(f" {len(shells):4}")
            f.write("\n")

            for i in atom_idxs:
                coord = atoms[i]["Coords"]
                f.write(f"{symbol} ")
                for x in coord:
                    f.write(f" {x:18.10f}")
                f.write("\n")

            for shells in basis.values():
                for exps, coeffs in shells:
                    n_prim = len(exps)
                    n_cont = len(coeffs)

                    f.write(f"{n_prim:9} {n_cont:4}\n")
                    for i in range(n_prim):
                        f.write(f"{exps[i]:15.7f}")
                        for j in range(n_cont):
                            f.write(f" {coeffs[j][i]:11.8f}")
                        f.write("\n")


def organize_basis(atoms):
    # Basis as key, atom index as value
    basis_sets = {}

    for atom in atoms:
        basis_hash = str((atom["ElementLabel"], atom["Basis"]))
        if basis_hash not in basis_sets:
            basis_sets[basis_hash] = (
                atom["ElementLabel"], atom["NuclearCharge"],
                consolidate_basis(atom["Basis"]), [atom["Idx"]]
            )
        else:
            basis_sets[basis_hash][3].append(atom["Idx"])

    basis_sets = list(basis_sets.values())

    # print(basis_sets)

    print_mrchem_bas_file(atoms, basis_sets)

    atoms_perm = []

    for _, _, _, idxs in basis_sets:
        for i in idxs:
            atoms_perm.append(i)

    print("Atom perm:", atoms_perm)


if __name__ == "__main__":
    with open(sys.argv[1]) as orca_file:
        orca_json = json.load(orca_file)

    organize_basis(orca_json["Molecule"]["Atoms"])

    # atoms = read_geometry(orca_file)
    # basis = read_basis(orca_file)

    # print_mrchem_bas_file(atoms, basis)

    # mo_perm = make_basis_permutation(atoms, basis)

    # mo = read_occ_mo(orca_file, mo_perm)

    # print_mrchem_mo_files(mo, len(mo_perm))
