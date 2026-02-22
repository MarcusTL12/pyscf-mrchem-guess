import sys
from pprint import pprint


def read_geometry(orca_file):
    for l in orca_file:
        if l.strip() == "CARTESIAN COORDINATES (A.U.)":
            break

    orca_file.readline()
    orca_file.readline()

    atoms = []

    for l in orca_file:
        if len(l.strip()) == 0:
            break

        _, sym, q, _, _, x, y, z = l.split()

        atoms.append((sym, float(q), [float(x), float(y), float(z)]))

    return atoms


def read_basis(orca_file):
    for l in orca_file:
        if l.strip() == "BASIS SET IN INPUT FORMAT":
            break

    orca_file.readline()

    basis = {}

    for l in orca_file:
        if l.startswith('-'):
            break
        elif l.strip().startswith("NewGTO"):
            _, element = l.split()

            basis[element] = {}

            for l in orca_file:
                if l.strip() == "end;":
                    break
                angmom, nprim = l.split()
                nprim = int(nprim)

                if angmom not in basis[element]:
                    basis[element][angmom] = []

                exps = []
                coeffs = []

                for _ in range(nprim):
                    _, e, c = orca_file.readline().split()
                    exps.append(float(e))
                    coeffs.append(float(c))

                if len(basis[element][angmom]) != 0 and \
                        basis[element][angmom][-1][0] == exps:
                    basis[element][angmom][-1][1].append(coeffs)
                else:
                    basis[element][angmom].append((exps, [coeffs]))

    return basis


def print_mrchem_bas_file(atoms, basis):
    with open("mrchem.bas", "w") as f:
        f.write("Gaussian basis from ORCA output file\n")
        f.write(f"{len(basis):9}\n")

        atom_coord_dict = {}

        for symbol, q, coord in atoms:
            if symbol in atom_coord_dict:
                atom_coord_dict[symbol][1].append(coord)
            else:
                atom_coord_dict[symbol] = (q, [coord])

        for symbol in atom_coord_dict:
            q, coords = atom_coord_dict[symbol]
            f.write(f"{q:9}. {len(coords):4}")
            funcs_per_shell = {}
            for angmom in basis[symbol]:
                if angmom in funcs_per_shell:
                    funcs_per_shell[angmom] += 1
                else:
                    funcs_per_shell[angmom] = 1
            f.write(f" {len(basis[symbol]):4}")
            for angmom in basis[symbol]:
                f.write(f" {len(basis[symbol][angmom]):4}")
            f.write("\n")

            for coord in coords:
                f.write(f"{symbol} ")
                for x in coord:
                    f.write(f" {x:18.10f}")
                f.write("\n")

            for angmom in basis[symbol]:
                for exps, coeffs in basis[symbol][angmom]:
                    n_prim = len(exps)
                    n_cont = len(coeffs)

                    f.write(f"{n_prim:9} {n_cont:4}\n")
                    for i in range(n_prim):
                        f.write(f"{exps[i]:15.7f}")
                        for j in range(n_cont):
                            f.write(f" {coeffs[j][i]:11.8f}")
                        f.write("\n")


def shell_permutation(angmom):
    if angmom == 'S':
        return [0]
    elif angmom == 'P':
        return [1, 2, 0]  # zxy -> xyz
    elif angmom == 'D':
        # (z2, xz, yz, x2y2, xy) -> (xy, yz, z2, xz, x2y2)
        return [4, 2, 0, 1, 3]
    else:
        raise NotImplementedError(f"{angmom} orbitals not implemented")


def atom_basis_permutation(atom_basis):
    perm = []

    for angmom in atom_basis:
        shell_perm = shell_permutation(angmom)

        num_funcs = 0
        for _, coeffs in atom_basis[angmom]:
            num_funcs += len(coeffs)

        for _ in range(num_funcs):
            offset = len(perm)
            for i in shell_perm:
                perm.append(i + offset)

    return perm


def make_basis_permutation(atoms, basis):
    atom_type_order = {}

    for s, _, _ in atoms:
        if s not in atom_type_order:
            atom_type_order[s] = len(atom_type_order)

    atoms_ordered = [(atom_type_order[s], i)
                     for i, (s, _, _) in enumerate(atoms)]
    atoms_ordered.sort()

    atoms_perm = [i for _, i in atoms_ordered]

    atom_basis_perms = {s: atom_basis_permutation(basis[s]) for s in basis}

    atoms_offsets = []
    offset = 0
    for s, _, _ in atoms:
        atoms_offsets.append(offset)
        offset += len(atom_basis_perms[s])

    perm = []

    for i in atoms_perm:
        s = atoms[i][0]
        for j in atom_basis_perms[s]:
            perm.append(j + atoms_offsets[i])

    return perm


if __name__ == "__main__":
    with open(sys.argv[1]) as orca_file:
        atoms = read_geometry(orca_file)
        basis = read_basis(orca_file)

        print_mrchem_bas_file(atoms, basis)

        mo_perm = make_basis_permutation(atoms, basis)

        print(mo_perm)

        # for l in orca_file:
        #     if l.strip() == "BASIS SET IN INPUT FORMAT":
        #         print("found 1")
        #         break

        # for l in orca_file:
        #     if l.strip() == "AUXILIARY/J BASIS SET INFORMATION":
        #         print("found 2")
        #         break
