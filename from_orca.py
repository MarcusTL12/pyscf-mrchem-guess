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

    basis = {}

    


if __name__ == "__main__":
    with open(sys.argv[1]) as orca_file:
        atoms = read_geometry(orca_file)

        pprint(atoms)

        # for l in orca_file:
        #     if l.strip() == "BASIS SET IN INPUT FORMAT":
        #         print("found 1")
        #         break

        # for l in orca_file:
        #     if l.strip() == "AUXILIARY/J BASIS SET INFORMATION":
        #         print("found 2")
        #         break
