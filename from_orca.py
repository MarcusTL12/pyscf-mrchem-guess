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


if __name__ == "__main__":
    with open(sys.argv[1]) as orca_file:
        atoms = read_geometry(orca_file)

        pprint(atoms)

        basis = read_basis(orca_file)

        print(basis)

        # for l in orca_file:
        #     if l.strip() == "BASIS SET IN INPUT FORMAT":
        #         print("found 1")
        #         break

        # for l in orca_file:
        #     if l.strip() == "AUXILIARY/J BASIS SET INFORMATION":
        #         print("found 2")
        #         break
