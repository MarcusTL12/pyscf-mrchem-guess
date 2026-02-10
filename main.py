import pyscf
from pprint import pprint


def print_mrchem_bas_file(mol, filename):
    with open(filename, "w") as f:
        # Write comment line
        f.write(f"Gaussian basis {mol.basis} from PySCF\n")

        bas_dict = mol._basis
        f.write(f"{len(bas_dict):9}\n")

        atom_coord_dict = {}

        for i in range(mol.natm):
            symbol = mol.atom_symbol(i)
            coord = mol.atom_coord(i)

            if symbol in atom_coord_dict:
                atom_coord_dict[symbol][1].append(coord)
            else:
                atom_coord_dict[symbol] = (mol.atom_charge(i), [coord])

        for symbol in atom_coord_dict:
            q, coords = atom_coord_dict[symbol]
            f.write(f"{q:9}. {len(coords):4}")
            funcs_per_shell = []
            for func in bas_dict[symbol]:
                l = func[0]
                if len(funcs_per_shell) > l:
                    funcs_per_shell[l] += 1
                else:
                    funcs_per_shell.append(1)
            f.write(f" {len(funcs_per_shell):4}")
            for nfuncs in funcs_per_shell:
                f.write(f" {nfuncs:4}")
            f.write("\n")

            for coord in coords:
                f.write(f"{symbol} ")
                for x in coord:
                    f.write(f" {x:18.10f}")
                f.write("\n")

            for func in bas_dict[symbol]:
                l = func[0]

                n_prim = len(func) - 1
                n_cont = len(func[1]) - 1

                f.write(f"{n_prim:9} {n_cont:4}\n")
                for i in range(1, n_prim+1):
                    f.write(f"{func[i][0]:15.7f}")
                    for j in range(n_cont):
                        f.write(f" {func[i][j + 1]:11.8f}")
                    f.write("\n")


mol = pyscf.M(atom="""
O       -2.61568        1.24534        0.01028
H       -1.62635        1.28445       -0.00726
H       -2.90504        2.10946       -0.37724
""", basis="cc-pVDZ")

print("\n" + "#" * 80)
pprint(mol._basis)
print("#" * 80 + "\n")

print_mrchem_bas_file(mol, "tmp.bas")
