import pyscf
from pprint import pprint

periodic_table_str = """
H                                                  He
Li Be                               B  C  N  O  F  Ne
Na Mg                               Al Si P  S  Cl Ar
K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
"""

periodic_table_dict = {}

for i, s in enumerate(periodic_table_str.split()):
    periodic_table_dict[s] = i + 1


def sort_atom_str(atom_str: str):
    lines = atom_str.splitlines()

    sortable = [(periodic_table_dict[l.split()[0]], l)
                for l in lines if len(l.strip()) > 0]
    sortable.sort()

    lines = [l for _, l in sortable]

    return '\n'.join(lines)


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


def print_mrchem_mo_file(mol, C, filename):
    with open(filename, "w") as f:
        nao = mol.nao
        f.write(f"{nao:12}\n")
        for x in C.T.flatten():
            f.write(f"{x:20.15f}\n")


def localize_occ_mo(mol, hf):
    no = mol.nelectron // 2
    C = hf.mo_coeff[:, 0:no]

    return pyscf.lo.boys.Boys(mol).kernel(C)


xyz = """
O         4.40988174540998    7.45211310544178    6.80383413247153
H         4.02905772874383    7.65799771302110    5.92947407439696
C         3.47744106818372    7.98095365314221    7.77472943007089
H         3.37682956954364    9.06072722191245    7.65247353819164
H         3.87283422362894    7.76389288840846    8.76739659140869
H         2.50240699396907    7.50433645096497    7.66056354455503
"""

xyz = sort_atom_str(xyz)
print(xyz)

mol = pyscf.M(atom=xyz, basis="sto-6G", charge=0, verbose=4)

print_mrchem_bas_file(mol, "mrchem.bas")

hf = mol.RKS()
hf.xc = "pbe"
# hf = mol.HF()
hf.run()

hf.analyze()

C = hf.mo_coeff # canonical
# C = localize_occ_mo(mol, hf) # localized
# C = pyscf.numpy.identity(mol.nao)

print_mrchem_mo_file(mol, C, "mrchem.mop")
