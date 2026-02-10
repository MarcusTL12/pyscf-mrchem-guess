import pyscf
from pprint import pprint
import sys

mol = pyscf.M(atom="""
              O 0 0 0
              H 1 0 0
              H 0 1 0
              """, basis="cc-pVDZ")

ovlp = mol.intor("int1e_ovlp")

pyscf.tools.dump_mat.dump_tri(sys.stdout, ovlp, label=mol.ao_labels(), ncol=6)
