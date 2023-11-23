from Bio import PDB
from Bio.PDB.DSSP import DSSP
from sys import argv
import numpy as np
import matplotlib.pyplot as plt


ramachandran_type = {
    "H": 1,
    "B": 2,
    "E": 3,
    "G": 4,
    "I": 5,
    "T": 6,
    "S": 7,
    "-": 8,
    "P": 8,
    "!": 8,
}

# check if the user has provided a PDB file
if len(argv) != 2:
    print("Usage: python main.py <pdb file>")
    exit(0)

# parse the PDB file and get the structure
parser = PDB.PDBParser(QUIET=True)
pdb_code = argv[1][-8:-4]
structure = parser.get_structure(pdb_code, argv[1])
ramachandran = np.zeros((360, 360))

for model in structure:
    dssp = DSSP(model, argv[1])
    for res in dssp:
        if int(res[4]) == 360 or int(res[5]) == 360:
            continue

        ramachandran[int(res[5]) + 180][int(res[4]) + 180] = ramachandran_type[res[2]]

fig, ax = plt.subplots()
types = [
    "Alpha-helix (4-12)",
    "Isolated beta-bridge residue",
    "Strand",
    "3-10 helix",
    "Pi helix",
    "Turn",
    "Bend",
    "None",
]
plt.title("Ramachandran plot for " + pdb_code)
for i in range(1, len(types) + 1):
    indices = np.where(ramachandran == i)
    ax.scatter(indices[1] - 180, indices[0] - 180, label=types[i - 1], marker=".")
ax.grid(True, linestyle="--", alpha=0.7)
plt.legend(bbox_to_anchor=(1.04, 1), loc="upper left")
plt.xlabel(r"$\phi$")
plt.ylabel(r"$\psi$")
ax.set_xlim(-180, 180)
ax.set_ylim(-180, 180)
plt.savefig(str(argv[1][-8:-4]) + "_ramachandran.png", bbox_inches="tight", dpi=300)
# plt.show()
