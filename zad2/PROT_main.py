from Bio import PDB
from sys import argv
import numpy as np
import matplotlib.pyplot as plt
import math


def degrees(rad_angle):
    if rad_angle is None:
        return None
    angle = rad_angle * 180 / math.pi
    while angle > 180:
        angle = angle - 360
    while angle < -180:
        angle = angle + 360
    return angle


def ramachandran_type(psi, phi):
    zakres = 45

    # alpha helix
    if -45 - zakres <= psi <= -45 + zakres and -60 - zakres <= phi <= -60 + zakres:
        return 1

    # beta sheet
    elif 105 - zakres <= psi <= 155 + zakres and -105 - zakres <= phi <= -100 + zakres:
        return 2

    # left handed alpha helix
    elif -50 <= psi <= 100 and 30 <= phi <= 130:
        return 3

    # inne
    else:
        return 4


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
    for chain in model:
        polypeptides = PDB.CaPPBuilder().build_peptides(chain)
        for poly_index, poly in enumerate(polypeptides):
            phi_psi = poly.get_phi_psi_list()
            for res_index, residue in enumerate(poly):
                phi, psi = phi_psi[res_index]
                if phi and psi:
                    ramachandran[int(degrees(psi)) + 180][
                        int(degrees(phi)) + 180
                    ] = ramachandran_type(degrees(psi), degrees(phi))


fig, ax = plt.subplots()
# Define colors for each value (0, 1, 2, 3, 4)
colors = ["red", "blue", "green", "grey"]
types = ["Alpha-helix", "Beta-sheet", "Left-handed Alpha-Helix", "Others"]
plt.title("Ramachandran plot for " + pdb_code)
for i in range(1, 5):
    indices = np.where(ramachandran == i)
    ax.scatter(indices[1] - 180, indices[0] - 180, c=colors[i - 1], label=types[i - 1])
ax.grid(True, linestyle="--", alpha=0.7)
ax.legend()
plt.xlabel(r"$\phi$")
plt.ylabel(r"$\psi$")
ax.set_xlim(-180, 180)
ax.set_ylim(-180, 180)
plt.savefig(str(argv[1][-8:-4]) + "_ramachandran.png", bbox_inches="tight")
plt.show()
