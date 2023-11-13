from Bio import PDB
from sys import argv
import matplotlib.pyplot as plt
import numpy as np

# check if the user has provided a PDB file
if len(argv) != 2:
    print("Usage: python main.py <pdb file>")
    exit(0)

# parse the PDB file and get the structure
parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure(argv[1][-8:-4:1], argv[1])
# setup variables and lists

threshold = 8
contact_map = np.zeros((10000, 10000))
residueList = []
x, y = 0, 0

# iterate over all residues in the structure
for residue1 in structure.get_residues():
    # skip if residue has no CA atom
    if not residue1.has_id("CA"):
        continue

    # restart the y index
    y = 0

    # add residue to list (for counting)
    residueList.append([str(x), residue1.get_resname(), residue1["CA"].get_coord()])

    # iterate over all residues again to compare with the first residue
    for residue2 in structure.get_residues():
        # skip if residue has no CA atom
        if not residue2.has_id("CA"):
            continue

        # get the distance between the two residues and update the contact map
        if residue1["CA"] - residue2["CA"] <= threshold:
            contact_map[x, y] = 1
        y += 1
    x += 1

# trim the array to the correct size
contact_map = contact_map[0 : len(residueList), 0 : len(residueList)]

# show the contact map
plt.imshow(contact_map, cmap="binary", interpolation="nearest")
# plt.gca().invert_yaxis()
plt.savefig(str(argv[1][-8:-4:1]) + "_contact_map.png", bbox_inches="tight")
plt.show()
