from Bio import PDB
from sys import argv
from Bio.PDB import calc_dihedral

# check if the user has provided a PDB file
if len(argv) != 2:
    print("Usage: python main.py <pdb file>")
    exit(0)

# parse the PDB file and get the structure
parser = PDB.PDBParser(QUIET=True)
structure = parser.get_structure(argv[1][-8:-4], argv[1])


def get_atom_vector(residue, atom_name):
    if atom_name in residue:
        return residue[atom_name].get_vector()
    else:
        return None


def calculate_dihedral(v1, v2, v3, v4):
    if v1 is not None and v2 is not None and v3 is not None and v4 is not None:
        return calc_dihedral(v1, v2, v3, v4)
    else:
        return "N/A"


headers = ["resname", "alpha", "beta", "gamma", "delta", "epsilon", "zeta", "chi"]
with open("dihedrals.csv", "w") as f:
    f.write(",".join(headers) + "\n")
    for chain in structure.get_chains():
        residues = list(chain)
        i = 0
        while i < len(residues):
            residue = residues[i]
            residue_name = residue.get_resname()
            # pomin wode
            if residue_name == "HOH":
                i += 1
                continue
            if i > 0:
                prev_residue = residues[i - 1]
                if "O3'" in prev_residue:
                    o3p = prev_residue["O3'"].get_vector()
                else:
                    o3p = None
            else:
                o3p = None

            if i < len(residues) - 1:
                next_residue = residues[i + 1]
                if "O5'" in next_residue:
                    o5n = next_residue["O5'"].get_vector()
                else:
                    o5n = None
                if "P" in next_residue:
                    pn = next_residue["P"].get_vector()
                else:
                    pn = None
            else:
                o5n = None
                pn = None

            p = get_atom_vector(residue, "P")
            o5 = get_atom_vector(residue, "O5'")
            c5 = get_atom_vector(residue, "C5'")
            c4 = get_atom_vector(residue, "C4'")
            c3 = get_atom_vector(residue, "C3'")
            o3 = get_atom_vector(residue, "O3'")
            o4 = get_atom_vector(residue, "O4'")
            c1 = get_atom_vector(residue, "C1'")
            c4w = get_atom_vector(residue, "C4")
            c2w = get_atom_vector(residue, "C2")
            n = get_atom_vector(
                residue, "N9" if residue_name == "A" or residue_name == "G" else "N1"
            )

            alpha = calculate_dihedral(o3p, p, o5, c5)
            beta = calculate_dihedral(p, o5, c5, c4)
            gamma = calculate_dihedral(o5, c5, c4, c3)
            delta = calculate_dihedral(c5, c4, c3, o3)
            epsilon = calculate_dihedral(c4, c3, o3, p)
            zeta = calculate_dihedral(c3, o3, p, o5)
            chi = calculate_dihedral(
                o4, c1, n, c4w if residue_name == "A" or residue_name == "G" else c2w
            )

            angles = [residue_name, alpha, beta, gamma, delta, epsilon, zeta, chi]
            line = ",".join(map(str, angles))
            f.write(line + "\n")
            i += 1
