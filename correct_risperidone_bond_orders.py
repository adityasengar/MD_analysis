from rdkit import Chem
from rdkit.Chem import AllChem

# The ligand's SMILES string, which defines its correct chemical bonds.
smiles_template = "Fc1ccc2c(c1)onc2C1CCN(CC1)CCc1c(C)nc2n(c1=O)CCCC2"

# Load the PDB file you extracted from your MD simulation.
pdb_file_path = "step5_input.pdb"
ligand_residue_name = "8NU"
temp_ligand_pdb_path = "risperidone_temp.pdb"

with open(pdb_file_path, 'r') as infile, open(temp_ligand_pdb_path, 'w') as outfile:
    for line in infile:
        if line.startswith("ATOM") or line.startswith("HETATM"):
            res_name = line[17:20].strip()
            if res_name == ligand_residue_name:
                outfile.write(line)

md_structure = Chem.MolFromPDBFile(temp_ligand_pdb_path, removeHs=True)

# Create the template molecule from the SMILES string.
template = Chem.MolFromSmiles(smiles_template)

# Use the template to assign the correct bond orders to your MD structure.
# This keeps your coordinates but fixes the chemistry.
corrected_ligand = Chem.AllChem.AssignBondOrdersFromTemplate(template, md_structure)

# Save the fixed structure to a new file (SDF is a good format).
writer = Chem.SDWriter("risperidone_corrected.sdf")
writer.write(corrected_ligand)
writer.close()

print("Corrected ligand structure saved to corrected_ligand.sdf")

