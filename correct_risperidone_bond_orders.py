from rdkit import Chem

# The ligand's SMILES string, which defines its correct chemical bonds.
# (Replace this with the actual SMILES for your ligand).
smiles_template = "YOUR_LIGANDS_SMILES_STRING"

# Load the PDB file you extracted from your MD simulation.
md_structure = Chem.MolFromPDBFile("your_md_snapshot.pdb", removeHs=False)

# Create the template molecule from the SMILES string.
template = Chem.MolFromSmiles(smiles_template)

# Use the template to assign the correct bond orders to your MD structure.
# This keeps your coordinates but fixes the chemistry.
corrected_ligand = Chem.AllChem.AssignBondOrdersFromTemplate(template, md_structure)

# Save the fixed structure to a new file (SDF is a good format).
writer = Chem.SDWriter("corrected_ligand.sdf")
writer.write(corrected_ligand)
writer.close()

print("Corrected ligand structure saved to corrected_ligand.sdf")
