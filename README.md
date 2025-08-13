# Ligand Bond Order Correction Scripts

This directory contains scripts for correcting the bond orders of ligand structures extracted from molecular dynamics (MD) simulations. MD simulations can sometimes produce PDB files where the chemical bond orders are not correctly assigned, especially for complex ligands. These scripts use RDKit to reassign correct bond orders based on a SMILES string template, while preserving the atomic coordinates from the MD snapshot.

## `correct_ligand_bond_orders_template.py`

This is a generic script that can be adapted to correct the bond orders for any ligand.

**Usage:**

1.  **Update `smiles_template`**: Replace `"YOUR_LIGANDS_SMILES_STRING"` with the actual SMILES string of your ligand. This SMILES string defines the correct chemical bonding.
2.  **Update `your_md_snapshot.pdb`**: Change `"your_md_snapshot.pdb"` to the path of the PDB file containing your ligand's structure from the MD simulation.
3.  **Run the script**:
    ```bash
    python correct_ligand_bond_orders_template.py
    ```

The script will output a new SDF file named `corrected_ligand.sdf` with the corrected bond orders.

## `correct_risperidone_bond_orders.py`

This script is specifically configured to correct the bond orders for the Risperidone ligand.

**Details:**

*   **SMILES Template**: The SMILES string for Risperidone is hardcoded within the script.
*   **Input PDB**: It reads the ligand structure from `step5_input.pdb`.
*   **Ligand Residue Name**: It specifically looks for a ligand with the residue name "8NU" within the PDB file.
*   **Output SDF**: The corrected Risperidone structure is saved to `risperidone_corrected.sdf`.

**Usage:**

Ensure `step5_input.pdb` is in the same directory and contains the Risperidone ligand with residue name "8NU".

```bash
python correct_risperidone_bond_orders.py
```

The script will output `risperidone_corrected.sdf` with the corrected bond orders for Risperidone.

