
import MDAnalysis as mda
from MDAnalysis.analysis import distances
import numpy as np
import os
from MDAnalysis.lib.distances import distance_array

def analyze_cc_bonds_in_trajectory(pdb_file, xtc_file, process_percentage=100.0):
    """
    Analyzes a GROMACS trajectory to calculate the mean and standard deviation
    of C-C bond lengths over time.

    Args:
        pdb_file (str): Path to the PDB file (for topology).
        xtc_file (str): Path to the XTC trajectory file.
        process_percentage (float): The percentage of frames to process.

    Returns:
        A dictionary containing mean and std dev for each C-C bond.
    """
    if not os.path.exists(pdb_file):
        print(f"Error: PDB file not found at {pdb_file}")
        return None
    if not os.path.exists(xtc_file):
        print(f"Error: XTC file not found at {xtc_file}")
        return None

    try:
        # Load the universe
        u = mda.Universe(pdb_file, xtc_file)
        
        # Select all carbon atoms
        # Select only carbon atoms belonging to the '8NU' ligand
        carbons = u.select_atoms("resname 8NU and type C")
        
        # Find all C-C bonds more efficiently using MDAnalysis.lib.distances.distance_array.
        # This calculates all pairwise distances between carbon atoms in the first frame
        # and identifies bonds based on a distance cutoff.
        
        # Get positions of all carbon atoms in the first frame
        carbon_positions = carbons.positions
        
        # Calculate distance array between all carbon atoms
        # This will be a symmetric matrix where dist_matrix[i, j] is distance between carbon i and carbon j
        dist_matrix = distance_array(carbon_positions, carbon_positions)
        
       	cc_bonds = []
        # Iterate through the upper triangle of the distance matrix to find unique pairs
        for i in range(len(carbons)):
            for j in range(i + 1, len(carbons)): # Only check upper triangle to avoid duplicates and self-distances
                if dist_matrix[i, j] < 1.7: # A reasonable cutoff for a C-C bond is ~1.7 Angstroms
                    cc_bonds.append((carbons[i], carbons[j]))

        if not cc_bonds:
            print("No C-C bonds found within the cutoff distance.")
            return None

        # Prepare for analysis
        num_frames = len(u.trajectory)
        bond_distances = [[] for _ in cc_bonds]

        # Determine the frame step based on the percentage
        if process_percentage >= 100.0 or process_percentage <= 0:
            frame_step = 1
        else:
            frame_step = int(100.0 / process_percentage)
        
       	print(f"Processing approximately {process_percentage:.2f}% of frames (1 in every {frame_step} frames).")

        # Iterate through the trajectory with the specified step
        for ts in u.trajectory[::frame_step]:
            for i, (atom1, atom2) in enumerate(cc_bonds):
                dist = np.linalg.norm(atom1.position - atom2.position)
                bond_distances[i].append(dist)

        # Calculate mean and std dev
        results = {}
        for i, (atom1, atom2) in enumerate(cc_bonds):
            bond_name = f"C{atom1.id}-C{atom2.id}"
            mean_dist = np.mean(bond_distances[i])
            std_dist = np.std(bond_distances[i])
            results[bond_name] = (mean_dist, std_dist)
            
	return results

    except Exception as e:
        print(f"An error occurred: {e}")
        return None

def print_results_table(results):
    """Prints the analysis results in a formatted table."""
    if not results:
        print("No results to display.")
        return
        
    print(f"{'Bond':<12} | {'Mean Distance (Å)':<20} | {'Std Dev (Å)':<15}")
    print("-" * 55)
    for bond, values in results.items():
        print(f"{bond:<12} | {values[0]:<20.4f} | {values[1]:<15.4f}")

if __name__ == "__main__":
    # --- Configuration ---
    # Edit these paths to point to your files
    PDB_PATH = "step5_input.pdb"
    XTC_PATH = "step7_prot_SOLU_5percent.xtc"  # Make sure this is the correct trajectory file
    FRAMES_TO_PROCESS_PERCENT = 5.0 # Process only 5% of the frames
    # ---------------------

    print(f"Analyzing trajectory for {PDB_PATH} and {XTC_PATH}...")
    analysis_results = analyze_cc_bonds_in_trajectory(PDB_PATH, XTC_PATH, FRAMES_TO_PROCESS_PERCENT)
