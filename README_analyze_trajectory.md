# `analyze_trajectory.py`

This script analyzes GROMACS molecular dynamics trajectories to calculate the mean and standard deviation of carbon-carbon (C-C) bond lengths over time, specifically for the ligand in the system.

## Features

-   Loads GROMACS topology (PDB) and trajectory (XTC) files.
-   Identifies C-C bonds within the specified ligand (currently `UNK`).
-   Calculates the mean and standard deviation of each identified C-C bond length across the trajectory.
-   Allows processing a percentage of frames to speed up analysis.
-   Prints results in a formatted table.

## How to Use

1.  **Place the script:** Ensure `analyze_trajectory.py` is in your working directory.
2.  **Update Configuration:** Open `analyze_trajectory.py` and modify the `PDB_PATH`, `XTC_PATH`, and `FRAMES_TO_PROCESS_PERCENT` variables in the `if __name__ == "__main__":` block to point to your specific files and desired processing percentage.

    ```python
    # --- Configuration ---
    PDB_PATH = "path/to/your/topology.pdb"
    XTC_PATH = "path/to/your/trajectory.xtc"
    FRAMES_TO_PROCESS_PERCENT = 10.0 # e.g., process 10% of frames
    # ---------------------
    ```

3.  **Run the script:** Execute the script from your terminal:

    ```bash
    python analyze_trajectory.py
    ```

## Output

The script will print a table to the console showing each identified C-C bond, its mean length, and the standard deviation of its length over the analyzed trajectory frames.

```
Bond         | Mean Distance (Å)    | Std Dev (Å)    
-------------------------------------------------------
C1-C2        | 1.5300               | 0.0200         
C2-C3        | 1.5450               | 0.0150         
...
```

## Future Enhancements

-   **Automatic Bond Type Assignment:** Implement functionality to automatically assign bond types (e.g., single, double, triple) based on the calculated mean bond lengths.
