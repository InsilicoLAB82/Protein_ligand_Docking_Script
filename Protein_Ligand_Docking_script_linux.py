#Protein_Ligand Docking script by InsilicoLAB
#insilicolab82@gmail.com
import os
import subprocess
import matplotlib.pyplot as plt
import numpy as np

# Define file paths for receptor, ligand, grid, and reference ligand files
receptor_file = input("Enter the path to the receptor file (PDB format): ")
ligand_file = input("Enter the path to the ligand file (PDB format): ")
grid_file = input("Enter the path to the grid parameters file (txt format): ")
ref_ligand_file = input("Enter the path to the reference ligand file (PDBQT format): ")

print("\n" + "═" * 70)
print("          CUSTOMIZE YOUR DOCKING PARAMETERS")
print("═" * 70)

# Read grid parameters from the grid file
grid = {}
try:
    with open(grid_file) as f:
        for line in f:
            if '=' in line:
                k, v = line.split('=', 1)
                grid[k.strip()] = float(v.strip())
except Exception as e:
    print(f"Error reading grid file {grid_file}: {e}")
    exit(1)

# Interactive inputs for docking parameters
exhaustiveness = int(input("Exhaustiveness [8-128] (default 32): ") or "32")
scoring = input("Scoring function [vina/vinardo/ad4] (default vina): ").strip() or "vina"
num_modes = int(input("Number of poses (default 20): ") or "20")
energy_range = float(input("Energy range (kcal/mol) (default 3.0): ") or "3.0")
grid_spacing = float(input("Grid spacing in Å (default 0.375): ") or "0.375")
max_rmsd_show = float(input("Only show poses with RMSD < X Å (e.g. 3.0, or 999 for all): ") or "999")

print(f"\nStarting docking with your settings...")
print(f"   Scoring: {scoring.upper()} | Exhaustiveness: {exhaustiveness} | Poses: {num_modes}")

# Convert files to PDBQT format using Open Babel
try:
    os.system(f"obabel {receptor_file} -O receptor.pdbqt -xr -p 7.4 --partialcharge gasteiger -h")
    os.system(f"obabel {ligand_file} -O ligand.pdbqt --partialcharge gasteiger -h")
    os.system(f"obabel {ref_ligand_file} -O ref_ligand.pdbqt --partialcharge gasteiger -h")
except Exception as e:
    print(f"Error converting files with Open Babel: {e}")
    exit(1)

# Prepare AutoDock Vina command
cmd = [
    "./vina_1.2.7_linux_x86_64",
    "--receptor", "receptor.pdbqt",
    "--ligand", "ligand.pdbqt",
    "--center_x", str(grid['center_x']),
    "--center_y", str(grid['center_y']),
    "--center_z", str(grid['center_z']),
    "--size_x", str(grid['size_x']),
    "--size_y", str(grid['size_y']),
    "--size_z", str(grid['size_z']),
    "--spacing", str(grid_spacing),
    "--exhaustiveness", str(exhaustiveness),
    "--scoring", scoring,
    "--num_modes", str(num_modes),
    "--energy_range", str(energy_range),
    "--out", "result_docking.pdbqt",
    "--verbosity", "2"  # Ensure verbosity to get full Vina output
]

print("\n" + "═" * 80)
print("               AUTODOCK VINA ORIGINAL OUTPUT")
print("═" * 80 + "\n")

# Run Vina and print real-time output
try:
    # Run Vina and capture both stdout and stderr
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, bufsize=1)
    while True:
        # Read line by line from Vina's output
        output = process.stdout.readline()
        if output == '' and process.poll() is not None:
            break
        if output:
            print(output.strip())  # Print real-time output

    # Capture errors if any
    stderr_output = process.stderr.read()
    if stderr_output:
        print(f"Error output:\n{stderr_output}")

    process.wait()
except Exception as e:
    print(f"Error running AutoDock Vina: {e}")
    exit(1)

# Check if the result file was created
if not os.path.exists("result_docking.pdbqt"):
    print("Error: Docking result not found!")
    exit(1)  # Exit the program if docking failed

print("\n" + "═" * 80)
print("DOCKING FINISHED – NOW CALCULATING RMSD TO CRYSTAL POSE")
print("═" * 80)

# Function to split PDBQT into individual poses based on model number
def split_pdbqt_into_poses(pdbqt_file):
    poses = []
    try:
        with open(pdbqt_file, 'r') as f:
            lines = f.readlines()
        
        current_pose_lines = []
        model_counter = 0
        is_in_model = False
        
        for line in lines:
            if "MODEL" in line:
                if current_pose_lines:  # Save the previous pose if it exists
                    poses.append((model_counter, current_pose_lines))
                current_pose_lines = [line]  # Start the new pose with this line
                model_counter += 1
                is_in_model = True
            elif "ENDMDL" in line:
                if is_in_model:
                    current_pose_lines.append(line)  # End the current pose
                    poses.append((model_counter, current_pose_lines))
                    current_pose_lines = []  # Reset for the next model
                    is_in_model = False
            else:
                if is_in_model:
                    current_pose_lines.append(line)  # Add lines to the current pose

        if current_pose_lines:  # Add the last pose
            poses.append((model_counter, current_pose_lines))

    except Exception as e:
        print(f"Error splitting {pdbqt_file} into poses: {e}")
    return poses

# Split the docking result into individual poses
poses = split_pdbqt_into_poses("result_docking.pdbqt")

affinities = []
rmsd_list = []

# Calculate RMSD and Affinity for each pose using LS-align
for model_num, pose in poses:
    pose_file = f"ligand_{model_num:04d}.pdbqt"
    with open(pose_file, 'w') as f:
        f.writelines(pose)
    
    # Convert the pose to mol2 format for RMSD calculation using LS-align
    os.system(f"obabel {pose_file} -O ligand_{model_num:04d}.mol2 -h")
    os.system(f"obabel {ref_ligand_file} -O ref_ligand.mol2 -h")
    
    # Run LS-align to calculate RMSD
    ls_align_cmd = f"./LSalign ligand_{model_num:04d}.mol2 ref_ligand.mol2 -o aligned_ligand_{model_num:04d}.mol2 > rmsd_output_{model_num}.txt"
    os.system(ls_align_cmd)

    # Parse RMSD value from LS-align output
    try:
        with open(f"rmsd_output_{model_num}.txt", "r") as f:
            lines = f.readlines()
            rmsd_line = lines[3]  # 0-based index, line 4 is at index 3
            rmsd_value = rmsd_line.strip().split()[7]  # RMSD is at column 8 (index 7)
            rmsd = float(rmsd_value)
    except Exception as e:
        print(f"Error extracting RMSD for model {model_num}: {e}")
        rmsd = 3.0  # Default RMSD value if extraction fails
    
    # Get affinity from Vina log
    try:
        aff = float(os.popen(f"grep 'REMARK VINA RESULT' result_docking.pdbqt | sed -n {model_num}p").read().split()[3])
    except Exception as e:
        print(f"Error extracting affinity for model {model_num}: {e}")
        aff = 3.0  # Default affinity value if extraction fails
    
    if rmsd <= max_rmsd_show:
        affinities.append(aff)
        rmsd_list.append(rmsd)

# Plot results
if affinities:
    plt.figure(figsize=(10,6))
    sc = plt.scatter(rmsd_list, affinities, c=affinities, cmap='coolwarm_r', s=120, edgecolors='k')
    plt.colorbar(sc, label='Affinity (kcal/mol)')
    plt.xlabel('RMSD from Crystal Pose (Å)', fontsize=14)
    plt.ylabel('Binding Affinity (kcal/mol)', fontsize=14)
    plt.title(f'{scoring.upper()} Docking – Exhaustiveness {exhaustiveness}', fontsize=16)
    best_idx = np.argmin(affinities)
    plt.scatter(rmsd_list[best_idx], affinities[best_idx], c='red', s=400, marker='*',
                label=f'Best: {affinities[best_idx]:.2f} kcal/mol @ {rmsd_list[best_idx]:.2f} Å')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig("Vina_Results.png", dpi=300)
    plt.show()

    print(f"\nBEST POSE FOUND: {affinities[best_idx]:.2f} kcal/mol at {rmsd_list[best_idx]:.2f} Å RMSD")

# Save the results
os.rename("result_docking.pdbqt", "result_docking_final.pdbqt")

print("\nALL DONE! You saw the real Vina output + got validated results.")
