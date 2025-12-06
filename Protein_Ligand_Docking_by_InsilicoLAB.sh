#!/bin/bash
# Protein–Ligand Docking Pipeline – InsilicoLAB 
# Tested on Ubuntu/WSL with example files

set -euo pipefail

clear
echo "══════════════════════════════════════════════════════════════════════"
echo "     PROTEIN-LIGAND DOCKING SCRIPT – InsilicoLAB (insilicolab82@gmail)"
echo "══════════════════════════════════════════════════════════════════════════════"

# -------------------------- INPUTS --------------------------
read -rp "Enter the path to the receptor file (PDB format, energy-minimized): " receptor_file
read -rp "Enter the path to the ligand file (PDB format, energy-minimized): " ligand_file
read -rp "Enter the path to the grid parameters file (txt format): " grid_file
read -rp "Enter the path to the reference ligand file (PDB/PDBQT): " ref_ligand_file

echo -e "\n══════════════════════════════════════════════════════════════════════\n"

# -------------------------- READ GRID FILE (ULTRA ROBUST) --------------------------
declare -A grid

while IFS='=' read -r rawkey rawvalue || [[ -n "$rawkey" ]]; do
    [[ -z "${rawkey// }" || "$rawkey" =~ ^[[:space:]]*# ]] && continue
    key=$(echo "$rawkey"  | xargs)     # trim spaces
    value=$(echo "$rawvalue" | xargs)
    case "$key" in
        center_x|center_y|center_z|size_x|size_y|size_z)
            grid["$key"]="$value"
            ;;
    esac
done < <(tr -d '\r' < "$grid_file")

# Check all 6 parameters exist
for k in center_x center_y center_z size_x size_y size_z; do
    [[ -z "${grid[$k]+x}" ]] && {
        echo "ERROR: Missing grid parameter '$k' in $grid_file"
        cat -n "$grid_file"
        exit 1
    }
done

echo "Grid loaded → Center (${grid[center_x]}, ${grid[center_y]}, ${grid[center_z]}) | Size (${grid[size_x]}, ${grid[size_y]}, ${grid[size_z]})"

# -------------------------- DOCKING PARAMETERS --------------------------
read -rp "Exhaustiveness [8-32] (default 32): " exhaustiveness; exhaustiveness=${exhaustiveness:-32}
read -rp "Scoring function [vina/vinardo/ad4] (default vina): " scoring; scoring=${scoring:-vina}
read -rp "Number of poses (default 20): " num_modes; num_modes=${num_modes:-20}
read -rp "Energy range (kcal/mol) (default 3.0): " energy_range; energy_range=${energy_range:-3.0}
read -rp "Grid spacing in Å (default 0.375): " grid_spacing; grid_spacing=${grid_spacing:-0.375}
read -rp "Show only poses with RMSD < X Å (999 = all, use any value to set cutoff ,ex: 3.0): " max_rmsd_show; max_rmsd_show=${max_rmsd_show:-999}

echo -e "\nStarting AutoDock Vina (exhaustiveness = $exhaustiveness, scoring = ${scoring^^}) ...\n"

# -------------------------- PREPARE PDBQT --------------------------
obabel "$receptor_file"   -O receptor.pdbqt -xr -p 7.4 --partialcharge gasteiger -h >/dev/null 2>&1
obabel "$ligand_file"     -O ligand.pdbqt   --partialcharge gasteiger -h >/dev/null 2>&1
obabel "$ref_ligand_file" -O ref_ligand.pdbqt --partialcharge gasteiger -h >/dev/null 2>&1

# -------------------------- RUN VINA --------------------------
./vina_1.2.7_linux_x86_64 \
    --receptor receptor.pdbqt \
    --ligand   ligand.pdbqt \
    --center_x "${grid[center_x]}" \
    --center_y "${grid[center_y]}" \
    --center_z "${grid[center_z]}" \
    --size_x   "${grid[size_x]}" \
    --size_y   "${grid[size_y]}" \
    --size_z   "${grid[size_z]}" \
    --spacing  "$grid_spacing" \
    --exhaustiveness "$exhaustiveness" \
    --scoring  "$scoring" \
    --num_modes "$num_modes" \
    --energy_range "$energy_range" \
    --out result_docking.pdbqt \
    --verbosity 2

[[ ! -f result_docking.pdbqt ]] && { echo "Vina failed!"; exit 1; }

# -------------------------- RMSD CALCULATION --------------------------
echo -e "\nCalculating RMSDs using LS-align...\n"

mkdir -p .tmp_poses && cd .tmp_poses
csplit -sz -f pose_ ../result_docking.pdbqt '/^MODEL/' '{*}' >/dev/null 2>&1

affinities=()
rmsd_vals=()
n=1

for f in pose_*; do
    [[ ! -f "$f" ]] && continue
    cp "$f" "ligand_$(printf "%04d" $n).pdbqt"

    obabel "ligand_$(printf "%04d" $n).pdbqt" -O "ligand_$n.mol2" -h >/dev/null 2>&1
    ((n==1)) && obabel ../ref_ligand.pdbqt -O ref.mol2 -h >/dev/null 2>&1

    ../LSalign "ligand_$n.mol2" ref.mol2 -o "al$n" > "rmsd$n.txt" 2>/dev/null || true

    rmsd=$(awk 'NR==4{print $8}' "rmsd$n.txt" 2>/dev/null || echo 999)
    affinity=$(grep "^REMARK VINA RESULT" ../result_docking.pdbqt | sed -n "${n}p" | awk '{print $4}')

    if (( $(echo "$rmsd <= $max_rmsd_show" | bc -l) )); then
        affinities+=("$affinity")
        rmsd_vals+=("$rmsd")
    fi
    ((n++))
done

cd ..
rm -rf .tmp_poses

# -------------------------- PLOT --------------------------
if (( ${#affinities[@]} > 0 )); then
    python3 - <<PY
import matplotlib.pyplot as plt, numpy as np, seaborn as sns
aff = [$(IFS=, ; echo "${affinities[*]}")]
rms = [$(IFS=, ; echo "${rmsd_vals[*]}")]
scoring = "${scoring^^}"
exh = $exhaustiveness

fig, ax = plt.subplots(1,2,figsize=(15,6))
ax[0].scatter(rms, aff, c=aff, cmap='coolwarm_r', s=120, edgecolor='k')
ax[0].set_xlabel('RMSD (Å)'); ax[0].set_ylabel('Affinity (kcal/mol)')
ax[0].set_title(f'{scoring} – Exhaustiveness {exh}')
best = np.argmin(aff)
ax[0].scatter(rms[best], aff[best], c='red', s=500, marker='*',
              label=f'Best: {aff[best]:.2f} kcal/mol @ {rms[best]:.2f} Å')
ax[0].legend(); ax[0].grid(alpha=0.3)

sns.kdeplot(rms, ax=ax[1], fill=True, color='skyblue', linewidth=2)
ax[1].hist(rms, bins=20, alpha=0.6, color='lightgray', edgecolor='black')
ax[1].set_xlabel('RMSD (Å)'); ax[1].set_title('RMSD Distribution')

plt.figtext(0.5, 0.01, 'Generated by InsilicoLAB • insilicolab82@gmail.com', ha='center', color='gray')
plt.tight_layout()
plt.savefig("Combined_Vina_Results.png", dpi=300)
plt.show()

print(f"BEST POSE: {aff[best]:.2f} kcal/mol at {rms[best]:.2f} Å RMSD")
PY
fi

# -------------------------- FINISH --------------------------
mv result_docking.pdbqt result_docking_final.pdbqt

echo -e "\nALL DONE!"
echo "→ result_docking_final.pdbqt"
echo "→ Combined_Vina_Results.png"
echo "Happy docking!"

exit 0