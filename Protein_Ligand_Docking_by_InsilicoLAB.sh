#!/bin/bash
# Protein–Ligand Docking Pipeline – InsilicoLAB (insilicoLAB82@gmail.com)

set -euo pipefail
clear
echo "══════════════════════════════════════════════════════════════════════"
echo " PROTEIN-LIGAND DOCKING SCRIPT – InsilicoLAB"
echo "                  ★ AutoDock Vina ★"
echo "══════════════════════════════════════════════════════════════════════"

# -------------------------- INPUTS --------------------------
read -rp "Enter the path to the receptor file (PDB format, energy-minimized): " receptor_file
read -rp "Enter the path to the ligand file (PDB format, energy-minimized): " ligand_file
read -rp "Enter the path to the grid parameters file (txt format): " grid_file
read -rp "Enter the path to the reference ligand file (PDB/PDBQT): " ref_ligand_file
echo -e "\n══════════════════════════════════════════════════════════════════════\n"

# -------------------------- READ GRID FILE --------------------------
declare -A grid
while IFS='=' read -r rawkey rawvalue || [[ -n "$rawkey" ]]; do
    [[ -z "${rawkey// }" || "$rawkey" =~ ^[[:space:]]*# ]] && continue
    key=$(echo "$rawkey" | xargs)
    value=$(echo "$rawvalue" | xargs)
    [[ -v grid["$key"] ]] || grid["$key"]="$value"
done < <(tr -d '\r' < "$grid_file")

for k in center_x center_y center_z size_x size_y size_z; do
    [[ -z "${grid[$k]+x}" ]] && { echo "ERROR: Missing '$k' in $grid_file"; exit 1; }
done
echo "Grid → Center (${grid[center_x]}, ${grid[center_y]}, ${grid[center_z]}) | Size (${grid[size_x]}, ${grid[size_y]}, ${grid[size_z]})"

# -------------------------- DOCKING PARAMETERS --------------------------
read -rp "Exhaustiveness [8-128] (default 32): " exhaustiveness; exhaustiveness=${exhaustiveness:-32}
read -rp "Scoring function [vina/vinardo/ad4] (default vina): " scoring; scoring=${scoring:-vina}
read -rp "Number of poses (default 20): " num_modes; num_modes=${num_modes:-20}
read -rp "Energy range (kcal/mol) (default 3.0): " energy_range; energy_range=${energy_range:-3.0}
read -rp "Grid spacing in Å (default 0.375): " grid_spacing; grid_spacing=${grid_spacing:-0.375}
read -rp "Show only poses with RMSD < X Å (999 = all): " max_rmsd_show; max_rmsd_show=${max_rmsd_show:-999}

# -------------------------- PREPARE PDBQT --------------------------
obabel "$receptor_file" -O receptor.pdbqt -xr -p 7.4 --partialcharge gasteiger -h >/dev/null 2>&1
obabel "$ligand_file"     -O ligand.pdbqt    --partialcharge gasteiger -h >/dev/null 2>&1
obabel "$ref_ligand_file" -O ref_ligand.pdbqt --partialcharge gasteiger -h >/dev/null 2>&1

# -------------------------- RUN VINA --------------------------
./vina_1.2.7_linux_x86_64 \
    --receptor receptor.pdbqt --ligand ligand.pdbqt \
    --center_x "${grid[center_x]}" --center_y "${grid[center_y]}" --center_z "${grid[center_z]}" \
    --size_x "${grid[size_x]}" --size_y "${grid[size_y]}" --size_z "${grid[size_z]}" \
    --spacing "$grid_spacing" --exhaustiveness "$exhaustiveness" --scoring "$scoring" \
    --num_modes "$num_modes" --energy_range "$energy_range" \
    --out result_docking.pdbqt --verbosity 2

[[ ! -f result_docking.pdbqt ]] && { echo "Vina failed!"; exit 1; }

# -------------------------- EXTRACT POSES & CALCULATE RMSD/AFFINITY --------------------------
mkdir -p .tmp_poses Best_pose && cd .tmp_poses
csplit -sz -f pose_ ../result_docking.pdbqt '/^MODEL/' '{*}' >/dev/null 2>&1

affinities=()
rmsd_vals=()
pose_files=()
n=1
best_affinity=999.0
best_rmsd=999.0
best_pose_num=0
best_pose_file=""

for f in pose_*; do
    [[ ! -f "$f" ]] && continue
    pose_pdbqt="ligand_$(printf "%04d" $n).pdbqt"
    pose_mol2="ligand_$n.mol2"
    cp "$f" "$pose_pdbqt"
    obabel "$pose_pdbqt" -O "$pose_mol2" -h >/dev/null 2>&1

    ((n==1)) && obabel ../ref_ligand.pdbqt -O ref.mol2 -h >/dev/null 2>&1

    ../LSalign "$pose_mol2" ref.mol2 -o "al$n" > "rmsd$n.txt" 2>/dev/null || true
    rmsd=$(awk 'NR==4{print $8}' "rmsd$n.txt" 2>/dev/null || echo 999)

    affinity=$(grep "^REMARK VINA RESULT" ../result_docking.pdbqt | sed -n "${n}p" | awk '{print $4}')

    # Store for plot (respecting RMSD filter)
    if (( $(echo "$rmsd <= $max_rmsd_show" | bc -l) )); then
        affinities+=("$affinity")
        rmsd_vals+=("$rmsd")
        pose_files+=("$pose_mol2")
    fi

    # Track absolute best pose by affinity
    if (( $(echo "$affinity < $best_affinity" | bc -l) )); then
        best_affinity="$affinity"
        best_rmsd="$rmsd"
        best_pose_num=$n
        best_pose_file="$pose_mol2"
    fi
    ((n++))
done
cd ..

# -------------------------- SAVE BEST POSE --------------------------
if [[ -n "$best_pose_file" && -f ".tmp_poses/$best_pose_file" ]]; then
    best_filename="Best_Pose_$(printf "%.2f" "$best_affinity")kcal_$(printf "%.2f" "$best_rmsd")A_RMSD_Pose${best_pose_num}.mol2"
    cp ".tmp_poses/$best_pose_file" "Best_pose/$best_filename"
    echo -e "\nBEST POSE SAVED → Best_pose/$best_filename"
    echo "   Affinity: $best_affinity kcal/mol | RMSD: $best_rmsd Å (LS-align)"
fi

# -------------------------- PLOT (FIXED f-STRING) --------------------------
if (( ${#affinities[@]} > 0 )); then
    python3 - <<PY
import matplotlib.pyplot as plt, numpy as np, seaborn as sns

aff = [$(IFS=, ; echo "${affinities[*]}")]
rms = [$(IFS=, ; echo "${rmsd_vals[*]}")]
scoring = "${scoring^^}"
exh = $exhaustiveness
max_rmsd = float($max_rmsd_show)
n_poses = len(aff)
best_idx = np.argmin(aff)
best_aff = aff[best_idx]
best_rms = rms[best_idx]

fig, ax = plt.subplots(1, 2, figsize=(17, 7), dpi=120)

# Scatter
sc = ax[0].scatter(rms, aff, c=aff, cmap='viridis_r', s=350, edgecolor='white', linewidth=1.5, alpha=0.98, zorder=3)
ax[0].scatter(rms[best_idx], aff[best_idx], c='yellow', s=1000, marker='*',
              edgecolor='red', linewidth=2.5, zorder=10)

for i in range(n_poses):
    ax[0].text(rms[i], aff[i], f"P{i+1}", fontsize=10, fontweight='bold',
               color='black', ha='center', va='bottom', zorder=11)

ax[0].text(rms[best_idx] + 0.2, aff[best_idx],
           f"BEST POSE (Pose{best_idx+1})\n{best_aff:.2f} kcal/mol\n{best_rms:.2f} Å",
           fontsize=8, fontweight='bold', color='darkred', ha='left', va='center',
           bbox=dict(boxstyle="round,pad=0.5", facecolor='gold', alpha=0.9, edgecolor='red', linewidth=1.5))

ax[0].set_xlabel('RMSD from Reference (LS-align, Å)', fontsize=13)
ax[0].set_ylabel('Binding Affinity (kcal/mol)', fontsize=13)

# Fixed title – no nested braces!
title = f"{scoring} Scoring • Exhaustiveness = {exh} • {n_poses} Poses"
if max_rmsd < 999:
    title += f" (RMSD < {max_rmsd} Å)"
ax[0].set_title(title, fontsize=14, pad=20)

ax[0].invert_yaxis()
ax[0].grid(alpha=0.3, linestyle='--')
ax[0].legend(['Best Pose'], loc='upper right', markerscale=1.2, fontsize=8)

# RMSD distribution
sns.kdeplot(rms, ax=ax[1], fill=True, color='#4a90e2', linewidth=3, alpha=0.7)
ax[1].hist(rms, bins=min(20, n_poses), color='#a0d8ff', edgecolor='black', alpha=0.8)
ax[1].set_xlabel('RMSD (Å)', fontsize=13)
ax[1].set_ylabel('Density', fontsize=13)
ax[1].set_title('RMSD Distribution (LS-align)', fontsize=14)

plt.figtext(0.5, 0.02,
            'Generated by InsilicoLAB/',
            ha='center', fontsize=11, style='italic', color='gray')

plt.tight_layout(rect=[0, 0.04, 1, 0.96])
plt.savefig("Combined_Vina_Results_Final.png", dpi=600, bbox_inches='tight')
plt.show()

print(f"\nBEST POSE: Pose {best_idx+1} → {best_aff:.2f} kcal/mol | {best_rms:.2f} Å RMSD")
print("Plot saved → Combined_Vina_Results_Final.png")
PY
fi

# -------------------------- FINALIZE --------------------------
mv result_docking.pdbqt result_docking_final.pdbqt
rm -rf .tmp_poses

echo -e "\n══════════════════════════════════════════════════════════════════════"
echo "                         ALL DONE – SUCCESS!"
echo "══════════════════════════════════════════════════════════════════════"
echo "→ result_docking_final.pdbqt"
echo "→ Combined_Vina_Results_Final.png"
echo "→ Best pose saved in: Best_pose/$best_filename"
echo -e "\nHappy docking! See you in the next time!\n"
exit 0
