#!/bin/bash
# Protein–Ligand Docking Pipeline – InsilicoLAB (Insilicolab82@gmail.com)

set -euo pipefail
clear
echo "══════════════════════════════════════════════════════════════════════════════════════════════════"
echo " PROTEIN-LIGAND DOCKING SCRIPT – by InsilicoLAB"
echo " Powered by AutoDock Vina "
echo " Best pose rule: lowest RMSD first(LS-align) → best affinity on tie"
echo "══════════════════════════════════════════════════════════════════════════════════════════════════"

# -------------------------- INPUTS --------------------------
read -rp "Enter the path to the receptor file (PDB format, energy-minimized): " receptor_file
read -rp "Enter the path to the ligand file (PDB format, energy-minimized): " ligand_file
read -rp "Enter the path to the grid parameters file (txt format): " grid_file
read -rp "Enter the path to the reference ligand file (PDB/PDBQT, for RMSD): " ref_ligand_file
echo -e "\n══════════════════════════════════════════════════════════════════════════════════════════════════\n"

# -------------------------- READ GRID FILE --------------------------
declare -A grid
while IFS='=' read -r rawkey rawvalue || [[ -n "$rawkey" ]]; do
    [[ -z "${rawkey// }" || "$rawkey" =~ ^[[:space:]]*# ]] && continue
    key=$(echo "$rawkey" | xargs)
    value=$(echo "$rawvalue" | xargs)
    grid["$key"]="$value"
done < <(tr -d '\r' < "$grid_file")

for k in center_x center_y center_z size_x size_y size_z; do
    [[ -z "${grid[$k]+x}" ]] && { echo "ERROR: Missing '$k' in $grid_file"; exit 1; }
done
echo "Grid → Center (${grid[center_x]}, ${grid[center_y]}, ${grid[center_z]}) | Size (${grid[size_x]}, ${grid[size_y]}, ${grid[size_z]})"

# -------------------------- DOCKING PARAMETERS --------------------------
echo "══════════════════════════════════════════════════════════════════════════════════════════════════"
echo " PROTEIN-LIGAND DOCKING PARAMETERS"
read -rp "Exhaustiveness [8-32] (default 8): " exhaustiveness; exhaustiveness=${exhaustiveness:-8}
read -rp "Scoring function [vina/vinardo/ad4] (default vina): " scoring; scoring=${scoring:-vina}
read -rp "Number of poses (default 10): " num_modes; num_modes=${num_modes:-10}
read -rp "Energy range (kcal/mol) (default 3.0): " energy_range; energy_range=${energy_range:-3.0}
read -rp "Grid spacing in Å (default 0.375): " grid_spacing; grid_spacing=${grid_spacing:-0.375}
read -rp "Show only poses with RMSD < X Å (999 = all): " max_rmsd_show; max_rmsd_show=${max_rmsd_show:-999}
echo "══════════════════════════════════════════════════════════════════════════════════════════════════"
# -------------------------- PREPARE PDBQT --------------------------
obabel "$receptor_file" -O receptor.pdbqt -xr -p 7.4 --partialcharge gasteiger -h >/dev/null 2>&1
obabel "$ligand_file" -O ligand.pdbqt --partialcharge gasteiger -h >/dev/null 2>&1
obabel "$ref_ligand_file" -O ref_ligand.pdbqt --partialcharge gasteiger -h >/dev/null 2>&1

# -------------------------- RUN VINA (PROPER LINE CONTINUATION) --------------------------
./vina_1.2.7_linux_x86_64 \
    --receptor receptor.pdbqt \
    --ligand ligand.pdbqt \
    --center_x "${grid[center_x]}" \
    --center_y "${grid[center_y]}" \
    --center_z "${grid[center_z]}" \
    --size_x "${grid[size_x]}" \
    --size_y "${grid[size_y]}" \
    --size_z "${grid[size_z]}" \
    --spacing "$grid_spacing" \
    --exhaustiveness "$exhaustiveness" \
    --scoring "$scoring" \
    --num_modes "$num_modes" \
    --energy_range "$energy_range" \
    --out result_docking.pdbqt \
    --verbosity 2

[[ ! -f result_docking.pdbqt ]] && { echo "Vina failed!"; exit 1; }

# -------------------------- EXTRACT POSES & CALCULATE RMSD/AFFINITY --------------------------
mkdir -p .tmp_poses Best_pose && cd .tmp_poses
csplit -sz -f pose_ ../result_docking.pdbqt '/^MODEL/' '{*}' >/dev/null 2>&1

affinities=()
rmsd_vals=()
pose_files=()
n=1

best_rmsd=999.0
best_affinity=999.0
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

    # NEW BEST POSE LOGIC: lowest RMSD wins → tie = best affinity
    if (( $(echo "$rmsd < $best_rmsd - 0.001" | bc -l) )); then
        best_rmsd="$rmsd"
        best_affinity="$affinity"
        best_pose_num=$n
        best_pose_file="$pose_mol2"
    elif (( $(echo "$rmsd >= $best_rmsd - 0.001 && $rmsd <= $best_rmsd + 0.001" | bc -l) )); then
        if (( $(echo "$affinity < $best_affinity" | bc -l) )); then
            best_affinity="$affinity"
            best_pose_num=$n
            best_pose_file="$pose_mol2"
        fi
    fi

    # Store for plotting (respect user RMSD filter)
    if (( $(echo "$rmsd <= $max_rmsd_show" | bc -l) )); then
        affinities+=("$affinity")
        rmsd_vals+=("$rmsd")
        pose_files+=("$pose_mol2")
    fi
    ((n++))
done
cd ..

# -------------------------- SAVE BEST POSE --------------------------
if [[ -n "$best_pose_file" && -f ".tmp_poses/$best_pose_file" ]]; then
    best_filename="Best_Pose_$(printf "%.3f" "$best_affinity")kcal_$(printf "%.3f" "$best_rmsd")A_RMSD_Pose${best_pose_num}.mol2"
    cp ".tmp_poses/$best_pose_file" "Best_pose/$best_filename"
    echo -e "\nBEST POSE SAVED → Best_pose/$best_filename"
    echo " Affinity: $best_affinity kcal/mol | RMSD: $best_rmsd Å (LS-align)"
    echo " Rule: lowest RMSD first → best affinity on tie"
fi

# -------------------------- 3×2 OPTIMIZED LAYOUT --------------------------
if (( ${#affinities[@]} > 0 )); then
    python3 - <<PY
import matplotlib.pyplot as plt
import numpy as np
import seaborn as sns
import sys
from matplotlib import colormaps
from scipy.stats import gaussian_kde

try:
    # Safely parse the arrays - handle empty or malformed data
    aff_raw = "$(printf "%s," "${affinities[@]}" | sed 's/,$//')"
    rms_raw = "$(printf "%s," "${rmsd_vals[@]}" | sed 's/,$//')"
    
    # Convert to numpy arrays, handling empty cases
    if aff_raw and rms_raw:
        aff = np.array([float(x) for x in aff_raw.split(',') if x.strip()])
        rms = np.array([float(x) for x in rms_raw.split(',') if x.strip()])
    else:
        print("ERROR: No affinity or RMSD data available", file=sys.stderr)
        sys.exit(1)
    
    n_poses = len(aff)
    
    if n_poses == 0:
        print("ERROR: No valid poses found", file=sys.stderr)
        sys.exit(1)
        
    best_pose_num = $best_pose_num
    best_aff = float("$best_affinity")
    best_rms = float("$best_rmsd")
    best_idx = best_pose_num - 1 if best_pose_num <= n_poses else 0
    
    # Adjust best_idx if out of bounds
    if best_idx >= n_poses:
        best_idx = 0
        best_aff = aff[0]
        best_rms = rms[0]
    
    mean_aff = np.mean(aff); std_aff = np.std(aff)
    mean_rms = np.mean(rms); std_rms = np.std(rms)
    
    # Optimized style for compact layout
    plt.rcParams.update({
        'font.size': 10,
        'font.family': 'Arial',
        'axes.titlesize': 12,
        'axes.labelsize': 11,
        'xtick.labelsize': 9,
        'ytick.labelsize': 9,
        'legend.fontsize': 9,
        'axes.linewidth': 1.0,
        'axes.titlepad': 10,
        'axes.labelpad': 6,
        'figure.constrained_layout.use': False,  # Disable for manual control
        'figure.titlesize': 16,
        'figure.titleweight': 'bold'
    })
    sns.set_style("whitegrid")
    
    # Create compact figure with better proportions
    fig = plt.figure(figsize=(19, 14), dpi=600)  # Slightly larger for better spacing
    
    # Create main title at the top with proper spacing
    fig.suptitle('AutoDock Vina Docking Results', fontsize=18, 
                 fontweight='bold', y=0.98)
    
    # Adjust grid spec with proper spacing - fixed values to prevent zero-size axes
    gs = fig.add_gridspec(3, 2, width_ratios=[1.6, 1], height_ratios=[1, 1, 1.1],
                         hspace=0.35, wspace=0.25,
                         top=0.92, bottom=0.07, left=0.08, right=0.95)
    
    # 1. Affinity vs RMSD - WITH MODEL NUMBERS
    ax1 = fig.add_subplot(gs[0, 0])
    scatter = ax1.scatter(rms, aff, c=aff, cmap='viridis_r', s=120, 
                          edgecolors='white', linewidth=0.8, zorder=5, alpha=0.8)
    ax1.scatter(best_rms, best_aff, facecolors='gold', edgecolor='red', 
                s=250, marker='*', linewidth=1.8, zorder=10, label='Best Pose')
    
    # Calculate dynamic offsets to prevent label overlap
    rms_range = np.ptp(rms) if len(rms) > 1 else 1.0
    aff_range = np.ptp(aff) if len(aff) > 1 else 1.0
    
    # Label only key poses to avoid clutter
    key_indices = set([best_idx, 0, n_poses-1])
    # Add a few more key poses based on affinity extremes
    if n_poses > 5:
        key_indices.add(np.argmin(aff))
        key_indices.add(np.argmax(aff))
        key_indices.add(np.argmin(rms))
        key_indices.add(np.argmax(rms))
    
    for i in sorted(key_indices):
        if i >= n_poses:
            continue
            
        label_text = f"MODEL {i+1}"
        
        # Calculate offset based on data position to avoid overlap
        if rms_range > 0 and aff_range > 0:
            x_norm = (rms[i] - rms.min()) / rms_range
            y_norm = (aff[i] - aff.min()) / aff_range
        else:
            x_norm = 0.5
            y_norm = 0.5
        
        # Set offsets based on quadrant
        offset_factor = 0.05
        if x_norm < 0.5 and y_norm < 0.5:
            x_offset = -offset_factor * rms_range
            y_offset = -offset_factor * aff_range
            ha = 'right'
            va = 'top'
        elif x_norm >= 0.5 and y_norm < 0.5:
            x_offset = offset_factor * rms_range
            y_offset = -offset_factor * aff_range
            ha = 'left'
            va = 'top'
        elif x_norm < 0.5 and y_norm >= 0.5:
            x_offset = -offset_factor * rms_range
            y_offset = offset_factor * aff_range
            ha = 'right'
            va = 'bottom'
        else:
            x_offset = offset_factor * rms_range
            y_offset = offset_factor * aff_range
            ha = 'left'
            va = 'bottom'
        
        # Special formatting for best pose
        if i == best_idx:
            bbox_props = dict(boxstyle="round,pad=0.15", facecolor="yellow", 
                             edgecolor="red", linewidth=1.0, alpha=0.9)
            ax1.text(rms[i] + x_offset, aff[i] + y_offset, label_text,
                    fontsize=8, fontweight='bold', color='darkred',
                    ha=ha, va=va, bbox=bbox_props)
        else:
            ax1.text(rms[i] + x_offset, aff[i] + y_offset, label_text,
                    fontsize=7, fontweight='normal', color='black',
                    ha=ha, va=va, alpha=0.7)
    
    ax1.set_xlabel('RMSD from Reference (LS-align, Å)', labelpad=8)
    ax1.set_ylabel('Binding Affinity (kcal/mol)', labelpad=8)
    ax1.set_title('Affinity vs RMSD', fontweight='bold', pad=12)
    ax1.invert_yaxis()
    ax1.grid(alpha=0.25, linestyle='--', linewidth=0.5)
    ax1.legend(loc='upper right', framealpha=0.9, fontsize=8)
    
    # Add compact colorbar
    cbar = plt.colorbar(scatter, ax=ax1, pad=0.01, shrink=0.8)
    cbar.set_label('Affinity', rotation=270, labelpad=10, fontsize=9)
    cbar.ax.tick_params(labelsize=8)
    
    # 2. Affinity Distribution with gradient colors
    ax2 = fig.add_subplot(gs[0, 1])
    bins = max(5, min(12, n_poses // 2))
    
    # Create gradient colormap for histogram - FIXED DEPRECATION
    blues_cmap = colormaps['Blues']
    hist_counts, hist_bins, patches = ax2.hist(aff, bins=bins, edgecolor='black', 
                                                linewidth=0.7, alpha=0.85)
    
    # Color bars with gradient
    if len(hist_counts) > 0:
        bin_centers = 0.5 * (hist_bins[:-1] + hist_bins[1:])
        if len(bin_centers) > 0:
            col = (bin_centers - min(bin_centers))
            if max(col) > 0:
                col = col / max(col)
            else:
                col = np.ones_like(bin_centers)
            
            for c, p in zip(col, patches):
                plt.setp(p, 'facecolor', blues_cmap(c))
    
    # Add KDE curve
    if n_poses > 1:
        kde = gaussian_kde(aff)
        x_range = np.linspace(aff.min(), aff.max(), 200)
        ax2.plot(x_range, kde(x_range) * len(aff) * (hist_bins[1] - hist_bins[0]), 
                 color='darkblue', linewidth=2.0, label='KDE', alpha=0.9)
    
    ax2.axvline(best_aff, color='red', linewidth=2.0, linestyle='--', 
                label=f'Best: {best_aff:.2f}')
    ax2.set_xlabel('Affinity (kcal/mol)', labelpad=8)
    ax2.set_ylabel('Frequency', labelpad=8)
    ax2.set_title('Affinity Distribution', fontweight='bold', pad=12)
    ax2.legend(loc='upper right', fontsize=8)
    ax2.grid(alpha=0.2, linestyle=':', axis='y')
    
    # 3. Pose Ranking by Affinity with gradient colors
    ax3 = fig.add_subplot(gs[1, 0])
    sorted_indices = np.argsort(aff)
    sorted_aff = aff[sorted_indices]
    
    # Create gradient colors for bars - FIXED DEPRECATION
    reds_cmap = colormaps['Reds']
    bar_norm = plt.Normalize(sorted_aff.min(), sorted_aff.max())
    
    bars = ax3.bar(range(1, n_poses+1), sorted_aff, edgecolor='black', 
                   linewidth=0.5, alpha=0.9, width=0.6)
    
    # Apply gradient colors
    for i, (bar, val) in enumerate(zip(bars, sorted_aff)):
        color = reds_cmap(bar_norm(val))
        bar.set_facecolor(color)
    
    # Highlight best pose with special style
    if best_idx in sorted_indices:
        best_bar_idx = np.where(sorted_indices == best_idx)[0][0]
        if best_bar_idx < len(bars):
            bars[best_bar_idx].set_facecolor('gold')
            bars[best_bar_idx].set_edgecolor('red')
            bars[best_bar_idx].set_linewidth(1.5)
            bars[best_bar_idx].set_alpha(1.0)
            
            # Add MODEL label on the best bar
            ax3.text(best_bar_idx + 1, sorted_aff[best_bar_idx] + 0.02, f"MODEL {best_pose_num}",
                    fontsize=8, fontweight='bold', color='red',
                    ha='center', va='bottom')
    
    ax3.set_xlabel('Pose Rank (best to worst)', labelpad=8)
    ax3.set_ylabel('Affinity (kcal/mol)', labelpad=8)
    ax3.set_title('Pose Ranking by Affinity', fontweight='bold', pad=12)
    ax3.invert_yaxis()
    # Set x-ticks
    if n_poses > 10:
        ax3.set_xticks(range(1, n_poses+1, max(1, n_poses//8)))
    else:
        ax3.set_xticks(range(1, n_poses+1))
    ax3.grid(axis='y', alpha=0.2, linestyle=':')
    
    # 4. RMSD Distribution with gradient colors
    ax4 = fig.add_subplot(gs[1, 1])
    bins = max(5, min(12, n_poses // 2))
    
    # Create gradient colormap for histogram - FIXED DEPRECATION
    greens_cmap = colormaps['Greens']
    hist_counts_rms, hist_bins_rms, patches_rms = ax4.hist(rms, bins=bins, 
                                                           edgecolor='black', 
                                                           linewidth=0.7, alpha=0.85)
    
    # Color bars with gradient
    if len(hist_counts_rms) > 0:
        bin_centers_rms = 0.5 * (hist_bins_rms[:-1] + hist_bins_rms[1:])
        if len(bin_centers_rms) > 0:
            col_rms = (bin_centers_rms - min(bin_centers_rms))
            if max(col_rms) > 0:
                col_rms = col_rms / max(col_rms)
            else:
                col_rms = np.ones_like(bin_centers_rms)
            
            for c, p in zip(col_rms, patches_rms):
                plt.setp(p, 'facecolor', greens_cmap(c))
    
    # Add KDE curve
    if n_poses > 1:
        kde_rms = gaussian_kde(rms)
        x_range_rms = np.linspace(rms.min(), rms.max(), 200)
        ax4.plot(x_range_rms, kde_rms(x_range_rms) * len(rms) * (hist_bins_rms[1] - hist_bins_rms[0]), 
                 color='darkgreen', linewidth=2.0, label='KDE', alpha=0.9)
    
    ax4.axvline(best_rms, color='red', linewidth=2.0, linestyle='--',
                label=f'Best: {best_rms:.2f} Å')
    ax4.set_xlabel('RMSD (Å)', labelpad=8)
    ax4.set_ylabel('Frequency', labelpad=8)
    ax4.set_title('RMSD Distribution', fontweight='bold', pad=12)
    ax4.legend(loc='upper right', fontsize=8)
    ax4.grid(alpha=0.2, linestyle=':', axis='y')
    
    # 5. Affinity by Pose Index
    ax5 = fig.add_subplot(gs[2, 0])
    ax5.plot(range(1, n_poses+1), aff, color='gray', linewidth=1.0, 
             alpha=0.6, marker='o', markersize=4, zorder=1)
    
    # Color points with gradient based on affinity - FIXED DEPRECATION
    coolwarm_cmap = colormaps['coolwarm']
    point_norm = plt.Normalize(aff.min(), aff.max())
    
    for i in range(n_poses):
        color = coolwarm_cmap(point_norm(aff[i]))
        edgecolor = 'red' if i == best_idx else 'white'
        ax5.scatter(i+1, aff[i], c=[color], s=50, edgecolor=edgecolor, 
                   linewidth=1.0, zorder=2)
    
    # Highlight best pose
    ax5.scatter(best_pose_num, best_aff, c='gold', s=120, marker='*', 
                edgecolor='red', linewidth=1.2, zorder=10, label='Best Pose')
    
    # Add MODEL label for best pose only to avoid clutter
    ax5.text(best_pose_num, best_aff + 0.1, f"MODEL {best_pose_num}",
            fontsize=8, fontweight='bold', color='red',
            ha='center', va='bottom', alpha=0.9)
    
    ax5.set_xlabel('Pose Index (Vina order)', labelpad=8)
    ax5.set_ylabel('Affinity (kcal/mol)', labelpad=8)
    ax5.set_title('Affinity by Pose Index', fontweight='bold', pad=12)
    ax5.invert_yaxis()
    # Set x-ticks
    if n_poses > 10:
        ax5.set_xticks(range(1, n_poses+1, max(1, n_poses//8)))
    else:
        ax5.set_xticks(range(1, n_poses+1))
    ax5.grid(alpha=0.25, linestyle='--', linewidth=0.5)
    ax5.legend(loc='upper right', fontsize=8)
    
    # 6. DOCKING SUMMARY
    ax6 = fig.add_subplot(gs[2, 1])
    ax6.axis('off')
    
    # Get scoring parameters from shell
    scoring = "${scoring:-vina}"
    exhaustiveness = "${exhaustiveness:-8}"
    grid_spacing = "${grid_spacing:-0.375}"
    energy_range = "${energy_range:-3}"
    
    summary = f"""
                  DOCKING SUMMARY
──────────────────────────────────────────
Parameters
  • Scoring        : {scoring.upper()}
  • Exhaustiveness : {exhaustiveness}
  • Grid spacing   : {grid_spacing} Å
  • Energy range   : {energy_range} kcal/mol
  • Poses analyzed : {n_poses}

Results
  • Best affinity  : {best_aff:.3f} kcal/mol
  • Avg affinity   : {mean_aff:.3f} ± {std_aff:.2f} kcal/mol
  • Best RMSD      : {best_rms:.3f} Å
  • Avg RMSD       : {mean_rms:.2f} ± {std_rms:.2f} Å

Best Pose
  • MODEL number   : {best_pose_num}
  • Affinity       : {best_aff:.3f} kcal/mol
  • RMSD           : {best_rms:.3f} Å

Selection Rule
  Lowest RMSD first → best affinity on tie
──────────────────────────────────────────
In silico LAB – insilicolab82@gmail.com
"""
    ax6.text(0.05, 0.98, summary, transform=ax6.transAxes, fontsize=9.5,
             va='top', fontfamily='monospace', linespacing=1.5,
             bbox=dict(boxstyle="round,pad=0.6", facecolor='#f8f9fa', 
                      edgecolor='#2c3e50', linewidth=1.2, alpha=0.95))
    
    # Use tight_layout instead of constrained_layout for better control
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    
    # Save
    plt.savefig("AutoDock_Vina_Results_combine.png", dpi=600, 
                bbox_inches='tight', facecolor='white', pad_inches=0.2)
    plt.savefig("AutoDock_Vina_Results_combine.pdf", dpi=600, 
                bbox_inches='tight', pad_inches=0.2)
    
    print("\nOPTIMIZED 3×2 LAYOUT READY!")
    print("   → AutoDock_Vina_Results_combine.png  (600 DPI)")
    print("   → AutoDock_Vina_Results_combine.pdf   (vector)")
    print(f"   Best pose: MODEL {best_pose_num} | {best_aff:.3f} kcal/mol | {best_rms:.3f} Å")

except Exception as e:
    print(f"ERROR in plotting: {str(e)}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)
    print("Affinities:", aff_raw, file=sys.stderr)
    print("RMSDs:", rms_raw, file=sys.stderr)
    sys.exit(1)
PY
fi
# -------------------------- CLEANUP --------------------------
mv result_docking.pdbqt result_docking_final.pdbqt
rm -rf .tmp_poses receptor.pdbqt ligand.pdbqt ref_ligand.pdbqt
echo -e "\n══════════════════════════════════════════════════════════════════════"
echo " ALL DONE – SUCCESS!"
echo " → result_docking_final.pdbqt"
echo " → Combined_Vina_Results_Final.png"
echo " → Best pose in: Best_pose/$best_filename"
echo -e "\nHappy docking! Your best pose is the most native-like one.\n"
exit 0
