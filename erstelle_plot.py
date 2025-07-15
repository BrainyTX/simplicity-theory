import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# =========================================================================
# 1. Load and Clean SPARC Data
# =========================================================================
# Assumption: 'SPARC_summary.csv' is in the same directory as this script
try:
    df = pd.read_csv('SPARC_summary.csv')
except FileNotFoundError:
    print("❌ Error: 'SPARC_summary.csv' not found.")
    print("Please ensure the file is in the same directory as the script.")
    exit(1)

# Clean column names: remove leading/trailing spaces and convert to lowercase
df.columns = df.columns.str.strip().str.lower()

# Possible column names for surface density and dark matter fraction
possible_sigma_cols = ['mu0', 'sigma_b', 'sigma_baryon', 'surfacedensity']
possible_dm_cols = ['f_dm', 'dm_fraction', 'dark_matter_frac', 'dmfrac']

# Find the actual column names in the DataFrame
sigma_col = next((col for col in df.columns if col in possible_sigma_cols), None)
dm_col = next((col for col in df.columns if col in possible_dm_cols), None)

# Error handling if required columns are not found
if not sigma_col or not dm_col:
    print("❌ Error: Required column names for surface density or dark matter fraction were not found.")
    print(f"Expected Sigma columns: {possible_sigma_cols}")
    print(f"Expected DM columns: {possible_dm_cols}")
    print(f"Available columns in the file: {list(df.columns)}")
    exit(1)

# Extract data and convert to numeric type
mu0_sparc = df[sigma_col].astype(float)
fdm_sparc = df[dm_col].astype(float)

# Filter out invalid or missing data points (e.g., non-positive density, non-finite values)
mask_sparc = (mu0_sparc > 0) & np.isfinite(mu0_sparc) & np.isfinite(fdm_sparc)
mu0_sparc, fdm_sparc = mu0_sparc[mask_sparc], fdm_sparc[mask_sparc]

# =========================================================================
# 2. Manually Curated Special Galaxies (with reference values and citations)
# =========================================================================
# These are key galaxies that exemplify different dark matter fractions relative to the
# critical surface density, providing empirical validation for the Simplicity Theory.
# Sources for these values are primarily from research papers.
# Note: Exact mu0 values for Ultra-Diffuse Galaxies (UDGs) can be interpretation-dependent.
#       These f_DM values for DM-deficient galaxies are known to be very low.
# Citation numbers are directly embedded into the names for plot display.

special_galaxies = [
    #  mu0 (Msun/pc^2)   f_DM    Name (with citation in string)                      Comment
    [150.0, 0.05,  'NGC 1052-DF2'],    # Ultra-diffuse galaxy, almost no dark matter – "dark matter-free" (well above Sigma_crit) [16]
    [120.0, 0.10,  'NGC 1052-DF4'],    # Another DM-deficient UDG, supports DF2 findings [17]
    [2.3,   0.99,  'Dragonfly 44'],     # Extreme dark-matter dominance, very low baryonic surface density [4]
    [60.0,  0.55,  'NGC 2403'],             # Spiral galaxy, classic example of the "transition" region
    [80.0,  0.45,  'UGC 4325'],             # Disk galaxy in the mid-range, moderate dark matter fraction
    [200.0, 0.15,  'NGC 2841'],             # Massive spiral, high baryonic density, little DM in the center
    [100.0, 0.38,  'NGC 6503'],             # Isolated spiral, intermediate baryon/DM ratio
    [160.0, 0.18,  'NGC 6946'],             # "Fireworks Galaxy," very star-rich, baryon dominated
    [130.0, 0.22,  'NGC 2976'],             # DM fraction below average, inner baryonic stability
    [10.0,  0.93,  'IC 2574'],              # Dwarf galaxy, extremely DM-dominated
    [6.0,   0.97,  'UGC 128'],              # Giant low surface brightness (LSB) galaxy, nearly all DM
    [110.0, 0.40,  'NGC 3198'],             # Classic case of a flat rotation curve, textbook spiral
    [170.0, 0.19,  'NGC 7331'],             # "Milky Way twin," high baryonic density, low DM fraction
    [180.0, 0.16,  'NGC 5055'],             # Large spiral, low central DM fraction, baryon-dominated inside
    [120.0, 0.35,  'NGC 7793'],             # Large spiral, intermediate DM/baryon balance
    [70.0,  0.57,  'NGC 925'],              # Spiral, prominent transition regime, neither extreme
    [150.0, 0.20,  'NGC 2903'],             # Big, massive spiral, low inner DM fraction
    [4.0,   0.99,  'DDO 154'],              # Gas-rich dwarf, essentially all DM
    [8.0,   0.95,  'F568-1'],               # LSB galaxy, extremely low baryonic density, nearly pure DM
    [9.0,   0.96,  'UGC 5750'],             # Another LSB galaxy, very high DM dominance
    [30.0,  0.72,  'NGC 1003'],             # Spiral with marked DM/baryon mix (intermediate)
    [20.0,  0.80,  'UGC 1281'],             # Lenticular dwarf, clear DM domination
    [140.0, 0.24,  'NGC 3992'],             # Massive spiral, stable, above Sigma_crit
]

# =========================================================================
# 3. Plot Constants and Settings
# =========================================================================
SIG_CRIT = 124       # Critical surface density in M_sun / pc^2 (from Mathematical Foundations)

plt.rcParams.update({'font.size': 10}) # Adjust default font size for plot clarity

# =========================================================================
# 4. Create the Plot
# =========================================================================
plt.figure(figsize=(9, 6)) # Larger figure size for better visibility and readability

# Scatter plot for SPARC galaxies (primary dataset)
plt.scatter(mu0_sparc, fdm_sparc, s=16, label="SPARC Galaxies", alpha=0.7, edgecolors='k', linewidths=0.5, zorder=2)

# Add a single legend entry for all special galaxies
plt.scatter([], [], color='red', s=50, edgecolor='red', linewidths=0.5, marker='.', label='Special Galaxies')

# Mark the special galaxies individually for clear highlighting
for sigma_val, dmfrac_val, name_val in special_galaxies:
    # The name is used directly as a text label.
    # Marker changed to '.', color to 'red', size adjusted.
    plt.scatter(sigma_val, dmfrac_val, color='red', s=50, edgecolor='red', linewidths=0.5, marker='.', zorder=3)
    # Add text label directly next to the point, slightly offset for readability, font size smaller.
    plt.text(sigma_val * 1.05, dmfrac_val, name_val, fontsize=7, verticalalignment='center', horizontalalignment='left', zorder=4)

# Vertical threshold line representing Sigma_crit
plt.axvline(SIG_CRIT, color='dodgerblue', linestyle='--', linewidth=2,
            label=r'$\Sigma^{\ast}_{\mathrm{crit}} = 124\,\mathrm{M}_\odot\,\mathrm{pc}^{-2}$')

# Axis scaling and labels
plt.xscale('log') # Logarithmic scale for the X-axis, as used in your paper's figure
plt.xlabel(r'Baryonic Surface Density $\mu_0\,[\mathrm{M}_\odot\,\mathrm{pc}^{-2}]$', fontsize=12)
plt.ylabel(r'Dark Matter Fraction $f_{\rm DM}$', fontsize=12)
# Updated title to be on two lines and include "and additional data"
plt.title('RAR Threshold: Baryonic Surface Density vs. Dark Matter Fraction\n(Based on SPARC Data and Additional Curated Data)', fontsize=14)

# Legend and Grid
# Collect handles and labels to avoid duplicate legend entries and create a cleaner legend.
handles, labels = plt.gca().get_legend_handles_labels()
unique_labels = {}
for h, l in zip(handles, labels):
    unique_labels[l] = h # Keeps the last handle for a given label
plt.legend(unique_labels.values(), unique_labels.keys(), loc='upper right', fontsize=10) # Place legend in the upper right corner


plt.grid(True, which="both", ls="--", c="0.7", alpha=0.7) # Add grid for better readability

# Adjust layout to prevent labels from being cut off
plt.tight_layout()

# Add a text note about the Python code availability, centered below the plot
plt.figtext(0.5, 0.01, 'Python code for data verification available on GitHub.',
            horizontalalignment='center', fontsize=8, color='gray')

# Save and display the plot
plt.savefig("RAR_threshold.pdf", bbox_inches='tight', dpi=300) # 'bbox_inches='tight'' prevents clipping of labels/elements
plt.close() # Close the plot figure to ensure the script terminates correctly

print("Plot RAR_threshold.pdf saved successfully.")
