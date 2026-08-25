import numpy as np
import matplotlib.pyplot as plt

# Set font settings suitable for publication/EPS export (Type 42 vector fonts)
plt.rcParams['font.sans-serif'] = 'DejaVu Sans'
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['ps.fonttype'] = 42
plt.rcParams['pdf.fonttype'] = 42

# ==========================================
# 1. SIMULATION PARAMETERS
# ==========================================
N_events = 10000
np.random.seed(42)

# Unshifted reference integral value I_1 (Arbitrary scale)
# For this, just used the decay fit amplitude from the data
#I1_true = 574.815 # alphas
I1_true = 800.556 # protons

# Centroids and standard deviations (Spread)
# The lambda_mean comes from an exponential fit of the data,
# and the other values are estimated from a visual inspection
# of the data.

# For alphas
#dx_mean     = 100.0       # Center time shift (ns)
#dx_std      = 100.0       # Spread in time shifts (ns)
#lambda_mean = 0.000407085 # Mean decay constant (1/ns)
#lambda_std  = 0.0001      # Spread in decay constants across particles (1/ns)

# For protons
dx_mean     = 0.0         # Center time shift (ns)
dx_std      = 100.0       # Spread in time shifts (ns)
lambda_mean = 0.00031785  # Mean decay constant (1/ns)
lambda_std  = 0.0001      # Spread in decay constants across particles (1/ns)

# Generate distributions of true event parameters
dx_samples = np.random.normal(dx_mean, dx_std, N_events)
lambda_true_samples = np.random.normal(lambda_mean, lambda_std, N_events)

# ==========================================
# 2. INTEGRAL GENERATION AND CORRECTION
# ==========================================
# Measured shifted integral I_2 = I_1 * exp(-lambda_true * dx)
I2_measured = I1_true * np.exp(-lambda_true_samples * dx_samples)

# Ideal correction (using exact per-event lambda_true and dx)
I1_ideal_corrected = I2_measured * np.exp(lambda_true_samples * dx_samples)

# Approximate correction (using fixed mean lambda_bar and exact dx from TDC)
I1_approx_corrected = I2_measured * np.exp(lambda_mean * dx_samples)

# Residual error when using mean lambda instead of true lambda
residual_error = I1_approx_corrected - I1_ideal_corrected

# Percentage error relative to true integral
percentage_error = (residual_error / I1_true) * 100.0

# Uncorrected error for baseline comparison
uncorrected_error = I2_measured - I1_true

# Print summary metrics for the fractional correction
mean_pct_error = np.mean(percentage_error)
abs_mean_pct_error = np.mean(np.abs(percentage_error))

print(f"Mean Percent Error:          {mean_pct_error:+.3f}%")
print(f"Mean Absolute Percent Error: {abs_mean_pct_error:.3f}%")

# ==========================================
# 3. PLOTTING RESIDUAL ERRORS (Opaque colors for clean EPS rendering)
# ==========================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# Plot 1: Histogram of Uncorrected vs. Corrected Integrals
axes[0].hist(uncorrected_error, bins=60, alpha=0.5, color='red', label='Uncorrected ($I_2 - I_1$)')
axes[0].hist(residual_error, bins=60, alpha=0.7, color='blue', label=r'Corrected using $\bar{\lambda}$')
axes[0].axvline(0, color='black', linestyle='--', alpha=0.7)
axes[0].set_title("Integral Error Distribution")
axes[0].set_xlabel(r"Integral Offset ($I_{\mathrm{corr}} - I_1$) [a.u.]")
axes[0].set_ylabel("Counts")
axes[0].legend()
axes[0].grid(True, linestyle=':', alpha=0.6)

# Plot 2: Residual Error vs. Time Shift (dx)
sc1 = axes[1].scatter(dx_samples, percentage_error, c=lambda_true_samples, cmap='viridis', s=6, edgecolors='none')
axes[1].axhline(0, color='black', linestyle='--')
axes[1].set_title(r"Error vs. Time Shift ($\Delta x$)")
axes[1].set_xlabel(r"Time Shift $\Delta x$ (ns)")
axes[1].set_ylabel("Correction Error (%)")
cbar1 = fig.colorbar(sc1, ax=axes[1])
cbar1.set_label(r"True $\lambda$ ($\mathrm{ns}^{-1}$)")
axes[1].grid(True, linestyle=':', alpha=0.6)

# Plot 3: Residual Error vs. True Lambda
sc2 = axes[2].scatter(lambda_true_samples, percentage_error, c=np.abs(dx_samples), cmap='plasma', s=6, edgecolors='none')
axes[2].axvline(lambda_mean, color='#e41a1c', linestyle='--', label=r'Average $\bar{\lambda}$')
axes[2].axhline(0, color='black', linestyle='--')
axes[2].set_title(r"Error vs. True Decay Constant ($\lambda$)")
axes[2].set_xlabel(r"True $\lambda$ ($\mathrm{ns}^{-1}$)")
axes[2].set_ylabel("Correction Error (%)")
cbar2 = fig.colorbar(sc2, ax=axes[2])
cbar2.set_label(r"$|\Delta x|$ (ns)")
axes[2].legend()
axes[2].grid(True, linestyle=':', alpha=0.6)

plt.tight_layout()

# Save output as EPS vector graphic file
output_filename = "csi_integral_shift_error.eps"
plt.savefig(output_filename, format='eps', bbox_inches='tight')