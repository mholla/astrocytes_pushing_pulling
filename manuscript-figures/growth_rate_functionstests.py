import numpy as np
import matplotlib.pyplot as plt

# Parameters
delta, alpha = 1.0, 0.3
N = 4
threshold = 0.8
steepness = 20
dpi = 600
# Create radial profile for comparison
r_values = np.linspace(0, 2, 500)
psi_fixed = np.pi/4  # Fixed angle for 1D plot

# Original sharp transition
def phi_original(r, psi):
    gaussian = (1 / (2 * np.sqrt(2 * np.pi * alpha**2))) * \
               np.exp(-(r - delta)**2 / (2 * alpha**2))
    sine_part = np.sin(4 * psi * (N - 0.5)) + 1
    return np.where(r < threshold, 1.0, gaussian * sine_part)

# Various smooth transitions
def phi_sigmoid(r, psi):
    gaussian = (1 / (2 * np.sqrt(2 * np.pi * alpha**2))) * \
               np.exp(-(r - delta)**2 / (2 * alpha**2))
    sine_part = np.sin(4 * psi * (N - 0.5)) + 1
    growth = gaussian * sine_part
    transition = 1 / (1 + np.exp(-steepness * (r - threshold)))
    return 1 * (1 - transition) + growth * transition

def phi_smoothstep(r, psi, width=0.15):
    gaussian = (1 / (2 * np.sqrt(2 * np.pi * alpha**2))) * \
               np.exp(-(r - delta)**2 / (2 * alpha**2))
    sine_part = np.sin(4 * psi * (N - 0.5)) + 1
    growth = gaussian * sine_part
    t = (r - (threshold - width/2)) / width
    t = np.clip(t, 0, 1)
    t_smooth = t**3 * (t * (t * 6 - 15) + 10)
    return 1 * (1 - t_smooth) + growth * t_smooth

# Plot comparison
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# 1D profiles
ax = axes[0, 0]
ax.plot(r_values, phi_original(r_values, psi_fixed), 'k--', label='Sharp threshold', linewidth=2)
ax.plot(r_values, phi_sigmoid(r_values, psi_fixed), 'b-', label='Sigmoid (steepness=20)', linewidth=2)
ax.axvline(threshold, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel(r'$\tilde{r}$')
ax.set_ylabel(r'$\phi$')
ax.set_title('1D Profile Comparison')
ax.legend()
ax.grid(True, alpha=0.3)

ax = axes[0, 1]
ax.plot(r_values, phi_original(r_values, psi_fixed), 'k--', label='Sharp threshold', linewidth=2)
ax.plot(r_values, phi_smoothstep(r_values, psi_fixed, width=0.1), 'r-', label='Smoothstep (width=0.1)', linewidth=2)
ax.plot(r_values, phi_smoothstep(r_values, psi_fixed, width=0.2), 'g-', label='Smoothstep (width=0.2)', linewidth=2)
ax.axvline(threshold, color='gray', linestyle=':', alpha=0.5)
ax.set_xlabel(r'$\tilde{r}$')
ax.set_ylabel(r'$\phi$')
ax.set_title('Smoothstep with Different Widths')
ax.legend()
ax.grid(True, alpha=0.3)

# 2D visualization of smoothstep
ax = axes[1, 0]
X, Y = np.meshgrid(np.linspace(-6, 6, 100), np.linspace(-6, 6, 100))
r_grid = np.sqrt(X**2 + Y**2)  # Simplified for demo
psi_grid = np.arctan2(Y, X)
phi_grid = phi_smoothstep(r_grid, psi_grid, width=0.15)
im = ax.imshow(phi_grid, extent=[-6, 6, -6, 6], origin='lower', cmap='hot')
ax.set_title('2D Smooth Transition (Smoothstep)')
ax.set_xlabel('X')
ax.set_ylabel('Y')
plt.colorbar(im, ax=ax)

# 2D visualization of sigmoid
ax = axes[1, 1]
phi_grid_sigmoid = phi_sigmoid(r_grid, psi_grid)
im = ax.imshow(phi_grid_sigmoid, extent=[-6, 6, -6, 6], origin='lower', cmap='hot')
ax.set_title('2D Smooth Transition (Sigmoid)')
ax.set_xlabel('X')
ax.set_ylabel('Y')
plt.colorbar(im, ax=ax)

plt.tight_layout()
plt.savefig('smooth_transitions.pdf', dpi=600, bbox_inches='tight')