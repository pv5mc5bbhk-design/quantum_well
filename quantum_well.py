"""
Finite Quantum Well Solver
--------------------------
Computes bound-state energies and normalized wavefunctions
for a symmetric finite square well.

Input parameters:
  m  – particle mass (kg)
  V0 – barrier height (eV)
  L  – well width (nm)
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq
from scipy.integrate import quad

# ── Physical constants ────────────────────────────────────────────────────────
hbar = 1.054571817e-34   # J·s
eV   = 1.602176634e-19   # J per eV

# ── User input ────────────────────────────────────────────────────────────────
print("請輸入參數 / Enter parameters:")
m_kg  = float(input("  粒子質量 m  (kg)  [e.g. 9.11e-31]: "))
V0_eV = float(input("  勢壘高度 V0 (eV)  [e.g. 10]      : "))
L_nm  = float(input("  量子井寬度 L (nm) [e.g. 1]        : "))

V0 = V0_eV * eV          # barrier in Joules
L  = L_nm  * 1e-9        # well width in metres
a  = L / 2               # half-width  (well spans –a … +a)

# ── Transcendental equations for bound states ─────────────────────────────────
#
#  Inside  |x| < a : ψ ~ cos(kx)  [even]  or  sin(kx)  [odd]
#  Outside |x| > a : ψ ~ exp(–κ|x|)
#
#  Matching conditions give:
#    Even:  κ = k·tan(ka)
#    Odd :  κ = –k·cot(ka)   (equivalently k·cot(ka) + κ = 0 → κ = –k·cot)
#
#  with  k² + κ² = 2mV0/ħ²  (fixed circle)

def k_kappa(E_J):
    """Return (k, κ) for energy E (Joules, 0 < E < V0)."""
    k = np.sqrt(2 * m_kg * E_J) / hbar
    kappa = np.sqrt(2 * m_kg * (V0 - E_J)) / hbar
    return k, kappa

def f_even(E_J):
    k, kappa = k_kappa(E_J)
    return k * np.tan(k * a) - kappa

def f_odd(E_J):
    k, kappa = k_kappa(E_J)
    return k / np.tan(k * a) + kappa   # k·cot(ka) + κ = 0  ⟹  return value → 0

# ── Root finding ──────────────────────────────────────────────────────────────
E_min = 1e-6 * eV
E_max = V0 * (1 - 1e-9)

N_scan = 50_000
E_scan = np.linspace(E_min, E_max, N_scan)

even_roots, odd_roots = [], []

for fn, roots in [(f_even, even_roots), (f_odd, odd_roots)]:
    vals = fn(E_scan)
    for i in range(len(vals) - 1):
        # sign change AND not a pole (avoid tan/cot discontinuities)
        if vals[i] * vals[i+1] < 0 and abs(vals[i] - vals[i+1]) < 0.5 * abs(vals).max():
            try:
                root = brentq(fn, E_scan[i], E_scan[i+1], xtol=1e-30, rtol=1e-12)
                roots.append(root)
            except ValueError:
                pass

# Merge and sort by energy; label parity
states = [(E, 'even') for E in even_roots] + [(E, 'odd') for E in odd_roots]
states.sort(key=lambda x: x[0])

print(f"\n束縛態 / Bound states (V0 = {V0_eV} eV, L = {L_nm} nm):")
for idx, (E_J, parity) in enumerate(states, 1):
    print(f"  {idx}: {parity:4s},  E = {E_J/eV:.6f} eV")

# ── Wavefunction construction & normalization ─────────────────────────────────

def wavefunction(x_arr, E_J, parity):
    """
    Returns the (unnormalized) wavefunction on x_arr (metres).
    Amplitude inside the well is set to 1 before normalization.
    """
    k, kappa = k_kappa(E_J)
    psi = np.zeros_like(x_arr, dtype=float)

    inside  = np.abs(x_arr) <= a
    outside = ~inside

    if parity == 'even':
        psi[inside]  = np.cos(k * x_arr[inside])
        # Matching at x = +a: cos(ka) · exp(–κ(|x|–a))
        A_out = np.cos(k * a)
        psi[outside] = A_out * np.exp(-kappa * (np.abs(x_arr[outside]) - a))
    else:
        psi[inside]  = np.sin(k * x_arr[inside])
        # Matching at x = +a (right side positive, left side negative)
        A_out = np.sin(k * a)
        sign  = np.sign(x_arr[outside])
        psi[outside] = sign * A_out * np.exp(-kappa * (np.abs(x_arr[outside]) - a))

    # Normalize: ∫|ψ|² dx = 1
    norm_sq, _ = quad(lambda x: wavefunction_scalar(x, E_J, parity) ** 2,
                      -10 * a, 10 * a, limit=200)
    psi /= np.sqrt(norm_sq)
    return psi


def wavefunction_scalar(x, E_J, parity):
    """Scalar version used by quad for normalization."""
    k, kappa = k_kappa(E_J)
    if abs(x) <= a:
        return np.cos(k * x) if parity == 'even' else np.sin(k * x)
    else:
        A_out = np.cos(k * a) if parity == 'even' else np.sin(k * a)
        val   = A_out * np.exp(-kappa * (abs(x) - a))
        return val if parity == 'even' else np.sign(x) * val


# ── Plotting ──────────────────────────────────────────────────────────────────
x_nm   = np.linspace(-3 * L_nm, 3 * L_nm, 4000)
x_m    = x_nm * 1e-9

colors = plt.cm.tab10(np.linspace(0, 0.9, len(states)))
scale  = 0.8 * V0_eV          # display scale for wavefunctions (eV)

fig, ax = plt.subplots(figsize=(10, 6))

# Potential well outline
well_x = [-3*L_nm, -L_nm/2, -L_nm/2, L_nm/2, L_nm/2, 3*L_nm]
well_V = [V0_eV,   V0_eV,   0,        0,       V0_eV,  V0_eV]
ax.plot(well_x, well_V, 'k-', lw=2, label='V(x)')
ax.axhline(0, color='k', lw=0.5, ls='--')

for idx, ((E_J, parity), color) in enumerate(zip(states, colors), 1):
    E_eV = E_J / eV
    psi  = wavefunction(x_m, E_J, parity)

    # Offset wavefunction by its energy level for a clean display
    ax.axhline(E_eV, color=color, lw=0.6, ls=':', alpha=0.6)
    ax.plot(x_nm, E_eV + scale * psi, color=color, lw=1.6,
            label=f'{idx}: {parity} ({E_eV:.3f} eV)')

ax.set_xlabel('x (nm)', fontsize=12)
ax.set_ylabel('Energy (eV)', fontsize=12)
ax.set_title(f'Finite Quantum Well  |  V₀={V0_eV} eV, L={L_nm} nm, '
             f'm={m_kg:.2e} kg', fontsize=12)
ax.legend(fontsize=9, loc='upper right')
ax.set_xlim(x_nm[0], x_nm[-1])
ax.set_ylim(-0.5, V0_eV * 1.3)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig('quantum_well_plot.png', dpi=150)
plt.show()
print("\nPlot saved as quantum_well_plot.png")
