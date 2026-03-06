import matplotlib.pyplot as plt
import numpy as np
import os

# --- Measured data ---
n = [50, 200, 1000]

naive_ms   = [25.0,  404.0, 10100.0]   # n=1000 extrapolated via O(n^2): 404*(1000/200)^2
spatial_ms = [ 2.3,    7.3,    52.0]

naive_pairs   = [1_225,  19_900, 499_500]
spatial_pairs = [  402,   1_387,  13_314]   # peak (first interval, objects still falling)

out_dir = os.path.dirname(os.path.abspath(__file__))

# ── Figure 4: timing comparison ──────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(6, 4))
ax.plot(n, naive_ms,   'r-o', label=r'Naïve $O(n^2)$')
ax.plot(n, spatial_ms, 'b-o', label='Spatial hash')
ax.set_xlabel('$n$ (objects)')
ax.set_ylabel('Broad-phase time (ms / step)')
ax.set_title('Broad-phase performance vs. scene size')
ax.legend()
ax.set_xticks(n)
ax.set_xticklabels([str(x) for x in n])
ax.grid(True, linestyle='--', alpha=0.4)
# Annotate n=1000 naive as extrapolated
ax.annotate('†extrapolated', xy=(1000, naive_ms[2]),
            xytext=(700, naive_ms[2] * 0.75),
            arrowprops=dict(arrowstyle='->', color='red'),
            color='red', fontsize=8)
fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'fig4_perf.pdf'))
fig.savefig(os.path.join(out_dir, 'fig4_perf.png'), dpi=150)
print('Saved fig4_perf.pdf / .png')

# ── Figure 5: candidate pairs comparison ─────────────────────────────────────
fig, axes = plt.subplots(1, 2, figsize=(10, 4))

# Left: bar chart of pairs at peak load
x = np.arange(len(n))
width = 0.35
axes[0].bar(x - width/2, naive_pairs,   width, color='red',  alpha=0.7, label=r'Naïve $O(n^2)$')
axes[0].bar(x + width/2, spatial_pairs, width, color='blue', alpha=0.7, label='Spatial hash')
axes[0].set_xticks(x)
axes[0].set_xticklabels([f'$n={v}$' for v in n])
axes[0].set_ylabel('Candidate pairs checked')
axes[0].set_title('Candidate pairs at peak load')
axes[0].legend()
axes[0].grid(True, axis='y', linestyle='--', alpha=0.4)
# Annotate reduction at n=1000
axes[0].annotate(
    f'{naive_pairs[2]//spatial_pairs[2]}× fewer\npairs',
    xy=(x[2] + width/2, spatial_pairs[2]),
    xytext=(x[2] - 0.6, spatial_pairs[2] + 40000),
    arrowprops=dict(arrowstyle='->', color='black'),
    fontsize=8)

# Right: pairs vs n (log scale to show asymptotic difference)
axes[1].plot(n, naive_pairs,   'r-o', label=r'Naïve $O(n^2)$')
axes[1].plot(n, spatial_pairs, 'b-o', label='Spatial hash')
axes[1].set_xlabel('$n$ (objects)')
axes[1].set_ylabel('Candidate pairs')
axes[1].set_title('Pairs vs. $n$ (log scale)')
axes[1].set_yscale('log')
axes[1].set_xticks(n)
axes[1].set_xticklabels([str(x) for x in n])
axes[1].legend()
axes[1].grid(True, which='both', linestyle='--', alpha=0.4)

fig.tight_layout()
fig.savefig(os.path.join(out_dir, 'fig5_pairs.pdf'))
fig.savefig(os.path.join(out_dir, 'fig5_pairs.png'), dpi=150)
print('Saved fig5_pairs.pdf / .png')
