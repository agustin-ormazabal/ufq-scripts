import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
from brokenaxes import brokenaxes
from matplotlib.patches import Rectangle

p2 = np.loadtxt("entropias.p2")
r  = np.loadtxt("entropias.r")
p1 = np.loadtxt("entropias.p1")
x1 = np.loadtxt("entropias.r.sl2_sl2")
x2 = np.loadtxt("entropias.i.sl2_sl2")

fig = plt.figure(figsize=(12,10))
gs = plt.GridSpec(2, 1, height_ratios=[1,1])

bax = brokenaxes(
    xlims=((10,40),),
    ylims=((1800,2100),(6200,6700),(7500,8400)),
    hspace=0.09,
    fig=fig,
    subplot_spec=gs[0],
    d=.00000015
)
bax.tick_params(axis='x', which='both', bottom=False, labelsize=0, labelbottom=False)
bax.tick_params(axis='y', which='both', bottom=False, labelsize=16, labelbottom=False)

bax.plot(r[:,0],  r[:,1],  '-o', label='RsmE-(SL2)$_2$', color='#332d63', lw=2)
bax.plot(p1[:,0], p1[:,1], '-o', label='RsmE-(SL2)$_1$', color='#c87d7d', lw=2)
bax.plot(p2[:,0], p2[:,1], '-o', label='SL2', color='red', lw=2)
bax.legend(fontsize=14, loc='center right', bbox_to_anchor=(0.99, 0.64))
bax.set_ylabel("TΔS (kcal/mol)", fontsize=20, labelpad=70)
bax.grid(axis='y', lw=1, linestyle='--', alpha=0.6)

ax2 = fig.add_subplot(gs[1])
ax2.plot(x1[:,0], x1[:,1], '-o', color='#007BB8', lw=2, ms=6, label='RsmE-(SL-l-SL)')
ax2.plot(x2[:,0], x2[:,1], '-o', color='#5BA15F', lw=2, ms=6, label='RsmE-(SL-l*-SL*)')
ax2.spines['top'].set_visible(False)
ax2.spines['right'].set_visible(False)
ax2.legend(loc='lower right')
ax2.set_xlabel("Number of trajectories", fontsize=20, labelpad=20)
ax2.set_ylabel("TΔS (kcal/mol)", fontsize=20, labelpad=20)
ax2.grid(axis='y', lw=1, linestyle='--', alpha=0.6)
ax2.set_xlim(10, 40)
ax2.set_ylim(7300, 8400)

ax2.tick_params(axis='both', labelsize=16)
ax2.legend(fontsize=14, loc='lower right', bbox_to_anchor=(0.99, 0.085))

plt.tight_layout()
plt.savefig("broken_multiplot.png", dpi=200)

