#!/usr/bin/env python3
"""
Jet Shape Cut Comparison: E_T > 17 GeV vs Dijet (10, 7 GeV)
============================================================
Tight 2x2 shared-axes, smooth curves, publication-quality.

Author: Siddharth Singh
Output: jetshape_cut_comparison.pdf
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline

# ── Style ────────────────────────────────────────────────────────────────────
plt.style.use('default')
plt.rcParams.update({
    'font.size': 18,
    'font.family': 'serif',
    'text.usetex': False,          # flip to True on your machine
    'mathtext.fontset': 'cm',
    'axes.linewidth': 1.4,
    'axes.spines.top': True,
    'axes.spines.right': True,
    'xtick.major.size': 7,
    'xtick.minor.size': 4,
    'ytick.major.size': 7,
    'ytick.minor.size': 4,
    'xtick.major.width': 1.0,
    'ytick.major.width': 1.0,
    'xtick.minor.width': 0.7,
    'ytick.minor.width': 0.7,
    'lines.linewidth': 2.5,
    'legend.frameon': False,
    'figure.facecolor': 'white',
    'axes.facecolor': 'white',
})

# ── Data ─────────────────────────────────────────────────────────────────────
R_VALUES = np.array([0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])

ET17_QQ = {
    'eta1': np.array([0.3066, 0.5830, 0.7249, 0.8076, 0.8597, 0.8975, 0.9261, 0.9487, 0.9665, 0.9762]),
    'eta2': np.array([0.3125, 0.5816, 0.7249, 0.8074, 0.8600, 0.8972, 0.9270, 0.9506, 0.9692, 0.9791]),
    'eta3': np.array([0.3303, 0.6029, 0.7428, 0.8245, 0.8780, 0.9145, 0.9425, 0.9644, 0.9822, 0.9933]),
    'eta4': np.array([0.3315, 0.6074, 0.7465, 0.8269, 0.8801, 0.9168, 0.9441, 0.9650, 0.9826, 0.9966]),
}
ET17_GG = {
    'eta1': np.array([0.1082, 0.2799, 0.4352, 0.5780, 0.6791, 0.7626, 0.8319, 0.8856, 0.9289, 0.9491]),
    'eta2': np.array([0.0994, 0.2775, 0.4404, 0.5718, 0.6769, 0.7625, 0.8313, 0.8897, 0.9379, 0.9593]),
    'eta3': np.array([0.1191, 0.3172, 0.4870, 0.6189, 0.7221, 0.8054, 0.8717, 0.9275, 0.9687, 0.9876]),
    'eta4': np.array([0.1195, 0.3189, 0.4934, 0.6301, 0.7352, 0.8198, 0.8873, 0.9391, 0.9760, 0.9913]),
}

DIJET_QQ = {
    'eta1': np.array([0.2012, 0.4647, 0.6392, 0.7485, 0.8202, 0.8718, 0.9106, 0.9407, 0.9642, 0.9735]),
    'eta2': np.array([0.2065, 0.4658, 0.6376, 0.7450, 0.8162, 0.8676, 0.9069, 0.9384, 0.9627, 0.9727]),
    'eta3': np.array([0.2161, 0.4805, 0.6531, 0.7614, 0.8327, 0.8839, 0.9224, 0.9533, 0.9786, 0.9909]),
    'eta4': np.array([0.2124, 0.4772, 0.6515, 0.7614, 0.8350, 0.8867, 0.9254, 0.9556, 0.9803, 0.9952]),
}
DIJET_GG = {
    'eta1': np.array([0.0817, 0.2469, 0.4121, 0.5528, 0.6642, 0.7537, 0.8273, 0.8890, 0.9334, 0.9511]),
    'eta2': np.array([0.0768, 0.2358, 0.3978, 0.5368, 0.6509, 0.7443, 0.8222, 0.8871, 0.9364, 0.9563]),
    'eta3': np.array([0.0796, 0.2453, 0.4138, 0.5591, 0.6785, 0.7754, 0.8544, 0.9197, 0.9668, 0.9858]),
    'eta4': np.array([0.0805, 0.2459, 0.4182, 0.5679, 0.6901, 0.7896, 0.8700, 0.9315, 0.9731, 0.9865]),
}

ETA_LABELS = {
    'eta1': r'$-1 < \eta^{\mathrm{jet}} < 0$',
    'eta2': r'$0 < \eta^{\mathrm{jet}} < 1$',
    'eta3': r'$1 < \eta^{\mathrm{jet}} < 1.5$',
    'eta4': r'$1.5 < \eta^{\mathrm{jet}} < 2.5$',
}


def smooth(x, y, num=200):
    spl = make_interp_spline(x, y, k=3)
    x_new = np.linspace(x.min(), x.max(), num)
    return x_new, spl(x_new)


def create_plot():
    fig, axes = plt.subplots(2, 2, figsize=(10, 9),
                              sharex=True, sharey=True,
                              gridspec_kw={'hspace': 0.0, 'wspace': 0.0})

    eta_keys = ['eta1', 'eta2', 'eta3', 'eta4']
    quark_color = '#d62728' # Standard MPL Red
    gluon_color = '#1f77b4' # Standard MPL Blue

    for idx, eta_key in enumerate(eta_keys):
        row, col = divmod(idx, 2)
        ax = axes[row, col]

        # Smooth all curves
        r_s, et17_qq_s = smooth(R_VALUES, ET17_QQ[eta_key])
        r_s, et17_gg_s = smooth(R_VALUES, ET17_GG[eta_key])
        r_s, dj_qq_s   = smooth(R_VALUES, DIJET_QQ[eta_key])
        r_s, dj_gg_s   = smooth(R_VALUES, DIJET_GG[eta_key])

        # E_T > 17 GeV — solid
        ax.plot(r_s, et17_qq_s, color=quark_color, ls='-', lw=2.5,
                label=r'Quark jets ($E_T > 17$ GeV)')
        ax.plot(r_s, et17_gg_s, color=gluon_color, ls='-', lw=2.5,
                label=r'Gluon jets ($E_T > 17$ GeV)')

        # Dijet — dashed
        ax.plot(r_s, dj_qq_s, color=quark_color, ls='--', lw=2.5,
                label=r'Quark jets ($E_T > 10$ GeV)')
        ax.plot(r_s, dj_gg_s, color=gluon_color, ls='--', lw=2.5,
                label=r'Gluon jets ($E_T > 10$ GeV)')

        # Limits & ticks
        ax.set_xlim(0.05, 1.05)
        ax.set_ylim(0.0, 1.05)
        ax.minorticks_on()
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.tick_params(which='both', direction='in', top=True, right=True)

        # η label
        ax.text(0.06, 0.95, ETA_LABELS[eta_key],
                transform=ax.transAxes, fontsize=18, verticalalignment='top')

        # Legend — top-right panel only
        if eta_key == 'eta2':
            ax.legend(loc='lower right', fontsize=12, frameon=False)

    # Shared axis labels
    fig.supxlabel(r'$r$', fontsize=30, y=0.035)
    fig.supylabel(r'$\langle\psi(r)\rangle$', fontsize=30, x=0.015)

    out = 'jetshape_cut_comparison.pdf'
    fig.savefig(out, format='pdf', bbox_inches='tight', dpi=300)
    print(f"Saved → {out}")
    plt.close(fig)


if __name__ == '__main__':
    create_plot()