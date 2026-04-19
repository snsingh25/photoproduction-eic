#!/usr/bin/env python3
"""
Integrated Jet Shape Comparison: PYTHIA 8 vs ZEUS Data
=======================================================
Tight 2x2 grid with shared axes.

Author: Siddharth Singh
Usage: python jetrecointrange_plot.py <logfile.log>
"""

import sys
import re
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import make_interp_spline
from pathlib import Path

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

# =============================================================================
# ZEUS DATA (Table 1 of ZEUS paper)
# =============================================================================
ZEUS_DATA = {
    'eta1': {  # -1 < eta < 0
        'r': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        'psi': [0.3135, 0.5870, 0.7322, 0.8165, 0.8720, 0.9106, 0.9382, 0.9577, 0.9719, 0.9817],
        'stat': [0.0025, 0.0025, 0.0020, 0.0015, 0.0011, 0.0008, 0.0006, 0.0004, 0.0003, 0.0002],
    },
    'eta2': {  # 0 < eta < 1
        'r': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        'psi': [0.3040, 0.5730, 0.7189, 0.8075, 0.8649, 0.9048, 0.9335, 0.9541, 0.9690, 0.9797],
        'stat': [0.0012, 0.0012, 0.0010, 0.0007, 0.0006, 0.0004, 0.0003, 0.0002, 0.0002, 0.0001],
    },
    'eta3': {  # 1 < eta < 1.5
        'r': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        'psi': [0.2678, 0.5238, 0.6740, 0.7721, 0.8378, 0.8848, 0.9191, 0.9436, 0.9619, 0.9749],
        'stat': [0.0016, 0.0018, 0.0016, 0.0012, 0.0009, 0.0007, 0.0005, 0.0004, 0.0003, 0.0002],
    },
    'eta4': {  # 1.5 < eta < 2.5
        'r': [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
        'psi': [0.2236, 0.4581, 0.6121, 0.7180, 0.7936, 0.8488, 0.8903, 0.9215, 0.9457, 0.9637],
        'stat': [0.0012, 0.0013, 0.0012, 0.0010, 0.0008, 0.0006, 0.0005, 0.0004, 0.0003, 0.0002],
    },
}

ETA_LABELS = {
    'eta1': r'$-1 < \eta^{\mathrm{jet}} < 0$',
    'eta2': r'$0 < \eta^{\mathrm{jet}} < 1$',
    'eta3': r'$1 < \eta^{\mathrm{jet}} < 1.5$',
    'eta4': r'$1.5 < \eta^{\mathrm{jet}} < 2.5$',
}


# =============================================================================
# PARSING
# =============================================================================

def parse_log_file(filepath):
    with open(filepath, 'r') as f:
        content = f.read()

    data = {
        'QQ': {'eta1': [], 'eta2': [], 'eta3': [], 'eta4': []},
        'GG': {'eta1': [], 'eta2': [], 'eta3': [], 'eta4': []},
        'Combined': {'eta1': [], 'eta2': [], 'eta3': [], 'eta4': []},
    }

    row_pattern = re.compile(
        r'^\s*([\d.]+)\s+'
        r'([\d.]+)\s*\+/-\s*([\d.]+)\s+'
        r'([\d.]+)\s*\+/-\s*([\d.]+)\s+'
        r'([\d.]+)\s*\+/-\s*([\d.]+)\s+'
        r'([\d.]+)\s*\+/-\s*([\d.]+)',
        re.MULTILINE
    )

    sections = {
        'QQ': re.search(r'>>> QQ_Events <<<.*?(?=>>> [A-Z]|$)', content, re.DOTALL),
        'GG': re.search(r'>>> GG_Events <<<.*?(?=>>> [A-Z]|$)', content, re.DOTALL),
        'Combined': re.search(r'>>> Combined_Events <<<.*?(?=>>> [A-Z]|$|={5,})', content, re.DOTALL),
    }

    for event_type, match in sections.items():
        if match:
            section_text = match.group(0)
            rows = row_pattern.findall(section_text)
            for row in rows:
                r_val = float(row[0])
                data[event_type]['eta1'].append({'r': r_val, 'psi': float(row[1]), 'err': float(row[2])})
                data[event_type]['eta2'].append({'r': r_val, 'psi': float(row[3]), 'err': float(row[4])})
                data[event_type]['eta3'].append({'r': r_val, 'psi': float(row[5]), 'err': float(row[6])})
                data[event_type]['eta4'].append({'r': r_val, 'psi': float(row[7]), 'err': float(row[8])})

    return data


def extract_arrays(data_list):
    r   = np.array([d['r']   for d in data_list])
    psi = np.array([d['psi'] for d in data_list])
    err = np.array([d['err'] for d in data_list])
    return r, psi, err


def smooth(x, y, num=200):
    spl = make_interp_spline(x, y, k=3)
    x_new = np.linspace(x.min(), x.max(), num)
    return x_new, spl(x_new)


# =============================================================================
# PLOTTING
# =============================================================================

def create_plot(sim_data, output_path):
    """Tight 2x2 with shared axes."""

    fig, axes = plt.subplots(2, 2, figsize=(10, 9),
                              sharex=True, sharey=True,
                              gridspec_kw={'hspace': 0.0, 'wspace': 0.0})

    eta_keys = ['eta1', 'eta2', 'eta3', 'eta4']
    quark_color = '#D62728'
    gluon_color = '#1F77B4'

    for idx, eta_key in enumerate(eta_keys):
        row, col = divmod(idx, 2)
        ax = axes[row, col]

        # ZEUS data — stat error bars only
        zeus_r   = np.array(ZEUS_DATA[eta_key]['r'])
        zeus_psi = np.array(ZEUS_DATA[eta_key]['psi'])
        zeus_stat = np.array(ZEUS_DATA[eta_key]['stat'])

        ax.errorbar(zeus_r, zeus_psi, yerr=zeus_stat,
                     fmt='o', color='black', markersize=6,
                     capsize=2.5, capthick=1, elinewidth=1,
                     label='ZEUS Data', zorder=10)

        # Simulation (smooth)
        qq_r, qq_psi, _ = extract_arrays(sim_data['QQ'][eta_key])
        gg_r, gg_psi, _ = extract_arrays(sim_data['GG'][eta_key])

        r_s, qq_s = smooth(qq_r, qq_psi)
        r_s, gg_s = smooth(gg_r, gg_psi)

        ax.plot(r_s, qq_s, color=quark_color, ls='--', lw=2.5,
                label=r'Quark jets ($q\bar{q}$)')
        ax.plot(r_s, gg_s, color=gluon_color, ls='-.', lw=2.5,
                label=r'Gluon jets ($gg$)')

        # Limits
        ax.set_xlim(0.05, 1.05)
        ax.set_ylim(0.0, 1.05)
        ax.minorticks_on()
        ax.tick_params(axis='both', which='major', labelsize=16)
        ax.tick_params(which='both', direction='in', top=True, right=True)

        # η label
        ax.text(0.06, 0.95, ETA_LABELS[eta_key],
                transform=ax.transAxes, fontsize=17, verticalalignment='top')

        # Legend — top-right panel only
        if eta_key == 'eta2':
            ax.legend(loc='lower right', fontsize=14, frameon=False)

    # Shared axis labels — closer to axes
    fig.supxlabel(r'$r$', fontsize=24, y=0.035)
    fig.supylabel(r'$\langle\psi(r)\rangle$', fontsize=24, x=0.015)

    fig.savefig(output_path, format='pdf', bbox_inches='tight', dpi=300)
    print(f"Saved → {output_path}")
    plt.close(fig)


# =============================================================================
# MAIN
# =============================================================================

def main():
    if len(sys.argv) < 2:
        print("Usage: python jetrecointrange_plot.py <logfile.log>")
        sys.exit(1)

    log_file = Path(sys.argv[1])
    if not log_file.exists():
        print(f"Error: File not found: {log_file}")
        sys.exit(1)

    print(f"Parsing: {log_file}")
    sim_data = parse_log_file(log_file)

    for evt in ['QQ', 'GG', 'Combined']:
        n = len(sim_data[evt]['eta1'])
        print(f"  {evt}: {n} r-points per eta bin")

    output_path = Path('./') / (log_file.stem + '.pdf')
    print("Creating plot...")
    create_plot(sim_data, output_path)
    print("Done!")


if __name__ == '__main__':
    main()