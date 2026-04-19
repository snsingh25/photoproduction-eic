#!/usr/bin/env python3
"""
Integrated Jet Shape Plotting Script
Reads histogram data from C++ analysis output and creates publication-quality plots.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.patches import Patch
import sys
import os

# ============================================
# STANDARD FORMATTING
# ============================================

# LaTeX rendering
rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman', 'DejaVu Serif', 'Times']
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb} \usepackage{bm}'

# Font sizes
rcParams['font.size'] = 32
rcParams['axes.labelsize'] = 30
rcParams['axes.titlesize'] = 28
rcParams['xtick.labelsize'] = 30
rcParams['ytick.labelsize'] = 30
rcParams['legend.fontsize'] = 36

# Figure quality
rcParams['figure.dpi'] = 100
rcParams['savefig.dpi'] = 300

# Axis and tick styling
rcParams['axes.linewidth'] = 2
rcParams['xtick.major.size'] = 10
rcParams['ytick.major.size'] = 10
rcParams['xtick.minor.size'] = 6
rcParams['ytick.minor.size'] = 6
rcParams['xtick.major.width'] = 2
rcParams['ytick.major.width'] = 2
rcParams['xtick.minor.width'] = 1.5
rcParams['ytick.minor.width'] = 1.5
rcParams['xtick.direction'] = 'in'
rcParams['ytick.direction'] = 'in'

# Hatching
rcParams['hatch.linewidth'] = 1.5

# ============================================
# STANDARD COLORS
# ============================================
color_QQ = '#CC0000'  # Red
color_GG = '#0066CC'  # Blue
color_GQ = '#00AA00'  # Green


def read_int_jet_shape_data(filename):
    """Read histogram data from C++ output file."""
    metadata = {}
    bin_centers = []
    qq_vals = []
    gg_vals = []
    gq_vals = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('# Mean_QQ:'):
                metadata['mean_qq'] = float(line.split(':')[1])
            elif line.startswith('# Mean_GG:'):
                metadata['mean_gg'] = float(line.split(':')[1])
            elif line.startswith('# Mean_GQ:'):
                metadata['mean_gq'] = float(line.split(':')[1])
            elif line.startswith('# Entries_QQ:'):
                metadata['entries_qq'] = int(line.split(':')[1])
            elif line.startswith('# Entries_GG:'):
                metadata['entries_gg'] = int(line.split(':')[1])
            elif line.startswith('# Entries_GQ:'):
                metadata['entries_gq'] = int(line.split(':')[1])
            elif line.startswith('# etaMin'):
                metadata['etaMin'] = float(line.split('=')[1].strip())
            elif line.startswith('# etaMax'):
                metadata['etaMax'] = float(line.split('=')[1].strip())
            elif line.startswith('# leading_jet_etMin'):
                metadata['leading_et'] = float(line.split('=')[1].strip())
            elif line.startswith('# subleading_jet_etMin'):
                metadata['subleading_et'] = float(line.split('=')[1].strip())
            elif line.startswith('# Integrated radius'):
                metadata['int_radius'] = float(line.split('=')[1].strip())
            elif line.startswith('# R ='):
                metadata['R'] = float(line.split('=')[1].strip())
            elif line.startswith('#'):
                continue
            elif line:
                parts = line.split()
                bin_centers.append(float(parts[0]))
                qq_vals.append(float(parts[1]))
                gg_vals.append(float(parts[2]))
                gq_vals.append(float(parts[3]))
    
    return np.array(bin_centers), np.array(qq_vals), np.array(gg_vals), np.array(gq_vals), metadata


def plot_int_jet_shape(input_file, output_file=None, include_gq=False):
    """Create publication-quality integrated jet shape plot."""
    
    # Read data
    bin_centers, qq, gg, gq, meta = read_int_jet_shape_data(input_file)
    
    # Determine output filename
    if output_file is None:
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_file = f'{base_name}.pdf'
    
    # Create figure
    fig, ax = plt.subplots(figsize=(13, 13))
    
    # Calculate bin edges for step histogram
    bin_width = bin_centers[1] - bin_centers[0]
    bin_edges = np.append(bin_centers - bin_width/2, bin_centers[-1] + bin_width/2)
    
    # Plot as step histograms - outline
    ax.stairs(qq, bin_edges, color=color_QQ, linewidth=2.5, label=r'\textbf{QQ}')
    ax.stairs(gg, bin_edges, color=color_GG, linewidth=2.5, label=r'\textbf{GG}')
    if include_gq:
        ax.stairs(gq, bin_edges, color=color_GQ, linewidth=2.5, label=r'\textbf{GQ}')
    
    # Add hatching
    ax.stairs(qq, bin_edges, edgecolor=color_QQ, linewidth=0, fill=True, facecolor='none', hatch='//')
    ax.stairs(gg, bin_edges, edgecolor=color_GG, linewidth=0, fill=True, facecolor='none', hatch='\\\\')
    if include_gq:
        ax.stairs(gq, bin_edges, edgecolor=color_GQ, linewidth=0, fill=True, facecolor='none', hatch='..')
    
    # Labels
    int_radius = meta.get('int_radius', 0.3)
    ax.set_xlabel(r'$\mathbf{\Psi(r=' + f'{int_radius}' + r')}$')
    ax.set_ylabel(r'$\mathbf{N_{\mathrm{jets}}}$')
    
    # Axis limits
    ax.set_xlim(0, 1)
    max_val = max(np.max(qq), np.max(gg))
    if include_gq:
        max_val = max(max_val, np.max(gq))
    ax.set_ylim(1, max_val * 1.35)

    # Grey shaded region for ambiguous zone (0.6 < psi < 0.8)
    ax.axvspan(0.6, 0.8, alpha=0.3, color='grey', edgecolor='grey', linewidth=2)

    # Arrows and labels for jet classification regions
    # Thin jets (right side, psi > 0.8)
    # Thin jets (right side, psi > 0.8) - arrow pointing right, text above
    ax.annotate('', xy=(0.99, max_val * 1.2), xytext=(0.80, max_val * 1.2),arrowprops=dict(arrowstyle='->', color='black', lw=2))
    ax.text(0.90, max_val * 1.135, r'\textbf{Thin jets}', fontsize=28, ha='center', va='bottom')

    # Thick jets (left side, psi < 0.6) - arrow pointing left, text above
    ax.annotate('', xy=(0.45, max_val * 1.2), xytext=(0.60, max_val * 1.2), arrowprops=dict(arrowstyle='->', color='black', lw=2))
    ax.text(0.50, max_val * 1.135, r'\textbf{Thick jets}', fontsize=28, ha='center', va='bottom')
            
    # Legend with hatch boxes
    legend_qq = Patch(facecolor='none', edgecolor=color_QQ, hatch='//', label=r'\textbf{QQ}', linewidth=1.5)
    legend_gg = Patch(facecolor='none', edgecolor=color_GG, hatch='\\\\', label=r'\textbf{GG}', linewidth=1.5)
    # if include_gq:
    #     legend_gq = Patch(facecolor='none', edgecolor=color_GQ, hatch='..', label=r'\textbf{GQ}', linewidth=1.5)
    #     ax.legend(handles=[legend_qq, legend_gg, legend_gq], loc='upper right', frameon=False, fontsize=24, handlelength=2, handleheight=2)
    # else:
    #     ax.legend(handles=[legend_qq, legend_gg], loc='upper right', frameon=False, fontsize=24, handlelength=2, handleheight=2)
    if include_gq:
        legend_gq = Patch(facecolor='none', edgecolor=color_GQ, hatch='..', label=r'\textbf{GQ}', linewidth=1.5)
        ax.legend(handles=[legend_qq, legend_gg, legend_gq], loc='upper left', bbox_to_anchor=(0.03, 0.85), frameon=False, fontsize=30, handlelength=2, handleheight=2.5)
    else:
        ax.legend(handles=[legend_qq, legend_gg], loc='upper left', bbox_to_anchor=(0.03, 0.95), frameon=False, fontsize=42, handlelength=2, handleheight=2.5)
    
    # Text annotations
    etaMin = -1.0 #int(meta.get('etaMin', -4))
    etaMax = 2 # int(meta.get('etaMax', 4))
    
    ax.text(0.05, 0.95, rf'$\mathbf{{{etaMin} < \eta_{{\mathrm{{jet}}}} < {etaMax}}}$', transform=ax.transAxes, fontsize=40, verticalalignment='top')
    ax.text(0.05, 0.89, r'ep, 300 GeV', transform=ax.transAxes, fontsize=40, verticalalignment='top')
    
    # Ticks on all sides
    ax.minorticks_on()
    ax.tick_params(which='both', direction='in', top=True, right=True)
    
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight', pad_inches=0.15)
    print(f'Plot saved: {output_file}')
    # plt.show()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python plot_int_jet_shape.py <input_file.txt> [output_file.pdf] [--include-gq]")
        print("Example: python plot_int_jet_shape.py int_jet_shape_-4_4.txt")
        print("         python plot_int_jet_shape.py int_jet_shape_-4_4.txt --include-gq")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = None
    include_gq = True
    
    for arg in sys.argv[2:]:
        if arg == '--include-gq':
            include_gq = True
        elif not arg.startswith('--'):
            output_file = arg
    
    plot_int_jet_shape(input_file, output_file, include_gq)