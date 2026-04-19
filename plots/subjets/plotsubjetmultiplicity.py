#!/usr/bin/env python3
"""
Subjet Multiplicity Plotting Script
Reads histogram data from C++ analysis output and creates publication-quality plots.
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib import rcParams
import sys
import re

# Configure matplotlib for LaTeX-style fonts and publication quality
rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb} \usepackage{bm} \renewcommand{\familydefault}{\sfdefault} \boldmath'

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman', 'DejaVu Serif', 'Times']
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb} \usepackage{bm}'
rcParams['font.size'] = 28

rcParams['axes.labelsize'] = 26
rcParams['axes.titlesize'] = 30
rcParams['xtick.labelsize'] = 24
rcParams['ytick.labelsize'] = 24

rcParams['legend.fontsize'] = 22
rcParams['figure.dpi'] = 100
rcParams['savefig.dpi'] = 300
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
rcParams['hatch.linewidth'] = 1.5  # Default is 1.0

# Colors
color_QQ = '#CC0000'  # Red
color_GG = '#0066CC'  # Blue ## '#00AA00'  # Green


def read_subjet_data(filename):
    """Read histogram data from C++ output file."""
    metadata = {}
    bin_centers = []
    qq_vals = []
    gg_vals = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('# Mean_QQ:'):
                metadata['mean_qq'] = float(line.split(':')[1])
            elif line.startswith('# Mean_GG:'):
                metadata['mean_gg'] = float(line.split(':')[1])
            elif line.startswith('# Entries_QQ:'):
                metadata['entries_qq'] = int(line.split(':')[1])
            elif line.startswith('# Entries_GG:'):
                metadata['entries_gg'] = int(line.split(':')[1])
            elif line.startswith('# Experiment:'):
                metadata['experiment'] = line.split(':')[1].strip()
            elif line.startswith('# ycut'):
                metadata['ycut'] = line.split('=')[1].strip()
            elif line.startswith('# etaMin'):
                # Use regex to extract value, preserving decimals
                match = re.search(r'=\s*([-+]?\d*\.?\d+)', line)
                if match:
                    metadata['etaMin'] = float(match.group(1))
            elif line.startswith('# etaMax'):
                # Use regex to extract value, preserving decimals
                match = re.search(r'=\s*([-+]?\d*\.?\d+)', line)
                if match:
                    metadata['etaMax'] = float(match.group(1))
            elif line.startswith('# Leading_Jet_ET_Min'):
                match = re.search(r'=\s*([-+]?\d*\.?\d+)', line)
                if match:
                    metadata['leading_et'] = float(match.group(1))
            elif line.startswith('# Subleading_Jet_ET_Min'):
                match = re.search(r'=\s*([-+]?\d*\.?\d+)', line)
                if match:
                    metadata['subleading_et'] = float(match.group(1))
            elif line.startswith('#'):
                continue
            elif line:
                parts = line.split()
                bin_centers.append(float(parts[0]))
                qq_vals.append(float(parts[1]))
                gg_vals.append(float(parts[2]))
    
    return np.array(bin_centers), np.array(qq_vals), np.array(gg_vals), metadata


def plot_subjet_multiplicity(input_file, output_file=None):
    """Create publication-quality subjet multiplicity plot."""
    
    # Read data
    bin_centers, qq, gg, meta = read_subjet_data(input_file)
    
    # Determine output filename
    if output_file is None:
      import os
      base_name = os.path.splitext(os.path.basename(input_file))[0]
      output_file = f'{base_name}.pdf'
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 8))
    
    # Calculate bin edges for step histogram
    bin_width = bin_centers[1] - bin_centers[0]
    bin_edges = np.append(bin_centers - bin_width/2, bin_centers[-1] + bin_width/2)
    
    # Plot as step histograms with hatching
    ax.stairs(qq, bin_edges, color=color_QQ, linewidth=2.5, label=f'QQ')
    ax.stairs(gg, bin_edges, color=color_GG, linewidth=2.5, label=f'GG')
    # Hatching with color
    ax.stairs(qq, bin_edges, linewidth=0, fill=True, facecolor='none', edgecolor=color_QQ, hatch='///')
    ax.stairs(gg, bin_edges, linewidth=0, fill=True, facecolor='none', edgecolor=color_GG, hatch='\\\\\\')

    legend_qq = Patch(facecolor='none', edgecolor=color_QQ, hatch='///', label='QQ', linewidth=1.5)
    legend_gg = Patch(facecolor='none', edgecolor=color_GG, hatch='\\\\\\', label='GG', linewidth=1.5)
    ax.legend(handles=[legend_qq, legend_gg], loc='upper right', frameon=False, fontsize=24, handlelength=2, handleheight=2)
    
    # Labels and formatting
    ax.set_xlabel(r'$n_{\mathrm{sbj}}$')
    ax.set_ylabel(r'$\frac{1}{N_{\mathrm{jet}}} \frac{dN_{\mathrm{jet}}}{dn_{\mathrm{sbj}}}$')
    
    # Axis limits
    ax.set_xlim(0, 10)
    max_val = max(np.max(qq), np.max(gg))
    ax.set_ylim(0.001, max_val * 1.35)

    # Extract eta values with proper decimal handling
    eta_min = meta.get('etaMin', 0)
    eta_max = meta.get('etaMax', 0)
    
    # Format eta display: use decimal if needed, otherwise show as integer
    if eta_min == int(eta_min):
        eta_min_str = f'{int(eta_min)}'
    else:
        eta_min_str = f'{eta_min:.1f}'
    
    if eta_max == int(eta_max):
        eta_max_str = f'{int(eta_max)}'
    else:
        eta_max_str = f'{eta_max:.1f}'

    # Add eta range and ep energy labels
    ax.text(0.05, 0.95, rf'${eta_min_str} < \eta_{{\mathrm{{jet}}}} < {eta_max_str}$', 
            transform=ax.transAxes, fontsize=28, verticalalignment='top')
    ax.text(0.05, 0.87, r'ep, 300 GeV', transform=ax.transAxes, fontsize=28, verticalalignment='top')
    
    # Minor ticks
    ax.minorticks_on()
    ax.tick_params(which='both', direction='in', top=True, right=True)
    
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight', pad_inches=0.15)
    print(f'Plot saved: {output_file}')
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python plot_subjet_multiplicity.py <input_file.txt> [output_file.pdf]")
        print("Example: python plot_subjet_multiplicity.py subjet_multiplicity_data_hera300_pT7.txt")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2] if len(sys.argv) > 2 else None
    
    plot_subjet_multiplicity(input_file, output_file)