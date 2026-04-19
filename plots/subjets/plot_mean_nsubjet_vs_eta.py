#!/usr/bin/env python3
"""
Mean Subjet Multiplicity vs Eta Plotting Script
Reads data from multiple C++ analysis outputs and creates publication-quality plots.
Uses spline interpolation for smooth curves.
Overlays multiple energies with different line styles.
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from matplotlib.lines import Line2D
from scipy.interpolate import make_interp_spline
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
rcParams['font.size'] = 30
rcParams['axes.labelsize'] = 40
rcParams['axes.titlesize'] = 34
rcParams['xtick.labelsize'] = 24
rcParams['ytick.labelsize'] = 24
rcParams['legend.fontsize'] = 20

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

# ============================================
# COLORS
# ============================================
color_QQ = '#d62728' # Standard MPL Red
color_GG = '#1f77b4' # Standard MPL Blue


def read_mean_nsubjet_data(filename):
    """Read mean subjet multiplicity data from C++ output file."""
    metadata = {}
    eta_low = []
    eta_high = []
    eta_center = []
    mean_qq = []
    mean_gg = []
    n_jets_qq = []
    n_jets_gg = []
    
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('# Experiment:'):
                metadata['experiment'] = line.split(':')[1].strip()
            elif line.startswith('# sqrt_s'):
                metadata['sqrt_s'] = float(line.split('=')[1].strip())
            elif line.startswith('# R ='):
                metadata['R'] = float(line.split('=')[1].strip())
            elif line.startswith('# ycut'):
                metadata['ycut'] = float(line.split('=')[1].strip())
            elif line.startswith('# leading_jet_etMin'):
                metadata['leading_et'] = float(line.split('=')[1].strip())
            elif line.startswith('# subleading_jet_etMin'):
                metadata['subleading_et'] = float(line.split('=')[1].strip())
            elif line.startswith('# dijet_events'):
                metadata['dijet_events'] = int(line.split('=')[1].strip())
            elif line.startswith('#'):
                continue
            elif line:
                parts = line.split()
                eta_low.append(float(parts[0]))
                eta_high.append(float(parts[1]))
                eta_center.append(float(parts[2]))
                mean_qq.append(float(parts[3]))
                mean_gg.append(float(parts[4]))
                n_jets_qq.append(int(parts[5]))
                n_jets_gg.append(int(parts[6]))
    
    return {
        'eta_low': np.array(eta_low),
        'eta_high': np.array(eta_high),
        'eta_center': np.array(eta_center),
        'mean_qq': np.array(mean_qq),
        'mean_gg': np.array(mean_gg),
        'n_jets_qq': np.array(n_jets_qq),
        'n_jets_gg': np.array(n_jets_gg),
        'metadata': metadata
    }


def plot_multiple_energies(input_files, output_file=None):
    """
    Create publication-quality plot of <n_subjet> vs eta for multiple energies.
    
    input_files: list of filenames
    """
    
    # Create figure
    fig, ax = plt.subplots(figsize=(10, 10))
    
    # Line styles for different energies (order: 300, 141, 105, 64)
    linestyles_dict = {
        300: (0, ()),               # Solid
        141: (0, (5, 5)),           # Long dash
        105: (0, (3, 5, 1, 5)),     # Dash-dot
        64:  (0, (1, 1))            # Dotted (tight)
    }
    
    markers_dict = {
        300: 'o',
        141: 's',
        105: '^',
        64: 'd'
    }
    
    linewidth = 3
    markersize = 14
    
    # Store legend handles
    legend_handles = []
    
    for filename in input_files:
        # Read data
        data = read_mean_nsubjet_data(filename)
        meta = data['metadata']
        
        # Get energy
        sqrt_s = int(meta.get('sqrt_s', 300))
        
        # Get line style and marker
        linestyle = linestyles_dict.get(sqrt_s, '-')
        marker = markers_dict.get(sqrt_s, 'o')
        
        # Create smooth curves using spline interpolation
        eta_smooth = np.linspace(data['eta_center'].min(), data['eta_center'].max(), 300)
        # eta_smooth = np.linspace(-1.0, 2.5, 300)
        
        # Spline for QQ (k=2 for quadratic, safer with 4 points)
        spline_qq = make_interp_spline(data['eta_center'], data['mean_qq'], k=2)
        mean_qq_smooth = spline_qq(eta_smooth)
        
        # Spline for GG
        spline_gg = make_interp_spline(data['eta_center'], data['mean_gg'], k=2)
        mean_gg_smooth = spline_gg(eta_smooth)
        
        # Plot smooth curves
        line_qq, = ax.plot(eta_smooth, mean_qq_smooth, 
                           color=color_QQ, linestyle=linestyle, linewidth=linewidth)
        line_gg, = ax.plot(eta_smooth, mean_gg_smooth, 
                           color=color_GG, linestyle=linestyle, linewidth=linewidth)
        
        # Plot data points on top
        ax.scatter(data['eta_center'], data['mean_qq'], 
                   color=color_QQ, s=markersize**2, marker=marker,
                   edgecolor='white', linewidth=1.5, zorder=5)
        ax.scatter(data['eta_center'], data['mean_gg'], 
                   color=color_GG, s=markersize**2, marker=marker,
                   edgecolor='white', linewidth=1.5, zorder=5)
        
        # Create legend entries for this energy
        legend_qq = Line2D([0], [0], color=color_QQ, linestyle=linestyle, 
                           linewidth=linewidth, marker=marker, markersize=10,
                           markerfacecolor=color_QQ, markeredgecolor='white',
                           label=rf'\textbf{{QQ ep, {sqrt_s} GeV}}')
        legend_gg = Line2D([0], [0], color=color_GG, linestyle=linestyle, 
                           linewidth=linewidth, marker=marker, markersize=10,
                           markerfacecolor=color_GG, markeredgecolor='white',
                           label=rf'\textbf{{GG ep, {sqrt_s} GeV}}')
        
        legend_handles.append(legend_qq)
        legend_handles.append(legend_gg)
    
    # ============================================
    # LABELS
    # ============================================
    
    ax.set_xlabel(r'$\mathbf{\eta_{jet}}$')
    ax.set_ylabel(r'$\mathbf{\langle n_{sbj} \rangle}$')
    
    # Axis limits
    ax.set_xlim(-1.25, 2.5)
    ax.set_ylim(1, 8)
    
    # ============================================
    # LEGEND
    # ============================================
    
    ax.legend(handles=legend_handles, loc='upper right', frameon=False, fontsize=28, ncol=1, handlelength=3)
    
    # ============================================
    # TEXT ANNOTATIONS
    # ============================================
    
    # # Add ycut info
    # ycut = meta.get('ycut', 0.0005)
    # ax.text(0.05, 0.95, rf'$\mathbf{{y_{{cut}} = {ycut}}}$', 
    #         transform=ax.transAxes, fontsize=22, verticalalignment='top')
    
    # ============================================
    # TICKS
    # ============================================
    
    ax.minorticks_on()
    ax.tick_params(which='both', direction='in', top=True, right=True)
    
    # Set the fixed x-axis limits and ticks
    ax.set_xlim(-1.25, 2.5)
    # ax.set_xticks([-1, -0.5, 0, 0.5, 1, 1.5, 2])
    ax.set_xticks([-1, 0, 1, 2])
    
    plt.tight_layout()
    
    # Save
    if output_file is None:
        output_file = 'mean_nsubjet_vs_eta_comparison.pdf'
    
    plt.savefig(output_file, bbox_inches='tight', pad_inches=0.1)
    print(f'Plot saved: {output_file}')
    plt.show()


def plot_single_file(input_file, output_file=None):
    """Simple wrapper to plot a single file with spline interpolation."""
    
    # Read data
    data = read_mean_nsubjet_data(input_file)
    meta = data['metadata']
    
    # Determine output filename
    if output_file is None:
        base_name = os.path.splitext(os.path.basename(input_file))[0]
        output_file = f'{base_name}.pdf'
    
    # Create figure
    fig, ax = plt.subplots(figsize=(16, 16))
    
    linewidth = 3
    markersize = 10
    
    # Create smooth curves using spline interpolation
    eta_smooth = np.linspace(data['eta_center'].min(), data['eta_center'].max(), 300)
    
    # Spline for QQ (k=2 for quadratic, safer with 4 points)
    spline_qq = make_interp_spline(data['eta_center'], data['mean_qq'], k=2)
    mean_qq_smooth = spline_qq(eta_smooth)
    
    # Spline for GG
    spline_gg = make_interp_spline(data['eta_center'], data['mean_gg'], k=2)
    mean_gg_smooth = spline_gg(eta_smooth)
    
    # Plot smooth curves
    ax.plot(eta_smooth, mean_qq_smooth, 
            color=color_QQ, linestyle='-', linewidth=linewidth,
            label=r'\textbf{QQ}')
    ax.plot(eta_smooth, mean_gg_smooth, 
            color=color_GG, linestyle='-', linewidth=linewidth,
            label=r'\textbf{GG}')
    
    # Plot data points on top
    ax.scatter(data['eta_center'], data['mean_qq'], 
               color=color_QQ, s=markersize**2, marker='o',
               edgecolor='white', linewidth=2, zorder=5)
    ax.scatter(data['eta_center'], data['mean_gg'], 
               color=color_GG, s=markersize**2, marker='s',
               edgecolor='white', linewidth=2, zorder=5)
    
    # Labels
    ax.set_xlabel(r'$\mathbf{\eta_{jet}}$')
    ax.set_ylabel(r'$\mathbf{\langle n_{subjet} \rangle}$')
    
    # Axis limits
    ax.set_xlim(-1.0, 3.0)
    y_max = max(np.max(data['mean_qq']), np.max(data['mean_gg']))
    ax.set_ylim(0, y_max * 1.2)
    
    # Legend
    ax.legend(loc='upper right', frameon=False, fontsize=24)
    
    # Text annotations
    sqrt_s = meta.get('sqrt_s', 300)
    ycut = meta.get('ycut', 0.0005)
    
    ax.text(0.05, 0.95, rf'\textbf{{ep, {int(sqrt_s)} GeV}}', transform=ax.transAxes, fontsize=30, verticalalignment='top')
    # ax.text(0.05, 0.88, rf'$\mathbf{{y_{{cut}} = {ycut}}}$', transform=ax.transAxes, fontsize=22, verticalalignment='top')
    
    # Ticks
    ax.minorticks_on()
    ax.tick_params(which='both', direction='in', top=True, right=True)
    ax.set_xticks([-1,-0.5,0,0.5,1,1.5,2,2.5])
    
    plt.tight_layout()
    plt.savefig(output_file, bbox_inches='tight', pad_inches=0.15)
    print(f'Plot saved: {output_file}')
    plt.show()


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage:")
        print("  Single file:    python plot_mean_nsubjet_vs_eta.py <input_file.txt>")
        print("  Multiple files: python plot_mean_nsubjet_vs_eta.py file1.txt file2.txt ...")
        print("")
        print("Examples:")
        print("  python plot_mean_nsubjet_vs_eta.py hera300_pT7_10_7.txt")
        print("  python plot_mean_nsubjet_vs_eta.py hera300_pT7_10_7.txt eic141_pT5_10_7.txt")
        sys.exit(1)
    
    if len(sys.argv) == 2:
        # Single file mode
        plot_single_file(sys.argv[1])
    else:
        # Multiple files mode
        input_files = sys.argv[1:]
        plot_multiple_energies(input_files, 'meannsbjetaminus1to2.pdf')