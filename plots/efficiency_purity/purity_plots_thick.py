#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rcParams
from scipy.interpolate import make_interp_spline

rcParams['font.family'] = 'serif'
rcParams['font.serif'] = ['Computer Modern Roman']
rcParams['text.usetex'] = True
rcParams['text.latex.preamble'] = r'\usepackage{amsmath} \usepackage{amssymb}'
rcParams['font.size'] = 18
rcParams['axes.labelsize'] = 20
rcParams['axes.titlesize'] = 22
rcParams['xtick.labelsize'] = 22
rcParams['ytick.labelsize'] = 22
# rcParams['legend.fontsize'] = 18
# rcParams['figure.dpi'] = 100
rcParams['savefig.dpi'] = 300
rcParams['axes.linewidth'] = 2
rcParams['xtick.major.width'] = 2.0
rcParams['ytick.major.width'] = 2.0
rcParams['xtick.minor.width'] = 1.0
rcParams['ytick.minor.width'] = 1.0

def create_purity_plot():
    
    # X-axis points: center of eta ranges
    eta_centers = np.array([-0.5, 0.5, 1.5, 2.5])  # Centers of -1to0, 0to1, 1to2, 2to3
    
    # Thick jet composition (percentages) HERA
    thick_QQ = np.array([42.5, 28.4, 19.0, 14.5])
    thick_GG = np.array([8.2, 20.5, 33.2, 38.6])
    thick_GQ = np.array([49.4, 51.2, 47.8, 46.9])
    
    # ==========================================================================
    # CREATE FIGURE AND AXES
    # ==========================================================================
    
    fig, ax = plt.subplots(figsize=(9, 9))
    
    # Set minimalistic style
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    
    # ==========================================================================
    # CREATE SMOOTH INTERPOLATION
    # ==========================================================================

    eta_smooth = np.linspace(eta_centers.min(), eta_centers.max(), 300)  # Stay within data bounds

    # Create spline interpolations (k=3 for cubic splines)
    spline_QQ = make_interp_spline(eta_centers, thick_QQ, k=3)
    spline_GG = make_interp_spline(eta_centers, thick_GG, k=3)
    spline_GQ = make_interp_spline(eta_centers, thick_GQ, k=3)

    # Generate smooth curves
    thick_QQ_smooth = spline_QQ(eta_smooth)
    thick_GG_smooth = spline_GG(eta_smooth)
    thick_GQ_smooth = spline_GQ(eta_smooth)
    
    # ==========================================================================
    # PLOT DATA
    # ==========================================================================
    
    # Define colors matching ROOT style
    color_QQ = '#CC0000'  # Red
    color_GG = '#00AA00'  # Green  
    color_GQ = '#0066CC'  # Blue
    
    # Plot smooth curves
    ax.plot(eta_smooth, thick_QQ_smooth, color=color_QQ, linewidth=4, label=r'$\mathrm{QQ}$',zorder=2)
    ax.plot(eta_smooth, thick_GG_smooth, color=color_GG, linewidth=4, label=r'$\mathrm{GG}$',zorder=2)
    ax.plot(eta_smooth, thick_GQ_smooth, color=color_GQ, linewidth=4, label=r'$\mathrm{GQ}$',zorder=2)
    
    # Add the original data points as markers
    ax.scatter(eta_centers, thick_QQ, 
              color=color_QQ, 
              s=150, 
              marker='o',
              edgecolors='black',
              linewidth=1,
              zorder=3)
    
    ax.scatter(eta_centers, thick_GG, 
              color=color_GG, 
              s=150, 
              marker='s',
              edgecolors='black',
              linewidth=1,
              zorder=3)
    
    ax.scatter(eta_centers, thick_GQ, 
              color=color_GQ, 
              s=150, 
              marker='^',
              edgecolors='black',
              linewidth=1,
              zorder=3)
    
    # ==========================================================================
    # AXIS CONFIGURATION
    # ==========================================================================
    
    # Set axis labels with LaTeX
    ax.set_xlabel(r'\textbf{Pseudorapidity} $\boldsymbol{\eta}$', fontsize=26)
    ax.set_ylabel(r'\textbf{Percentage (\%)}', fontsize=26)
    
    # Set axis ranges
    # ax.set_xlim(-1.2, 3.2)
    ax.set_xlim(-1.0, 3.0)
    ax.set_ylim(0, 100)
    
    # Set custom x-ticks to match eta bins
    ax.set_xticks([-1, 0, 1, 2, 3])
    ax.set_xticklabels([r'$-1$', r'$0$', r'$1$', r'$2$', r'$3$'])
    
    # Set y-ticks
    ax.set_yticks(np.arange(0, 101, 10))
    
    # Add minor ticks for better readability
    ax.minorticks_on()
    ax.tick_params(which='minor', length=5, width=1.0)
    ax.tick_params(which='major', length=8, width=1.5)
    
    # Add grid for better readability (optional - comment out for more minimal look)
    ax.grid(True, alpha=0.25, linestyle='--', linewidth=0.5, which='major')
    
    # ==========================================================================
    # ADD LEGEND (MINIMALISTIC)
    # ==========================================================================
    
    legend = ax.legend(loc='upper right', 
                      title=r'\textbf{Thick Jets' + '\n' + r'\textbf{ep 300 GeV}',
                      frameon=False,
                      numpoints=1,
                      handlelength=2.5,
                      columnspacing=1,
                      fontsize=26,
                      title_fontsize=26)
    
    # ==========================================================================
    # FINAL ADJUSTMENTS
    # ==========================================================================
    
    plt.tight_layout(pad=1.5)
    
    # Save figures
    plt.savefig('purity_thick_hera.pdf', bbox_inches='tight')
    
    return fig, ax

def main():
    fig, ax = create_purity_plot()
    plt.show()

if __name__ == "__main__":
    main()