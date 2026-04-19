#!/usr/bin/env python3
"""
Differential Jet Shapes Plotting Script
Pure matplotlib version (no ROOT required)
"""

import matplotlib.pyplot as plt
import numpy as np

def setup_physics_style():
    """Setup matplotlib to look like professional physics plots"""
    
    plt.rcParams.update({
        'font.size': 12,
        'font.family': 'serif',
        'font.serif': ['Times New Roman', 'Times', 'DejaVu Serif'],
        'mathtext.fontset': 'stix',  # For math symbols
        'axes.linewidth': 1.2,
        'xtick.major.size': 6,
        'xtick.minor.size': 3,
        'ytick.major.size': 6,
        'ytick.minor.size': 3,
        'xtick.direction': 'in',
        'ytick.direction': 'in',
        'xtick.top': True,
        'ytick.right': True,
        'legend.frameon': True,
        'legend.fancybox': False,
        'legend.edgecolor': 'black',
        'legend.framealpha': 1.0,
        'axes.grid': False,
        'figure.dpi': 150,
        'savefig.dpi': 300,
        'savefig.bbox': 'tight'
    })

def plot_differential_jet_subplots():
    """Create differential jet shape plots with multiple energy panels"""
    
    # Setup professional styling
    setup_physics_style()
    
    # Define r values
    r_bins = np.array([0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0])
    r_centers = (r_bins[:-1] + r_bins[1:]) / 2  # Bin centers for plotting
    
    # EIC 64 GeV data: -4 < eta < 4, ET_j > 10GeV
    quark_values = [4.0742, 3.0304, 1.8749, 1.1804, 0.7980, 0.5787, 0.4420, 0.3400, 0.2134, 0.0288]
    gluon_values = [1.8683, 2.2251, 1.8797, 1.5105, 1.1599, 0.9186, 0.7097, 0.5976, 0.3276, 0.0482]
    gluon_quark_values = [3.3582, 2.8521, 1.9543, 1.3377, 0.9575, 0.7191, 0.5425, 0.4236, 0.2703, 0.0368]
    
    # EIC 105 GeV data
    quark_values_2 = [4.2309, 3.1495, 1.8968, 1.1836, 0.7812, 0.5813, 0.4403, 0.3446, 0.2254, 0.0309]
    gluon_values_2 = [1.9396, 2.3131, 1.9770, 1.6216, 1.2180, 0.9463, 0.7692, 0.5830, 0.3885, 0.0570]
    gluon_quark_values_2 = [3.4005, 2.8527, 1.9469, 1.3452, 0.9682, 0.7254, 0.5670, 0.4470, 0.2935, 0.0418]
    
    # EIC 141 GeV data
    quark_values_3 = [4.3648, 3.2278, 1.9475, 1.2198, 0.8301, 0.6095, 0.4731, 0.3819, 0.2572, 0.0358]
    gluon_values_3 = [2.0503, 2.5228, 2.1952, 1.7422, 1.3845, 1.1140, 0.9263, 0.7585, 0.5138, 0.0828]
    gluon_quark_values_3 = [3.3891, 2.9446, 2.0608, 1.4483, 1.0586, 0.8130, 0.6512, 0.5243, 0.3566, 0.0528]
    
    # HERA 300 GeV data: ET_j > 17GeV
    quark_values_4 = [5.2709, 2.7920, 1.5292, 0.9040, 0.6041, 0.4580, 0.3398, 0.2784, 0.2252, 0.0410]
    gluon_values_4 = [2.5964, 2.4778, 1.8710, 1.3170, 1.0125, 0.7854, 0.6305, 0.5161, 0.3929, 0.0744]
    gluon_quark_values_4 = [4.1289, 2.7390, 1.7630, 1.1515, 0.8064, 0.6157, 0.4709, 0.3915, 0.2924, 0.0550]
    
    # Organize data by panel
    data_sets = [
        (quark_values, gluon_values, gluon_quark_values),
        (quark_values_2, gluon_values_2, gluon_quark_values_2),
        (quark_values_3, gluon_values_3, gluon_quark_values_3),
        (quark_values_4, gluon_values_4, gluon_quark_values_4)
    ]
    
    # Panel-specific info
    panel_info = [
        ("ep (64 GeV)", r"$E_T^{\mathrm{jet}} > 10$ GeV", r"$-4 < \eta_{\mathrm{jet}} < 4$"),
        ("ep (105 GeV)", r"$E_T^{\mathrm{jet}} > 10$ GeV", r"$-4 < \eta_{\mathrm{jet}} < 4$"),
        ("ep (141 GeV)", r"$E_T^{\mathrm{jet}} > 10$ GeV", r"$-4 < \eta_{\mathrm{jet}} < 4$"),
        ("ep (300 GeV)", r"$E_T^{\mathrm{jet}} > 17$ GeV", r"$-4 < \eta_{\mathrm{jet}} < 4$")
    ]
    
    # Modern color palette
    colors = {
        'gq': '#264653',      # Dark Teal
        'quark': '#2A9D8F',   # Teal
        'gluon': '#E76F51'    # Terracotta
    }
    
    # Create figure with subplots
    fig, axes = plt.subplots(2, 2, figsize=(12, 12))
    fig.subplots_adjust(hspace=0.05, wspace=0.05)
    
    # Process each subplot
    for i, ax in enumerate(axes.flat):
        
        # Get data for this panel
        quark_data, gluon_data, gq_data = data_sets[i]
        energy_label, et_label, eta_label = panel_info[i]
        
        # Create step plots (histogram style)
        ax.step(r_bins, [quark_data[0]] + quark_data, where='pre', 
                color=colors['quark'], linewidth=2, label='Quark')
        ax.step(r_bins, [gluon_data[0]] + gluon_data, where='pre', 
                color=colors['gluon'], linewidth=2, label='Gluon')
        ax.step(r_bins, [gq_data[0]] + gq_data, where='pre', 
                color=colors['gq'], linewidth=2, label='Quark+Gluon')
        
        # Set axis limits and styling
        ax.set_xlim(0, 1.0)
        ax.set_ylim(0, 5.5)
        
        # Configure ticks
        ax.set_xticks(np.arange(0, 1.2, 0.2))
        ax.tick_params(which='both', direction='in', top=True, right=True)
        ax.tick_params(which='major', length=6, width=1.2)
        ax.tick_params(which='minor', length=3, width=1)
        
        # Panel-specific configurations
        if i == 0:  # Top left
            ax.set_xticklabels([])  # Hide x labels
            ax.set_ylabel(r'$\rho$', fontsize=16)
            # Hide "0" on y-axis
            yticks = ax.get_yticks()
            yticklabels = ['' if tick == 0 else f'{tick:.0f}' for tick in yticks]
            ax.set_yticklabels(yticklabels)
            
        elif i == 1:  # Top right
            ax.set_xticklabels([])  # Hide x labels
            ax.set_yticklabels([])  # Hide y labels
            
            # Add legend
            legend = ax.legend(loc='upper right', bbox_to_anchor=(0.98, 0.98), 
                             frameon=True, fancybox=False, edgecolor='black',
                             fontsize=11)
            legend.get_frame().set_linewidth(1)
            
        elif i == 2:  # Bottom left
            ax.set_xlabel(r'$r$', fontsize=16)
            ax.set_ylabel(r'$\rho$', fontsize=16)
            # Custom x-axis: hide last tick (1.0)
            xticks = ax.get_xticks()
            xticklabels = [f'{tick:.1f}' if tick != 1.0 else '' for tick in xticks]
            ax.set_xticklabels(xticklabels)
            
        elif i == 3:  # Bottom right
            ax.set_xlabel(r'$r$', fontsize=16)
            ax.set_yticklabels([])  # Hide y labels
            # Custom x-axis: hide first tick (0.0)
            xticks = ax.get_xticks()
            xticklabels = ['' if tick == 0.0 else f'{tick:.1f}' for tick in xticks]
            ax.set_xticklabels(xticklabels)
        
        # Add text annotations
        ax.text(0.97, 0.95, energy_label, transform=ax.transAxes, 
                fontsize=12, ha='right', va='top', weight='bold')
        ax.text(0.97, 0.35, et_label, transform=ax.transAxes, 
                fontsize=10, ha='right', va='top')
        ax.text(0.97, 0.28, eta_label, transform=ax.transAxes, 
                fontsize=10, ha='right', va='top')
    
    # Add shared axis labels
    fig.text(0.5, 0.02, r'$r$', ha='center', va='bottom', fontsize=18)
    fig.text(0.02, 0.5, r'$\rho$', ha='center', va='center', fontsize=18, rotation=90)
    
    # Save the figure
    plt.savefig('diff_jet_shapes_sub_eta_all_matplotlib.pdf', 
                bbox_inches='tight', dpi=300)
    plt.savefig('diff_jet_shapes_sub_eta_all_matplotlib.png', 
                bbox_inches='tight', dpi=300)
    
    print("✓ PDF and PNG saved successfully!")
    print("✓ Files: diff_jet_shapes_sub_eta_all_matplotlib.pdf/png")
    
    # Show the plot
    plt.show()
    
    return fig, axes

def main():
    """Main function to run the plotting"""
    fig, axes = plot_differential_jet_subplots()

if __name__ == "__main__":
    main()