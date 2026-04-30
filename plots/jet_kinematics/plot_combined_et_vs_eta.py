#!/usr/bin/env python3
"""
Jet E_T vs eta analysis - Fixed version
=====================================
Clean implementation with Computer Modern fonts
"""

import glob
import os
import sys
from pathlib import Path

import uproot
import awkward as ak
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import LinearSegmentedColormap


# Repo root = parents[2] of this file (plots/jet_kinematics/<this>.py).
REPO_ROOT = Path(__file__).resolve().parents[2]


def find_dijet_root(sample_dir):
    """Return the first dijets_*.root or alljets_*.root in <repo>/data-jets/<sample_dir>/."""
    base = REPO_ROOT / "data-jets" / sample_dir
    for pat in ("dijets_*.root", "alljets_*.root"):
        hits = sorted(glob.glob(str(base / pat)))
        if hits:
            return hits[0]
    return None

# LaTeX Computer Modern fonts
plt.rcParams.update({
    'font.size': 20,          # Increase base size (was 14)
    'axes.labelsize': 28,     # Size for x and y labels
    'xtick.labelsize': 28,    # Size for x tick numbers
    'ytick.labelsize': 28,    # Size for y tick numbers
    'font.family': 'serif',
    'font.serif': ['Computer Modern'],
    'text.usetex': True,
    'axes.linewidth': 1.5,    # Thicken frame slightly for larger text
    'xtick.major.size': 8,    # Make ticks larger to match text
    'ytick.major.size': 8,
    'xtick.direction': 'in',
    'ytick.direction': 'in',
    'xtick.top': True,
    'ytick.right': True,
    'savefig.dpi': 300
})

def create_root_like_colormap():
    """
    Create a colormap similar to ROOT's rainbow palette
    This mimics the ROOT.gStyle.SetPalette(1) behavior
    """
    # ROOT rainbow palette colors
    # colors = [
    #     (0.00, (0.7, 0.5, 0.9)),   # Brighter purple (was dark blue)
    #     (0.15, (0.2, 0.0, 0.6)),   # Medium purple 
    #     (0.25, (0.0, 0.2, 1.0)),   # Bright blue
    #     (0.40, (0.0, 1.0, 1.0)),   # Cyan
    #     (0.55, (0.0, 1.0, 0.0)),   # Green
    #     (0.70, (1.0, 1.0, 0.0)),   # Yellow
    #     (0.85, (1.0, 0.5, 0.0)),   # Orange
    #     (1.00, (1.0, 0.0, 0.0))    # Red
    # ]
    colors = [
        (0.00, (0.5, 0.0, 1.0)),    # Bright magenta/purple
        (0.15, (0.3, 0.0, 0.9)),    # Deep purple
        (0.30, (0.0, 0.0, 0.8)),    # Dark purple-blue
        (0.45, (0.0, 0.5, 1.0)),    # Light blue
        (0.60, (0.0, 1.0, 1.0)),    # Cyan
        (0.75, (0.0, 1.0, 0.0)),    # Green
        (0.85, (1.0, 1.0, 0.0)),    # Yellow
        (0.95, (1.0, 0.5, 0.0)),    # Orange
        (1.00, (1.0, 0.0, 0.0))     # Red
    ]
    
    cmap_name = 'bright_root_rainbow'
    n_bins = 256
    cmap = LinearSegmentedColormap.from_list(cmap_name, colors, N=n_bins)
    return cmap

def read_jets_from_root(filepath, event_type):
    """Read jet data from ROOT file with error handling"""
    try:
        with uproot.open(filepath) as file:
            # Check if the tree exists
            tree_name = f"{event_type}/jets_{event_type}"
            if tree_name not in file:
                print(f"Warning: Tree {tree_name} not found in {filepath}")
                # Try alternative tree structure
                available_trees = [key for key in file.keys() if event_type in key]
                if available_trees:
                    tree_name = available_trees[0]
                    print(f"Using alternative tree: {tree_name}")
                else:
                    print(f"No trees containing '{event_type}' found")
                    return np.array([]), np.array([])
            
            tree = file[tree_name]
            
            # Check if branches exist
            if "jet_et" not in tree or "jet_eta" not in tree:
                print(f"Warning: Required branches not found in {tree_name}")
                print(f"Available branches: {tree.keys()}")
                return np.array([]), np.array([])
            
            jet_et = ak.to_numpy(ak.flatten(tree["jet_et"].array(library="ak")))
            jet_eta = ak.to_numpy(ak.flatten(tree["jet_eta"].array(library="ak")))
            
            return jet_et, jet_eta
            
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return np.array([]), np.array([])

def create_2d_histogram(jet_eta, jet_et):
    """Create 2D histogram"""
    if len(jet_eta) == 0 or len(jet_et) == 0:
        # Return empty histogram
        H = np.zeros((60, 60))
        xedges = np.linspace(-1.5, 4.0, 61)
        yedges = np.linspace(5.0, 70.0, 61)
        return H, xedges, yedges
    
    H, xedges, yedges = np.histogram2d(
        jet_eta, jet_et,
        bins=[60, 60],
        range=[[-1.5, 4.0], [5.0, 70.0]]
    )
    H = np.ma.masked_where(H == 0, H)  # Mask zero entries
    return H, xedges, yedges

def check_file_existence(filepaths):
    """Check if files exist and return list of existing files"""
    existing_files = []
    for filepath in filepaths:
        if os.path.exists(filepath['filepath']):
            existing_files.append(filepath)
        else:
            print(f"Warning: File not found: {filepath['filepath']}")
    return existing_files

def main():
    """Main analysis function"""
    
    # Check for required packages
    try:
        import uproot
        import awkward
        print("Dependencies check: OK")
        print(f"  uproot version: {uproot.__version__}")
        print(f"  awkward version: {awkward.__version__}")
        print()
    except ImportError as e:
        print("ERROR: Missing required packages!")
        print(f"Error: {e}")
        return None, None, None
    
    # Discover dijet ROOT files inside this repo's data-jets/.
    datasets = [
        {'sample': 'hera300_kt_dijets',     'label': r'300 GeV'},
        {'sample': 'eic141_antikt_dijets',  'label': r'141 GeV'},
        {'sample': 'eic105_antikt_dijets',  'label': r'105 GeV'},
        {'sample': 'eic64_antikt_dijets',   'label': r'64 GeV'},
    ]
    for d in datasets:
        d['filepath'] = find_dijet_root(d['sample']) or ""
    
    # Check which files actually exist
    existing_datasets = check_file_existence(datasets)
    if not existing_datasets:
        print("ERROR: No input files found!")
        print("Please check the file paths and ensure the ROOT files are available.")
        return None, None, None
    
    event_types = ['QQ_Events', 'GG_Events']
    event_labels = ['QQ', 'GG']
    
    # Adjust for available datasets
    n_datasets = len(existing_datasets)
    
    # Create figure
    fig = plt.figure(figsize=(14, 5 * n_datasets))
    gs = gridspec.GridSpec(n_datasets, 2, figure=fig, hspace=0.0, wspace=0.001)
    
    # ROOT-like colormap
    cmap = create_root_like_colormap()
    # cmap = plt.cm.jet # turbo, viridis, plasma, inferno, magma, jet
    
    # Store histograms and statistics
    histograms = {}
    stats = {}
    
    # Process existing datasets
    for file_idx, dataset in enumerate(existing_datasets):
        for event_idx, (event_type, event_label) in enumerate(zip(event_types, event_labels)):
            hist_idx = file_idx * 2 + event_idx
            
            # Create subplot
            ax = fig.add_subplot(gs[file_idx, event_idx])
            
            # Read and plot data
            print(f"Processing {dataset['filepath']} - {event_type}")
            jet_et, jet_eta = read_jets_from_root(dataset['filepath'], event_type)
            
            # Store statistics
            stats[f"{dataset['label']}_{event_label}"] = {
                'n_jets': len(jet_et),
                'mean_et': np.mean(jet_et) if len(jet_et) > 0 else 0,
                'mean_eta': np.mean(jet_eta) if len(jet_eta) > 0 else 0
            }
            
            H, xedges, yedges = create_2d_histogram(jet_eta, jet_et)
            histograms[f"{dataset['label']}_{event_label}"] = H
            
            X, Y = np.meshgrid(xedges, yedges)

            # Adaptive color limits: 99th-percentile of populated bins keeps
            # the dynamic range high without one hot bin saturating the rest.
            populated = np.asarray(H[H > 0]).ravel() if np.ma.isMaskedArray(H) else H[H > 0]
            vmin = 0
            vmax = float(np.percentile(populated, 99)) if populated.size else 1.0
            
            im = ax.pcolormesh(X, Y, H.T, cmap=cmap, vmin=vmin, vmax=vmax, rasterized=True)
            
            # Colorbar
            cbar = plt.colorbar(im, ax=ax, pad=0.02, fraction=0.15)
            cbar.ax.tick_params(labelsize=22)
            # Get current ticks and remove the last one
            ticks = cbar.get_ticks()
            cbar.set_ticks(ticks[:-1])  # Remove last tick
            
            # Axis settings
            ax.set_xlim(-1.5, 4.0)
            ax.set_ylim(5.0, 70.0)
            
            # Labels
            ax.text(0.05, 0.90, dataset['label'], transform=ax.transAxes, fontsize=28)
            ax.text(0.05, 0.80, event_label, transform=ax.transAxes, fontsize=28)
            
            # # Add statistics
            # if len(jet_et) > 0:
            #     ax.text(0.12, 0.70, f"N={len(jet_et)}", transform=ax.transAxes, fontsize=12)
            
            # Tick settings
            if file_idx < n_datasets - 1:
                ax.set_xticklabels([])
            if event_idx > 0:
                ax.set_yticklabels([])
    
    # Axis labels
    fig.text(0.48, 0.08, r'$\eta_{\mathrm{jet}}$', fontsize=32, ha='center')
    fig.text(0.06, 0.5, r'$E_T$ (GeV)', fontsize=32, va='center', rotation=90)
    
    # Save into plots/jet_kinematics/output/ alongside the rest of the suite.
    out_dir = Path(__file__).resolve().parent / "output"
    out_dir.mkdir(parents=True, exist_ok=True)
    out_pdf = out_dir / "jeteteta2Dplot.pdf"
    plt.savefig(out_pdf, bbox_inches='tight')
    print(f"Plot saved: {out_pdf}")
    
    return fig, histograms, stats

if __name__ == "__main__":
    # Run the analysis
    try:
        result = main()
        if result[0] is not None:
            fig, histograms, stats = result
            
            # Print statistics
            print("\nAnalysis Statistics:")
            print("=" * 50)
            for key, stat in stats.items():
                print(f"{key}:")
                print(f"  Number of jets: {stat['n_jets']}")
                if stat['n_jets'] > 0:
                    print(f"  Mean E_T: {stat['mean_et']:.2f} GeV")
                    print(f"  Mean eta: {stat['mean_eta']:.2f}")
                print()
        else:
            print("Analysis failed - check file paths and data structure")
    except Exception as e:
        print(f"\nError during analysis: {e}")
        import traceback
        traceback.print_exc()