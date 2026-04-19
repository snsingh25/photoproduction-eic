#!/usr/bin/env python3
"""
Jet E_T vs eta analysis for ep collisions
Comparing QQ and GG processes at different center-of-mass energies
"""

import ROOT
import numpy as np
from array import array
import sys

def setup_root_style():
    """Configure ROOT plotting style for publication quality"""
    
    # Enable LaTeX fonts and set publication quality style
    ROOT.gROOT.SetStyle("Plain")
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    
    # LaTeX font settings (42 = Helvetica, closest to Computer Modern in ROOT)
    ROOT.gStyle.SetTextFont(42)
    ROOT.gStyle.SetLabelFont(42, "XYZ")
    ROOT.gStyle.SetTitleFont(42, "XYZ")
    ROOT.gStyle.SetStatFont(42)
    ROOT.gStyle.SetTextSize(0.04)
    ROOT.gStyle.SetLabelSize(0.05, "XYZ")
    ROOT.gStyle.SetTitleSize(0.05, "XYZ")
    
    # Use ROOT's default jet colormap (rainbow palette)
    ROOT.gStyle.SetPalette(1)  # Rainbow palette
    ROOT.gStyle.SetNumberContours(100)
    
    # Clean appearance
    ROOT.gStyle.SetFrameBorderMode(0)
    ROOT.gStyle.SetCanvasBorderMode(0)
    ROOT.gStyle.SetPadBorderMode(0)
    ROOT.gStyle.SetPadColor(0)
    ROOT.gStyle.SetCanvasColor(0)
    ROOT.gStyle.SetFrameFillColor(0)


def create_dataset_config():
    """Define dataset configuration"""
    datasets = [
        {
            'filepath': '/Users/siddharthsingh/Analysis/ph-new/subprocjets/jets-basic/dijets/hera300_dijets_pT7/dijets_hera300_R10_EtMin0_10_7.root',
            'label': 'ep, 300 GeV',
            'shortname': 'hera300'
        },
        {
            'filepath': '/Users/siddharthsingh/Analysis/ph-new/subprocjets/jets-basic/dijets/eic141_dijets_pT5/dijets_eic141_R10_EtMin0_10_7.root',
            'label': 'ep, 141 GeV',
            'shortname': 'eic141'
        },
        {
            'filepath': '/Users/siddharthsingh/Analysis/ph-new/subprocjets/jets-basic/dijets/eic105_dijets_pT5/dijets_eic105_R10_EtMin0_10_7.root',
            'label': 'ep, 105 GeV',
            'shortname': 'eic105'
        },
        {
            'filepath': '/Users/siddharthsingh/Analysis/ph-new/subprocjets/jets-basic/dijets/eic64_dijets_pT5/dijets_eic64_R10_EtMin0_10_7.root',
            'label': 'ep, 64 GeV',
            'shortname': 'eic64'
        }
    ]
    
    event_types = ['QQ_Events', 'GG_Events']
    event_labels = ['QQ', 'GG']
    
    return datasets, event_types, event_labels


def get_color_scale_limits():
    """Define individual color scale limits for each histogram"""
    color_limits = [
        (0, 1200),  # HERA 300 QQ
        (0, 900),   # HERA 300 GG
        (0, 500),   # EIC 141 QQ
        (0, 300),   # EIC 141 GG
        (0, 500),   # EIC 105 QQ
        (0, 180),   # EIC 105 GG
        (0, 500),   # EIC 64 QQ
        (0, 80)     # EIC 64 GG
    ]
    return color_limits


def create_canvas_and_pads():
    """Create canvas with 4x2 grid of pads"""
    
    # Create canvas with precise dimensions
    canvas = ROOT.TCanvas("canvas", "Jet E_{T} vs #eta Analysis", 1200, 1200)
    canvas.SetFillColor(0)
    canvas.SetBorderMode(0)
    canvas.SetFrameFillColor(0)
    canvas.SetFrameBorderMode(0)
    
    # Layout parameters
    plot_width = 0.32    # Width for each plot
    cbar_width = 0.001   # Width for colorbar space
    plot_height = 0.22
    left_margin = 0.08
    top_margin = 0.04
    v_spacing = 0.0
    
    # Create 8 plotting pads (4 rows x 2 columns)
    pads = []
    pad_count = 0
    
    for row in range(4):
        for plot_col in range(2):  # Only 2 plotting columns
            
            # Calculate pad positions
            if plot_col == 0:  # Left plot
                x1 = left_margin
                x2 = x1 + plot_width
            else:  # Right plot
                x1 = left_margin + plot_width + cbar_width
                x2 = x1 + plot_width
            
            y1 = 1.0 - top_margin - (row + 1) * plot_height - row * v_spacing
            y2 = y1 + plot_height
            
            # Create pad
            pad_name = f"pad_{pad_count}"
            pad = ROOT.TPad(pad_name, "", x1, y1, x2, y2)
            pad.SetFillColor(0)
            pad.SetBorderMode(0)
            pad.SetFrameBorderMode(0)
            
            # Set margins for individual colorbars
            pad.SetLeftMargin(0.06)
            pad.SetRightMargin(0.15)  # Space for individual colorbars
            pad.SetBottomMargin(0.12 if row == 3 else 0.01)
            pad.SetTopMargin(0.01)
            
            # Enable ticks and linear scale
            pad.SetTickx(1)
            pad.SetTicky(1)
            pad.SetLogz(0)  # Linear scale for colors
            
            # Draw pad on canvas
            canvas.cd()
            pad.Draw()
            
            pads.append(pad)
            pad_count += 1
    
    return canvas, pads


def process_dataset(file_path, event_type, event_label, dataset_label, shortname):
    """Process a single dataset and create histogram"""
    
    print(f"  Processing {event_label} events...")
    
    # Open ROOT file
    root_file = ROOT.TFile.Open(file_path, "READ")
    if not root_file or root_file.IsZombie():
        print(f"Error: Cannot open file {file_path}")
        return None, None
    
    # Navigate to event directory
    event_dir = root_file.Get(event_type)
    if not event_dir:
        print(f"Error: Cannot find directory {event_type}")
        root_file.Close()
        return None, None
    
    # Get tree
    tree_name = f"jets_{event_type}"
    tree = event_dir.Get(tree_name)
    if not tree:
        print(f"Error: Cannot find tree {tree_name}")
        root_file.Close()
        return None, None
    
    # Create histogram
    hist_name = f"h_{shortname}_{event_label}"
    hist_title = f"{dataset_label} {event_label};#eta_{{jet}};E_{{T}} [GeV]"
    
    hist = ROOT.TH2F(hist_name, hist_title,
                      60, -1.5, 4.0,   # eta range (60 bins)
                      60, 5.0, 70.0)   # ET range (60 bins)
    hist.SetDirectory(0)  # Detach from file
    
    # Check if branches exist
    if not tree.GetBranch("jet_et") or not tree.GetBranch("jet_eta"):
        print("Error: Required branches not found")
        root_file.Close()
        return None, None
    
    # Process events
    n_entries = tree.GetEntries()
    print(f"    Processing {n_entries} events for {event_label}")
    
    # Statistics variables
    total_jets = 0
    sum_eta = 0.0
    sum_et = 0.0
    sum_eta2 = 0.0
    sum_et2 = 0.0
    
    # Process tree
    for i in range(n_entries):
        if i % 100000 == 0 and i > 0:
            print(f"      Event {i}/{n_entries} ({100.0 * i / n_entries:.1f}%)")
        
        tree.GetEntry(i)
        
        # Access jet vectors
        jet_et = tree.jet_et
        jet_eta = tree.jet_eta
        
        if not jet_et or not jet_eta or jet_et.size() != jet_eta.size():
            continue
        
        # Fill histogram and calculate statistics
        for j in range(jet_et.size()):
            eta_val = jet_eta[j]
            et_val = jet_et[j]
            
            # Check for finite values
            if np.isfinite(eta_val) and np.isfinite(et_val):
                hist.Fill(eta_val, et_val)
                
                # Update statistics
                sum_eta += eta_val
                sum_et += et_val
                sum_eta2 += eta_val * eta_val
                sum_et2 += et_val * et_val
                total_jets += 1
    
    # Calculate statistics
    stats = {
        'n_jets': total_jets,
        'mean_eta': 0.0,
        'mean_et': 0.0,
        'std_eta': 0.0,
        'std_et': 0.0
    }
    
    if total_jets > 0:
        stats['mean_eta'] = sum_eta / total_jets
        stats['mean_et'] = sum_et / total_jets
        stats['std_eta'] = np.sqrt(sum_eta2 / total_jets - stats['mean_eta']**2)
        stats['std_et'] = np.sqrt(sum_et2 / total_jets - stats['mean_et']**2)
    
    print(f"    Total jets: {total_jets}")
    
    root_file.Close()
    
    return hist, stats


def style_histogram(hist, hist_idx, row, col):
    """Apply styling to histogram"""
    
    # Set axis ranges
    hist.GetXaxis().SetRangeUser(-1.5, 4.0)
    hist.GetYaxis().SetRangeUser(5.0, 70.0)
    
    # Style the axes
    hist.GetXaxis().SetTitleFont(42)
    hist.GetXaxis().SetLabelFont(42)
    hist.GetYaxis().SetTitleFont(42)
    hist.GetYaxis().SetLabelFont(42)
    
    # Remove individual axis titles
    hist.GetXaxis().SetTitle("")
    hist.GetYaxis().SetTitle("")
    
    # Configure tick sizes
    hist.GetXaxis().SetTickLength(0.05)
    hist.GetYaxis().SetTickLength(0.03)
    
    # Show Y-axis values
    hist.GetYaxis().SetLabelSize(0.08)
    hist.GetYaxis().SetNdivisions(506)  # 5 major divisions, 6 minor
    
    # Show X-axis values only on bottom row
    if row == 3:
        hist.GetXaxis().SetLabelSize(0.08)
    else:
        hist.GetXaxis().SetLabelSize(0)


def add_labels(pad, dataset_label, event_label):
    """Add dataset and process labels to pad"""
    
    # Dataset label
    label1 = ROOT.TLatex()
    label1.SetNDC()
    label1.SetTextFont(42)
    label1.SetTextSize(0.08)
    label1.SetTextColor(ROOT.kBlack)
    label1.SetTextAlign(13)  # Right top alignment
    label1.DrawLatex(0.12, 0.90, dataset_label)
    
    # Event type label
    label2 = ROOT.TLatex()
    label2.SetNDC()
    label2.SetTextFont(42)
    label2.SetTextSize(0.08)
    label2.SetTextColor(ROOT.kBlack)
    label2.SetTextAlign(13)  # Right top alignment
    label2.DrawLatex(0.12, 0.80, event_label)
    
    return label1, label2


def add_axis_titles(canvas):
    """Add common axis titles to canvas"""
    
    canvas.cd()
    
    # Y-axis title
    y_title = ROOT.TLatex()
    y_title.SetNDC()
    y_title.SetTextFont(42)
    y_title.SetTextSize(0.035)
    y_title.SetTextAlign(22)
    y_title.SetTextAngle(90)
    y_title.DrawLatex(0.05, 0.5, "E_{T} [GeV]")
    
    # X-axis title
    x_title = ROOT.TLatex()
    x_title.SetNDC()
    x_title.SetTextFont(42)
    x_title.SetTextSize(0.035)
    x_title.SetTextAlign(22)
    x_title.DrawLatex(0.4, 0.05, "#eta_{jet}")
    
    return y_title, x_title


def print_statistics(datasets, event_labels, histograms, all_stats):
    """Print analysis summary"""
    
    print("\n" + "=" * 60)
    print("ANALYSIS SUMMARY")
    print("=" * 60)
    
    for file_idx, dataset in enumerate(datasets):
        print(f"\n{dataset['label']}:")
        
        for event_idx, event_label in enumerate(event_labels):
            hist_idx = file_idx * 2 + event_idx
            
            if histograms[hist_idx]:
                stats = all_stats[hist_idx]
                print(f"  {event_label}: {stats['n_jets']} jets, "
                      f"<eta> = {stats['mean_eta']:.2f}±{stats['std_eta']:.2f}, "
                      f"<E_T> = {stats['mean_et']:.1f}±{stats['std_et']:.1f} GeV")
            else:
                print(f"  {event_label}: No data available")


def main():
    """Main analysis function"""
    
    # Setup ROOT style
    setup_root_style()
    
    # Get configuration
    datasets, event_types, event_labels = create_dataset_config()
    color_limits = get_color_scale_limits()
    
    # Create canvas and pads
    canvas, pads = create_canvas_and_pads()
    
    # Storage for histograms and statistics
    histograms = [None] * 8
    all_stats = [None] * 8
    labels_list = []  # Keep references to prevent garbage collection
    
    # Process each dataset
    for file_idx, dataset in enumerate(datasets):
        print(f"Processing {dataset['label']}...")
        
        for event_idx, (event_type, event_label) in enumerate(zip(event_types, event_labels)):
            hist_idx = file_idx * 2 + event_idx
            
            # Process dataset
            hist, stats = process_dataset(
                dataset['filepath'],
                event_type,
                event_label,
                dataset['label'],
                dataset['shortname']
            )
            
            histograms[hist_idx] = hist
            all_stats[hist_idx] = stats
    
    # Draw histograms
    for hist_idx, hist in enumerate(histograms):
        if not hist:
            continue
        
        # Determine row and column
        row = hist_idx // 2
        col = hist_idx % 2
        
        # Switch to pad
        pads[hist_idx].cd()
        
        # Apply styling
        style_histogram(hist, hist_idx, row, col)
        
        # Set color scale limits
        zmin, zmax = color_limits[hist_idx]
        hist.SetMinimum(zmin)
        hist.SetMaximum(zmax)
        
        # Adjust Z-axis range for colorbar
        hist.GetZaxis().SetRangeUser(zmin + 1, zmax - 1)
        
        # Draw histogram with colorbar
        hist.Draw("COLZ")
        
        # Adjust colorbar font size
        pads[hist_idx].Update()  # Force update to create palette
        palette = hist.GetListOfFunctions().FindObject("palette")
        if palette:
            palette.SetLabelSize(0.06)
        
        # Add labels
        dataset_label = datasets[row]['label']
        event_label = event_labels[col]
        labels = add_labels(pads[hist_idx], dataset_label, event_label)
        labels_list.extend(labels)
        
        # Update pad
        pads[hist_idx].Modified()
        pads[hist_idx].Update()
    
    # Add common axis titles
    axis_labels = add_axis_titles(canvas)
    labels_list.extend(axis_labels)
    
    # Update canvas
    canvas.Update()
    
    # Save PDF
    output_file = "jet_et_vs_eta_python.pdf"
    canvas.SaveAs(output_file)
    print(f"\nSaved: {output_file}")
    
    # Print statistics
    print_statistics(datasets, event_labels, histograms, all_stats)
    
    print("\nAnalysis completed successfully!")
    
    # Keep canvas alive (for interactive mode)
    return canvas, pads, histograms, labels_list


if __name__ == "__main__":
    # Run analysis
    canvas, pads, histograms, labels = main()
    
    # If running interactively, keep the window open
    if hasattr(ROOT, 'gApplication'):
        ROOT.gApplication.Run()