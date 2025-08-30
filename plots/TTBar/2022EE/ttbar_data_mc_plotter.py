#!/usr/bin/env python3

data_directory = "/gv0/Users/achihwan/SKNanoRunlog/out/TTbar_test/2022EE"

import ROOT
import cmsstyle as CMS
import os
import glob
from array import array
import gc

# Prevent ROOT from owning Python objects
ROOT.SetOwnership(ROOT.gROOT, False)

# Colors for plotting
DATA_COLOR = ROOT.kBlack
SIGNAL_COLOR = ROOT.TColor.GetColor("#e42536")  # Red for TTLJ signal
BACKGROUND_COLORS = [
    ROOT.TColor.GetColor("#5790fc"),  # Blue
    ROOT.TColor.GetColor("#f89c20"),  # Orange  
    ROOT.TColor.GetColor("#964a8b"),  # Purple
    ROOT.TColor.GetColor("#CDDC39"),  # Lime
    ROOT.TColor.GetColor("#009688"),  # Teal
    ROOT.TColor.GetColor("#795548"),  # Brown
    ROOT.TColor.GetColor("#9c9ca1")   # Gray for "Others"
]

class SignalBackgroundCanvas():
    def __init__(self, combined_data_hist, signal_hists, background_hists, config):
        super().__init__()
        
        self.data_hist = combined_data_hist
        self.signal_hists = signal_hists
        self.background_hists = background_hists
        self.config = config
        
        # Keep references to prevent garbage collection
        self._objects_to_keep = []
        
        # Build background stack (now includes signal)
        self._build_background_stack()
        
        # Style the data histogram
        if self.data_hist:
            self.data_hist.SetMarkerStyle(20)
            self.data_hist.SetMarkerSize(1.0)
            self.data_hist.SetMarkerColor(DATA_COLOR)
            self.data_hist.SetLineColor(DATA_COLOR)
            self.data_hist.SetLineWidth(2)
            self.data_hist.SetStats(0)
            if not self.data_hist.GetSumw2N():
                self.data_hist.Sumw2()
        
        # Create ratio histogram (Data vs Total Stack)
        self.ratio = None
        if self.data_hist and self.total_background:
            self.ratio = self.data_hist.Clone("data_stack_ratio")
            self.ratio.SetDirectory(0)
            self.ratio.SetStats(0)
            self.ratio.Divide(self.total_background)
            self._objects_to_keep.append(self.ratio)
        
        # Set up canvas
        self._setup_canvas()
    
    def _build_background_stack(self):
        """Build the background stack with grouped backgrounds"""
        self.background_stack = ROOT.THStack("bg_stack", "Background + Signal")
        self._objects_to_keep.append(self.background_stack)
        self.total_background = None
        
        # Define background groups
        groups = {
            "TTLJ": [],
            "TTLL": [],
            "ST": [],
            "DYJets": [],
            "QCD" : [],
            "Others": []
        }
        
        # Combine all histograms (backgrounds + signal) for grouping
        all_hists = {}
        if self.background_hists:
            all_hists.update(self.background_hists)
        if self.signal_hists:
            all_hists.update(self.signal_hists)
        
        if not all_hists:
            return
        
        # Group samples by category
        for name, hist in all_hists.items():
            if name.startswith("TTLJ"):
                groups["TTLJ"].append((name, hist))
            elif name.startswith("TTLL"):
                groups["TTLL"].append((name, hist))
            elif name.startswith("ST"):
                groups["ST"].append((name, hist))
            elif name.startswith("DYJets"):
                groups["DYJets"].append((name, hist))
            elif name.startswith("QCD"):
                groups["QCD"].append((name, hist))
            else:
                groups["Others"].append((name, hist))
        
        # Create combined histograms for each group
        self.grouped_hists = {}
        group_colors = {
            "TTLJ": SIGNAL_COLOR,
            "TTLL": BACKGROUND_COLORS[0],
            "ST": BACKGROUND_COLORS[1], 
            "DYJets": BACKGROUND_COLORS[2],
            "QCD": BACKGROUND_COLORS[3],
            "Others": BACKGROUND_COLORS[6]  # Gray
        }
        
        group_integrals = []
        
        for group_name, samples in groups.items():
            if not samples:
                continue
                
            # Create combined histogram for this group
            combined_hist = None
            total_integral = 0
            
            for sample_name, hist in samples:
                if not hist.GetSumw2N():
                    hist.Sumw2()
                    
                if combined_hist is None:
                    combined_hist = hist.Clone(f"combined_{group_name}")
                    combined_hist.SetDirectory(0)
                    self._objects_to_keep.append(combined_hist)
                else:
                    combined_hist.Add(hist)
                    
                total_integral += hist.Integral()
            
            if combined_hist is not None:
                # Style the combined histogram
                combined_hist.SetFillColor(group_colors[group_name])
                combined_hist.SetLineColor(group_colors[group_name])
                combined_hist.SetLineWidth(1)
                combined_hist.SetFillStyle(1001)
                combined_hist.SetStats(0)
                
                self.grouped_hists[group_name] = combined_hist
                group_integrals.append((group_name, total_integral))
                
                print(f"Group {group_name}: {total_integral:.1f} events from {len(samples)} samples")
        
        # Sort groups by integral (highest first for legend)
        group_integrals.sort(key=lambda x: x[1], reverse=True)
        
        # Add grouped histograms to stack (reverse order so highest is on top)
        for group_name, integral in reversed(group_integrals):
            hist = self.grouped_hists[group_name]
            self.background_stack.Add(hist)
            
            # Build total background
            if self.total_background is None:
                self.total_background = hist.Clone("total_background")
                self.total_background.SetDirectory(0)
                self.total_background.SetStats(0)
                self._objects_to_keep.append(self.total_background)
            else:
                self.total_background.Add(hist)
        
        # Store the grouped samples for legend ordering
        self.all_samples = group_integrals
    
    def _setup_canvas(self):
        """Setup the canvas with two pads for main plot and ratio"""
        ROOT.gStyle.SetOptStat(0)
        ROOT.gStyle.SetOptTitle(0)
        
        # Calculate y-axis range
        if "ymax" in self.config.keys() and self.config["ymax"] is not None:
            ymax = self.config["ymax"]
        else:
            data_max = self.data_hist.GetMaximum() if self.data_hist else 0
            stack_max = self.total_background.GetMaximum() if self.total_background else 0
            ymax = max(data_max, stack_max) * 1.5
            if ymax <= 0:
                ymax = 1000  # Default fallback
        
        ymin = self.config.get("ymin", 0.1 if self.config.get("logy", False) else 0.)
        
        # Apply log scale defaults
        if self.config.get('logy', False):
            if "ymin" not in self.config.keys():
                ymin = 1e-1
            if ymax < 1000:
                ymax = 1000
        
        # Debug information
        print(f"Y-axis range: [{ymin}, {ymax}]")
        if self.data_hist:
            print(f"Combined Data histogram max: {self.data_hist.GetMaximum()}")
        if self.total_background:
            print(f"Total stack histogram max: {self.total_background.GetMaximum()}")
        
        # Set CMS style
        CMS.SetEnergy(13.6)
        CMS.SetLumi("7.9104 fb^{-1}")
        CMS.SetExtraText("Preliminary")
        
        # Create canvas
        canvas_name = f"canvas_{id(self)}"
        self.canv = ROOT.TCanvas(canvas_name, canvas_name, 800, 800)
        ROOT.SetOwnership(self.canv, False)
        self._objects_to_keep.append(self.canv)
        
        # Create pads
        pad1_name = f"pad1_{id(self)}"
        pad2_name = f"pad2_{id(self)}"
        
        # Upper pad for main plot
        pad1 = ROOT.TPad(pad1_name, pad1_name, 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0.02)
        pad1.SetTopMargin(0.08)
        pad1.SetLeftMargin(0.12)
        pad1.SetRightMargin(0.05)
        ROOT.SetOwnership(pad1, False)
        self._objects_to_keep.append(pad1)
        
        # Lower pad for ratio
        pad2 = ROOT.TPad(pad2_name, pad2_name, 0, 0.0, 1, 0.3)
        pad2.SetTopMargin(0.02)
        pad2.SetBottomMargin(0.35)
        pad2.SetLeftMargin(0.12)
        pad2.SetRightMargin(0.05)
        ROOT.SetOwnership(pad2, False)
        self._objects_to_keep.append(pad2)
        
        # Apply log scale if needed
        if self.config.get('logy', False):
            pad1.SetLogy()
        
        pad1.Draw()
        pad2.Draw()
        
        self.pad1 = pad1
        self.pad2 = pad2
        
        # Store axis ranges
        self.xmin = self.config["xRange"][0]
        self.xmax = self.config["xRange"][-1]
        self.ymin = ymin
        self.ymax = ymax
        self.ratio_ymin = self.config["yRange"][0]
        self.ratio_ymax = self.config["yRange"][1]
        
        # Create legend (adjust size based on number of group entries)
        # Maximum 5 groups (TTLJ, TTLL, ST, DYJets, Others) + Data
        n_entries = 6  # Conservative estimate for grouped entries
        leg_height = 0.05 * (n_entries + 1)
        self.leg = ROOT.TLegend(0.8, 0.9 - leg_height, 0.92, 0.9)
        self.leg.SetTextSize(0.03)
        self.leg.SetBorderSize(0)
        self.leg.SetFillStyle(0)
        ROOT.SetOwnership(self.leg, False)
        self._objects_to_keep.append(self.leg)
    
    def draw(self):
        """Draw the stacked signal+background vs data comparison"""
        
        # =====================================
        # Upper pad - Main plot
        # =====================================
        self.pad1.cd()
        
        # Draw combined stack (backgrounds + signal)
        first_drawn = False
        if self.background_stack:
            self.background_stack.Draw("HIST")
            first_drawn = True
            
            # Set axis properties on the stack
            self.background_stack.SetMinimum(self.ymin)
            self.background_stack.SetMaximum(self.ymax)
            
            # Get the histogram from the stack for axis formatting
            stack_hist = self.background_stack.GetHistogram()
            if stack_hist:
                stack_hist.GetYaxis().SetTitle(self.config["yTitle"])
                stack_hist.GetYaxis().SetTitleSize(0.05)
                stack_hist.GetYaxis().SetTitleOffset(1.2)
                stack_hist.GetYaxis().SetLabelSize(0.04)
                stack_hist.GetXaxis().SetLabelSize(0)  # Hide x-axis labels on upper pad
                stack_hist.GetXaxis().SetTickLength(0.03)
                ROOT.SetOwnership(stack_hist, False)
                self._objects_to_keep.append(stack_hist)
        
        # Draw combined data histogram (as points with error bars)
        if self.data_hist:
            if not first_drawn:
                self.data_hist.SetMinimum(self.ymin)
                self.data_hist.SetMaximum(self.ymax)
                self.data_hist.GetYaxis().SetTitle(self.config["yTitle"])
                self.data_hist.GetYaxis().SetTitleSize(0.05)
                self.data_hist.GetYaxis().SetTitleOffset(1.2)
                self.data_hist.GetYaxis().SetLabelSize(0.04)
                self.data_hist.GetXaxis().SetLabelSize(0)
                self.data_hist.Draw("E1")
            else:
                self.data_hist.Draw("E1 SAME")
        
        # Add legend entries for grouped backgrounds (highest integral first)
        if hasattr(self, 'all_samples') and hasattr(self, 'grouped_hists'):
            for group_name, integral in self.all_samples:
                if group_name in self.grouped_hists:
                    # Use clean group names for legend
                    legend_name = group_name
                    if group_name == "ST":
                        legend_name = "Single Top"
                    elif group_name == "DYJets":
                        legend_name = "Drell-Yan"
                    self.leg.AddEntry(self.grouped_hists[group_name], legend_name, "F")
        
        if self.data_hist:
            self.leg.AddEntry(self.data_hist, "Data", "PE")
        
        self.leg.Draw()
        
        # Add CMS labels
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.045)
        latex.SetTextFont(62)
        latex.DrawLatex(0.12, 0.93, "CMS")
        
        latex.SetTextSize(0.035)
        latex.SetTextFont(52)
        latex.DrawLatex(0.20, 0.93, "Preliminary")
        
        latex.SetTextFont(42)
        latex.DrawLatex(0.73, 0.93, "7.9104 fb^{-1} (13.6 TeV)")
        
        ROOT.SetOwnership(latex, False)
        self._objects_to_keep.append(latex)
        
        # Redraw axis
        self.pad1.RedrawAxis()
        
        # =====================================
        # Lower pad - Ratio plot
        # =====================================
        self.pad2.cd()
        
        if self.ratio:
            # Create dummy histogram for ratio pad axis
            nb = self.total_background.GetNbinsX() if self.total_background else 100
            ratio_dummy = ROOT.TH1F(f"ratio_dummy_{id(self)}", "", nb, self.xmin, self.xmax)
            ratio_dummy.SetMinimum(self.ratio_ymin)
            ratio_dummy.SetMaximum(self.ratio_ymax)
            ratio_dummy.GetXaxis().SetTitle(self.config["xTitle"])
            ratio_dummy.GetYaxis().SetTitle("Data / MC")
            
            # Adjust font sizes for ratio pad
            ratio_dummy.GetXaxis().SetTitleSize(0.12)
            ratio_dummy.GetXaxis().SetTitleOffset(1.0)
            ratio_dummy.GetXaxis().SetLabelSize(0.10)
            ratio_dummy.GetYaxis().SetTitleSize(0.08)
            ratio_dummy.GetYaxis().SetTitleOffset(0.5)
            ratio_dummy.GetYaxis().SetLabelSize(0.08)
            ratio_dummy.GetYaxis().SetNdivisions(505)
            ratio_dummy.SetDirectory(0)
            
            ROOT.SetOwnership(ratio_dummy, False)
            self._objects_to_keep.append(ratio_dummy)
            ratio_dummy.Draw()
            
            # Create MC error band for ratio (total stack uncertainty)
            if self.total_background:
                mc_ratio_error = ROOT.TH1F(f"mc_ratio_error_{id(self)}", "", nb, self.xmin, self.xmax)
                for i in range(1, nb + 1):
                    mc_content = self.total_background.GetBinContent(i)
                    mc_error = self.total_background.GetBinError(i)
                    if mc_content > 0:
                        ratio_error = mc_error / mc_content  # Relative error
                        mc_ratio_error.SetBinContent(i, 1.0)
                        mc_ratio_error.SetBinError(i, ratio_error)
                    else:
                        mc_ratio_error.SetBinContent(i, 1.0)
                        mc_ratio_error.SetBinError(i, 0.0)
                
                mc_ratio_error.SetFillColor(ROOT.kBlack)
                mc_ratio_error.SetFillStyle(3013)  # Crosshatch pattern
                mc_ratio_error.SetMarkerSize(0)
                ROOT.SetOwnership(mc_ratio_error, False)
                self._objects_to_keep.append(mc_ratio_error)
                mc_ratio_error.Draw("E2 SAME")
            
            # Draw horizontal line at y=1
            line = ROOT.TLine(self.xmin, 1., self.xmax, 1.)
            line.SetLineStyle(2)
            line.SetLineColor(ROOT.kBlack)
            line.SetLineWidth(1)
            ROOT.SetOwnership(line, False)
            self._objects_to_keep.append(line)
            line.Draw()
            
            # Draw combined data ratio histogram
            self.ratio.SetLineColor(DATA_COLOR)
            self.ratio.SetLineWidth(2)
            self.ratio.SetMarkerStyle(20)
            self.ratio.SetMarkerSize(0.8)
            self.ratio.SetMarkerColor(DATA_COLOR)
            self.ratio.Draw("E1 SAME")
            
            # Redraw axis for ratio pad
            self.pad2.RedrawAxis()
        
        # Update canvas
        self.canv.Update()
    
    def save_as(self, filename):
        """Save the plot"""
        self.canv.SaveAs(filename)
    
    def close(self):
        """Clean up resources"""
        if hasattr(self, 'canv') and self.canv:
            self.canv.Close()
        self._objects_to_keep.clear()
        gc.collect()
    
    def __del__(self):
        """Destructor"""
        try:
            self.close()
        except:
            pass

def load_histogram(file_path, hist_name, systematic="Central", silent=False):
    """Load histogram from ROOT file"""
    root_file = ROOT.TFile.Open(file_path)
    if not root_file or root_file.IsZombie():
        if not silent:
            print(f"Error: Cannot open file {file_path}")
        return None
    
    # Navigate to systematic directory
    directory = root_file.Get(systematic)
    if not directory:
        if not silent:
            print(f"Error: Cannot find directory {systematic} in {file_path}")
        root_file.Close()
        return None
    
    hist = directory.Get(hist_name)
    if not hist:
        if not silent:
            print(f"Error: Cannot find histogram {hist_name} in {file_path}:{systematic}")
        root_file.Close()
        return None
    
    # Clone to avoid issues when file is closed
    hist_clone = hist.Clone(f"{os.path.basename(file_path)}_{hist_name}")
    hist_clone.SetDirectory(0)
    ROOT.SetOwnership(hist_clone, False)
    root_file.Close()
    
    return hist_clone

def combine_muon_data(data_dir, hist_name, systematic="Central"):
    """Load and combine Muon C, D and SingleMuon data histograms"""
    print("Loading data from Muon_C, Muon_D, and SingleMuon files")
    
    # Get all ROOT files in directory
    root_files = glob.glob(os.path.join(data_dir, "*.root"))
    
    data_hists = []
    total_events = 0
    
    for file_path in root_files:
        filename = os.path.basename(file_path)
        
        # Include only Muon_C, Muon_D, and SingleMuon files for data
        if (filename.startswith("Muon_C.root") or 
            filename.startswith("Muon_D.root") or
            filename.startswith("SingleMuon")):
            
            hist = load_histogram(file_path, hist_name, systematic)
            if hist:
                events = hist.Integral()
                data_hists.append((filename, hist, events))
                print(f"Loaded {filename}: {events:.1f} events")
                total_events += events
    
    if not data_hists:
        print("Warning: No data histograms found (Muon_C, Muon_D, or SingleMuon)")
        return None
    
    # Create combined histogram starting with the first one
    combined_hist = data_hists[0][1].Clone("Combined_Muon_Data")
    combined_hist.SetDirectory(0)
    
    # Add the rest
    for i in range(1, len(data_hists)):
        combined_hist.Add(data_hists[i][1])
    
    ROOT.SetOwnership(combined_hist, False)
    print(f"Combined Data Total (Muon_C + Muon_D + SingleMuon): {total_events:.1f} events")
    print(f"Total data files combined: {len(data_hists)}")
    
    return combined_hist

def load_signal_histograms(data_dir, hist_name, systematic="Central"):
    """Load signal histogram from TTLJ file"""
    signal_hists = {}
    ttlj_path = os.path.join(data_dir, "TTLJ_powheg.root")
    
    if os.path.exists(ttlj_path):
        sig_hist = load_histogram(ttlj_path, hist_name, systematic)
        if sig_hist:
            signal_hists["TTLJ_powheg"] = sig_hist
            print(f"Loaded TTLJ Signal: {sig_hist.Integral():.1f} events")
    
    return signal_hists

def load_background_histograms(data_dir, hist_name, systematic="Central"):
    """Load background histograms from all ROOT files except TTLJ, TBChannel, and data files"""
    background_hists = {}
    
    # Get all ROOT files in directory
    root_files = glob.glob(os.path.join(data_dir, "*.root"))
    
    for file_path in root_files:
        filename = os.path.basename(file_path)
        
        # Skip TTLJ (signal), TBChannel files, and all data files (Muon_, EGamma_, MuonEG_, SingleMuon)
        if (filename.startswith("TTLJ") or 
            filename.startswith("TBChannel") or
            filename.startswith("Muon_") or
            filename.startswith("EGamma_") or
            filename.startswith("MuonEG_") or
            filename.startswith("SingleMuon") or
            filename.startswith("ZZTo") or
            filename.startswith("WZTo") or
            filename.startswith("DYJets10") or
            filename.startswith("DYJets_MG") or
            filename.startswith("DYG") or
            filename.startswith("WJets") or
            filename.startswith("TTG") or
            filename.endswith("bcToE.root")):
            continue
        
        hist = load_histogram(file_path, hist_name, systematic, silent=True)
        if hist:
            sample_name = filename.replace(".root", "")
            background_hists[sample_name] = hist
            print(f"Loaded background {sample_name}: {hist.Integral():.1f} events")
    
    print(f"Total background samples loaded: {len(background_hists)}")
    return background_hists

def plot_signal_background_comparison(data_dir, hist_name, config, output_name="signal_background_comparison", systematic="Central"):
    """Create Signal+Background vs Data comparison plot"""
    
    print(f"Creating Signal+Background vs Data comparison for histogram: {hist_name}")
    print(f"Data directory: {data_dir}")
    
    # Load and combine muon data (C and D only, excluding E)
    combined_data_hist = combine_muon_data(data_dir, hist_name, systematic)
    
    # Load signal histograms (TTLJ)
    signal_hists = load_signal_histograms(data_dir, hist_name, systematic)
    
    # Load background histograms (all others except TTLJ, TBChannel, and data)
    background_hists = load_background_histograms(data_dir, hist_name, systematic)
    
    if not combined_data_hist and not signal_hists and not background_hists:
        print("Error: No histograms could be loaded!")
        return None
    
    if not combined_data_hist:
        print("Warning: No data histograms found!")
    if not signal_hists:
        print("Warning: No signal histograms found!")
    if not background_hists:
        print("Warning: No background histograms found!")
    
    # Apply rebinning if requested
    if "rebin" in config and config["rebin"] > 1:
        if combined_data_hist:
            combined_data_hist.Rebin(config["rebin"])
        for hist in signal_hists.values():
            hist.Rebin(config["rebin"])
        for hist in background_hists.values():
            hist.Rebin(config["rebin"])
    
    # Create and draw the plot
    canvas = None
    try:
        canvas = SignalBackgroundCanvas(combined_data_hist, signal_hists, background_hists, config)
        canvas.draw()
        
        # Save the plot
        canvas.save_as(f"{output_name}.png")
        canvas.save_as(f"{output_name}.pdf")
        
        print(f"Plot saved as {output_name}.png and {output_name}.pdf")
        
        return canvas
        
    except Exception as e:
        print(f"Error creating plot: {e}")
        if canvas:
            canvas.close()
        return None

def main():
    """Main function to create TTbar stacked MC vs Data plots"""
    
    print("TTbar Stacked MC vs Data Comparison Plotter")
    print("=" * 50)
    
    # Check if data directory exists
    if not os.path.exists(data_directory):
        print(f"Error: Data directory {data_directory} does not exist!")
        return
    
    # Configuration for the plot
    plot_config = {
        "xRange": [0, 2000],        # TTbar transverse mass range
        "yRange": [0.5, 2.0],       # Ratio plot range
        "xTitle": "m_{T}^{t#bar{t}} [GeV]",
        "yTitle": "Events / bin",
        "logy": True,               # Log scale
        "rebin": 5,                 # Rebin factor
        "ymax": 100000,             # Higher max for log scale
        "ymin": 0.1                 # Small positive value for log scale
    }
    
    # Create the plot
    hist_name = "TTbarTransverseMass_v2"
    output_name = "TTbar_StackedMC_comparison_2022"
    
    try:
        canvas = plot_signal_background_comparison(
            data_directory, 
            hist_name, 
            plot_config, 
            output_name
        )
        
        if canvas:
            print("Plot creation successful!")
            # Keep the canvas alive for a moment to ensure saving completes
            input("Press Enter to exit...")
            canvas.close()
        else:
            print("Plot creation failed!")
            
    except Exception as e:
        print(f"Error in main: {e}")

if __name__ == "__main__":
    main()