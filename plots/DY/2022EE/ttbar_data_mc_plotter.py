#!/usr/bin/env python3
# python ttbar_data_mc_plotter.py --option2 --bin-edges "15,20,25,30,35,40,45,50,55,60,64,68,72,76,81,86,91,96,101,106,110,115,120,126,133,141,150,160,171,185,200,220,243,273,320,380,440,510,600,700,830,1000,1500,3000"
data_directory = "/gv0/Users/achihwan/SKNanoRunlog/out/DY/2022EE"

import ROOT
import cmsstyle as CMS
import os
import glob
from array import array
import gc
import argparse

# Prevent ROOT from owning Python objects
ROOT.SetOwnership(ROOT.gROOT, False)

# Colors for plotting
DATA_COLOR = ROOT.kBlack
SIGNAL_COLOR = ROOT.TColor.GetColor("#e42536")  # Red for DYJets signal
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
        
        # Define background groups for DY analysis
        groups = {
            "DYJets": [],
            "TT": [],
            "VV": [],
            "ST": [],
            "QCD": [],
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
            if name.startswith("DYJets"):
                groups["DYJets"].append((name, hist))
            elif name.startswith("TTLJ"):
                groups["TT"].append((name, hist))
            elif name.startswith("TTLL"):
                groups["TT"].append((name, hist))
            elif name.startswith("ST"):
                groups["ST"].append((name, hist))
            elif name.startswith("QCD"):
                groups["QCD"].append((name, hist))
            elif name.startswith("WZ"):
                groups["VV"].append((name, hist))
            elif name.startswith("WW"):
                groups["VV"].append((name, hist))
            elif name.startswith("ZZ"):
                groups["VV"].append((name, hist))
            else:
                groups["Others"].append((name, hist))
        
        # Create combined histograms for each group
        self.grouped_hists = {}
        group_colors = {
            "DYJets": SIGNAL_COLOR,
            "TT": BACKGROUND_COLORS[0],
            "VV": BACKGROUND_COLORS[1],
            "ST": BACKGROUND_COLORS[2],
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
    
    def _set_custom_tick_labels(self, x_axis, hide_labels=True):
        """Set custom tick labels for option2 (20, 30, 40, 100, 200, 1000, 2000)"""
        # Hide default labels for upper pad - we'll draw custom ones
        if hide_labels:
            x_axis.SetLabelSize(0)
            x_axis.SetTickLength(0.02)  # Keep ticks but hide labels
        else:
            # For ratio pad, we'll use custom labels only
            x_axis.SetLabelSize(0)  # Hide default labels completely
            x_axis.SetTickLength(0.02)
        x_axis.SetMoreLogLabels(False)
        x_axis.SetNoExponent(True)
    
    def _draw_custom_x_labels(self, pad, is_ratio_pad=False):
        """Draw custom x-axis labels for option2"""
        if not self.config.get('option2', False):
            return
        
        # Only draw custom labels on ratio pad to avoid overlap
        if not is_ratio_pad:
            return
            
        custom_ticks = [20, 30, 40, 100, 200, 1000, 2000]
        valid_ticks = [tick for tick in custom_ticks if self.xmin <= tick <= self.xmax]
        
        if not valid_ticks:
            return
            
        pad.cd()
        
        # Calculate positions in NDC coordinates
        x_ndc_positions = []
        for tick in valid_ticks:
            if pad.GetLogx():
                if tick > 0 and self.xmin > 0 and self.xmax > 0:
                    import math
                    x_ndc = (math.log10(tick) - math.log10(self.xmin)) / (math.log10(self.xmax) - math.log10(self.xmin))
                else:
                    continue
            else:
                x_ndc = (tick - self.xmin) / (self.xmax - self.xmin)
            
            # Convert to pad coordinates
            x_pad = pad.GetLeftMargin() + x_ndc * (1 - pad.GetLeftMargin() - pad.GetRightMargin())
            x_ndc_positions.append((x_pad, str(int(tick))))
        
        # Draw custom labels only on ratio pad
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextAlign(23)  # Center, bottom
        latex.SetTextSize(0.08)
        
        # Position labels just below the ratio pad (using ratio pad coordinates)
        y_pos = 0.33  # Bottom area of ratio pad for tick labels
        
        for x_pos, label in x_ndc_positions:
            if pad.GetLeftMargin() < x_pos < (1 - pad.GetRightMargin()):  # Only draw if within pad boundaries
                latex.DrawLatex(x_pos, y_pos, label)
        
        # Draw x-axis title at the very bottom
        latex.SetTextAlign(22)  # Center, center
        latex.SetTextSize(0.10)
        x_center = 0.9
        y_title = 0.2  # Very bottom for axis title
        latex.DrawLatex(x_center, y_title, self.config["xTitle"])
        
        ROOT.SetOwnership(latex, False)
        self._objects_to_keep.append(latex)
        
        print(f"Option2: Drew custom labels at {[pos[1] for pos in x_ndc_positions]}")
    
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
        CMS.SetLumi("26.6717 fb^{-1}")
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
        pad1.SetLeftMargin(0.15)
        pad1.SetRightMargin(0.08)
        ROOT.SetOwnership(pad1, False)
        self._objects_to_keep.append(pad1)
        
        # Lower pad for ratio
        pad2 = ROOT.TPad(pad2_name, pad2_name, 0, 0.0, 1, 0.3)
        pad2.SetTopMargin(0.02)
        pad2.SetBottomMargin(0.35)
        pad2.SetLeftMargin(0.15)
        pad2.SetRightMargin(0.08)
        ROOT.SetOwnership(pad2, False)
        self._objects_to_keep.append(pad2)
        
        # Apply log scale if needed
        if self.config.get('logy', False):
            pad1.SetLogy()
        
        # Apply option2 settings (log scale x-axis)
        if self.config.get('option2', False):
            pad1.SetLogx()
            pad2.SetLogx()
        
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
        # Maximum 6 groups (DYJets, TTLJ, TTLL, ST, QCD, Others) + Data
        n_entries = 7  # Conservative estimate for grouped entries
        leg_height = 0.05 * (n_entries + 1)
        self.leg = ROOT.TLegend(0.65, 0.9 - leg_height, 0.92, 0.9)
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
                # Force x-axis range from config
                stack_hist.GetXaxis().SetRangeUser(self.xmin, self.xmax)
                
                # Apply option2 custom tick labels (hide default labels for upper pad)
                if self.config.get('option2', False):
                    self._set_custom_tick_labels(stack_hist.GetXaxis())
                
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
                # Force x-axis range from config
                self.data_hist.GetXaxis().SetRangeUser(self.xmin, self.xmax)
                
                # Apply option2 custom tick labels (hide default labels for upper pad)
                if self.config.get('option2', False):
                    self._set_custom_tick_labels(self.data_hist.GetXaxis())
                
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
        latex.DrawLatex(0.73, 0.93, "26.6717 fb^{-1} (13.6 TeV)")
        
        ROOT.SetOwnership(latex, False)
        self._objects_to_keep.append(latex)
        
        # Redraw axis
        self.pad1.RedrawAxis()
        
        # Draw custom x-axis labels for option2
        self._draw_custom_x_labels(self.pad1, is_ratio_pad=False)
        
        # =====================================
        # Lower pad - Ratio plot
        # =====================================
        self.pad2.cd()
        
        if self.ratio:
            # Create dummy histogram for ratio pad axis using same binning as actual histograms
            if self.total_background:
                ratio_dummy = self.total_background.Clone(f"ratio_dummy_{id(self)}")
                ratio_dummy.Reset()  # Clear content but keep binning
            elif self.data_hist:
                ratio_dummy = self.data_hist.Clone(f"ratio_dummy_{id(self)}")
                ratio_dummy.Reset()  # Clear content but keep binning
            else:
                ratio_dummy = ROOT.TH1F(f"ratio_dummy_{id(self)}", "", 100, self.xmin, self.xmax)
            ratio_dummy.SetMinimum(self.ratio_ymin)
            ratio_dummy.SetMaximum(self.ratio_ymax)
            # Set axis titles
            if self.config.get('option2', False):
                # For option2, we'll draw custom title and labels
                ratio_dummy.GetXaxis().SetTitle("")  # Hide default title
                ratio_dummy.GetXaxis().SetLabelSize(0)  # Hide all default labels
            else:
                ratio_dummy.GetXaxis().SetTitle(self.config["xTitle"])
                ratio_dummy.GetXaxis().SetLabelSize(0.10)
            
            ratio_dummy.GetYaxis().SetTitle("Data / MC")
            
            # Adjust font sizes for ratio pad
            ratio_dummy.GetXaxis().SetTitleSize(0.12)
            ratio_dummy.GetXaxis().SetTitleOffset(1.0)
            ratio_dummy.GetYaxis().SetTitleSize(0.08)
            ratio_dummy.GetYaxis().SetTitleOffset(0.5)
            ratio_dummy.GetYaxis().SetLabelSize(0.08)
            ratio_dummy.GetYaxis().SetNdivisions(505)
            ratio_dummy.SetDirectory(0)
            # Force x-axis range from config
            ratio_dummy.GetXaxis().SetRangeUser(self.xmin, self.xmax)
            
            # Apply option2 custom tick labels (hide all default labels for ratio pad too)
            if self.config.get('option2', False):
                self._set_custom_tick_labels(ratio_dummy.GetXaxis(), hide_labels=True)
            
            ROOT.SetOwnership(ratio_dummy, False)
            self._objects_to_keep.append(ratio_dummy)
            ratio_dummy.Draw()
            
            # Create MC error band for ratio (total stack uncertainty)
            if self.total_background:
                mc_ratio_error = self.total_background.Clone(f"mc_ratio_error_{id(self)}")
                mc_ratio_error.Reset()  # Clear content but keep binning
                
                for i in range(1, mc_ratio_error.GetNbinsX() + 1):
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
                mc_ratio_error.SetDirectory(0)
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
            self.ratio.GetXaxis().SetRangeUser(self.xmin, self.xmax)
            self.ratio.SetLineColor(DATA_COLOR)
            self.ratio.SetLineWidth(2)
            self.ratio.SetMarkerStyle(20)
            self.ratio.SetMarkerSize(0.8)
            self.ratio.SetMarkerColor(DATA_COLOR)
            self.ratio.Draw("E1 SAME")
            
            # Redraw axis for ratio pad
            self.pad2.RedrawAxis()
            
            # Draw custom x-axis labels for option2
            self._draw_custom_x_labels(self.pad2, is_ratio_pad=True)
        
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

def rebin_histogram_custom(hist, bin_edges):
    """Rebin histogram with custom bin edges"""
    if not hist:
        return None
    
    # Convert to array for ROOT
    bin_edges_array = array('d', bin_edges)
    n_bins = len(bin_edges) - 1
    
    # Create new histogram with custom binning
    new_hist = ROOT.TH1F(f"{hist.GetName()}_rebinned", hist.GetTitle(), n_bins, bin_edges_array)
    new_hist.SetDirectory(0)
    
    # Fill new histogram by integrating over bin ranges
    for i in range(1, n_bins + 1):
        bin_low = bin_edges[i-1]
        bin_high = bin_edges[i]
        
        # Find corresponding bins in original histogram
        bin_low_orig = hist.FindBin(bin_low)
        bin_high_orig = hist.FindBin(bin_high)
        
        # Handle bin boundaries carefully
        content = 0.0
        error_sq = 0.0
        
        for j in range(bin_low_orig, bin_high_orig + 1):
            if j < 1 or j > hist.GetNbinsX():
                continue
                
            bin_center = hist.GetBinCenter(j)
            bin_width = hist.GetBinWidth(j)
            bin_low_edge = bin_center - bin_width/2
            bin_high_edge = bin_center + bin_width/2
            
            # Calculate overlap fraction
            overlap_low = max(bin_low, bin_low_edge)
            overlap_high = min(bin_high, bin_high_edge)
            
            if overlap_high > overlap_low:
                overlap_fraction = (overlap_high - overlap_low) / bin_width
                bin_content = hist.GetBinContent(j)
                bin_error = hist.GetBinError(j)
                
                content += bin_content * overlap_fraction
                error_sq += (bin_error * overlap_fraction) ** 2
        
        new_hist.SetBinContent(i, content)
        new_hist.SetBinError(i, error_sq ** 0.5)
    
    # Copy histogram properties
    new_hist.SetLineColor(hist.GetLineColor())
    new_hist.SetLineWidth(hist.GetLineWidth())
    new_hist.SetFillColor(hist.GetFillColor())
    new_hist.SetFillStyle(hist.GetFillStyle())
    new_hist.SetMarkerColor(hist.GetMarkerColor())
    new_hist.SetMarkerStyle(hist.GetMarkerStyle())
    new_hist.SetMarkerSize(hist.GetMarkerSize())
    
    ROOT.SetOwnership(new_hist, False)
    return new_hist

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
    print("Loading data from Muon_E, Muon_F, and SingleMuon files")
    
    # Get all ROOT files in directory
    root_files = glob.glob(os.path.join(data_dir, "*.root"))
    
    data_hists = []
    total_events = 0
    
    for file_path in root_files:
        filename = os.path.basename(file_path)
        
        # Include only Muon_E, Muon_F, and SingleMuon files for data
        if (filename.startswith("Muon_E")or
            filename.startswith("Muon_F.root") or 
            filename.startswith("Muon_G.root") ):
            
            hist = load_histogram(file_path, hist_name, systematic)
            if hist:
                events = hist.Integral()
                data_hists.append((filename, hist, events))
                print(f"Loaded {filename}: {events:.1f} events")
                total_events += events
    
    if not data_hists:
        print("Warning: No data histograms found (Muon_E, Muon_F, , Muon_G)")
        return None
    
    # Create combined histogram starting with the first one
    combined_hist = data_hists[0][1].Clone("Combined_Muon_Fata")
    combined_hist.SetDirectory(0)
    
    # Add the rest
    for i in range(1, len(data_hists)):
        combined_hist.Add(data_hists[i][1])
    
    ROOT.SetOwnership(combined_hist, False)
    print(f"Combined Data Total (Muon_E + Muon_F + SingleMuon): {total_events:.1f} events")
    print(f"Total data files combined: {len(data_hists)}")
    
    return combined_hist

def load_signal_histograms(data_dir, hist_name, systematic="Central"):
    """Load signal histogram from DYJets file"""
    signal_hists = {}
    
    # Get all ROOT files in directory
    root_files = glob.glob(os.path.join(data_dir, "*.root"))
    
    for file_path in root_files:
        filename = os.path.basename(file_path)
        
        # Load DYJets files as signal
        if filename.startswith("DYJets.root") or filename.startswith("DYJets10to50.root") :
            sig_hist = load_histogram(file_path, hist_name, systematic)
            if sig_hist:
                sample_name = filename.replace(".root", "")
                signal_hists[sample_name] = sig_hist
                print(f"Loaded DYJets Signal: {sig_hist.Integral():.1f} events")
    
    return signal_hists

def load_background_histograms(data_dir, hist_name, systematic="Central"):
    """Load background histograms from all ROOT files except DYJets, TBChannel, and data files"""
    background_hists = {}
    
    # Get all ROOT files in directory
    root_files = glob.glob(os.path.join(data_dir, "*.root"))
    
    for file_path in root_files:
        filename = os.path.basename(file_path)
        
        # Skip DYJets (signal), TBChannel files, and all data files (Muon_, EGamma_, MuonEG_, SingleMuon)
        if (filename.startswith("DYJet") or
            filename.startswith("TBChannel") or
            filename.startswith("Muon_") or
            filename.startswith("EGamma_") or
            filename.startswith("MuonEG_") or
            filename.startswith("SingleMuon") or
            filename.startswith("ZZTwo") or
            filename.startswith("WZTo")or
            filename.startswith("WJets_MG")or
            filename.startswith("DYG") ):
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
    
    # Load signal histograms (DYJets)
    signal_hists = load_signal_histograms(data_dir, hist_name, systematic)
    
    # Load background histograms (all others except DYJets, TBChannel, and data)
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
    if "custom_bin_edges" in config and config["custom_bin_edges"]:
        # Apply custom binning
        bin_edges = config["custom_bin_edges"]
        print(f"Applying custom binning with {len(bin_edges)-1} bins")
        
        if combined_data_hist:
            combined_data_hist = rebin_histogram_custom(combined_data_hist, bin_edges)
        
        # Rebin signal histograms
        rebinned_signal_hists = {}
        for name, hist in signal_hists.items():
            rebinned_hist = rebin_histogram_custom(hist, bin_edges)
            if rebinned_hist:
                rebinned_signal_hists[name] = rebinned_hist
        signal_hists = rebinned_signal_hists
        
        # Rebin background histograms
        rebinned_background_hists = {}
        for name, hist in background_hists.items():
            rebinned_hist = rebin_histogram_custom(hist, bin_edges)
            if rebinned_hist:
                rebinned_background_hists[name] = rebinned_hist
        background_hists = rebinned_background_hists
        
        # Update x-axis range to match custom bins
        config["xRange"] = [bin_edges[0], bin_edges[-1]]
        
    elif "rebin" in config and config["rebin"] > 1:
        # Apply uniform rebinning
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

def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="DY Dilepton Z Mass Stacked MC vs Data Comparison Plotter")
    
    # Binning options
    parser.add_argument("--xmin", type=float, default=15, help="X-axis minimum value (default: 15)")
    parser.add_argument("--xmax", type=float, default=3000, help="X-axis maximum value (default: 3000)")
    parser.add_argument("--rebin", type=int, default=2, help="Rebin factor (default: 2)")
    parser.add_argument("--bin-edges", type=str, help="Custom bin edges as comma-separated values (e.g., '15,20,60,64,106,110,3000'). Overrides xmin/xmax/rebin.")
    
    # Y-axis options
    parser.add_argument("--ymin", type=float, default=0.1, help="Y-axis minimum value (default: 0.1)")
    parser.add_argument("--ymax", type=float, default=1e6, help="Y-axis maximum value (default: 1e6)")
    parser.add_argument("--logy", action="store_true", default=True, help="Use logarithmic Y-axis (default: True)")
    parser.add_argument("--linear", action="store_true", help="Use linear Y-axis (overrides --logy)")
    
    # Ratio plot options
    parser.add_argument("--ratio-ymin", type=float, default=0.4, help="Ratio plot Y-axis minimum (default: 0.4)")
    parser.add_argument("--ratio-ymax", type=float, default=1.6, help="Ratio plot Y-axis maximum (default: 1.6)")
    
    # Output options
    parser.add_argument("--output", type=str, default="DY_DileptonZMass_comparison_2022", help="Output filename prefix")
    parser.add_argument("--hist-name", type=str, default="DileptonMass", help="Histogram name to plot")
    
    # Special options
    parser.add_argument("--option2", action="store_true", help="Apply log scale x-axis with custom tick labels (20,30,40,100,200,1000,2000)")
    
    return parser.parse_args()

def main():
    """Main function to create DY dilepton Z mass stacked MC vs Data plots"""
    
    # Parse command line arguments
    args = parse_args()
    
    print("DY Dilepton Z Mass Stacked MC vs Data Comparison Plotter")
    print("=" * 56)
    
    # Parse custom bin edges if provided
    custom_bin_edges = None
    if args.bin_edges:
        try:
            custom_bin_edges = [float(x.strip()) for x in args.bin_edges.split(',')]
            custom_bin_edges.sort()  # Ensure ascending order
            print(f"Custom bin edges: {custom_bin_edges}")
        except ValueError:
            print(f"Error: Invalid bin edges format '{args.bin_edges}'. Use comma-separated numbers.")
            return
    else:
        print(f"X-axis range: [{args.xmin}, {args.xmax}]")
        print(f"Rebin factor: {args.rebin}")
    
    if args.option2:
        print("Option2: Log scale X-axis with custom tick labels enabled")
    
    print(f"Y-axis: {'Log' if args.logy and not args.linear else 'Linear'}")
    print(f"Output: {args.output}")
    print("-" * 56)
    
    # Check if data directory exists
    if not os.path.exists(data_directory):
        print(f"Error: Data directory {data_directory} does not exist!")
        return
    
    # Configuration for the plot using command line arguments
    plot_config = {
        "xRange": [args.xmin, args.xmax],
        "yRange": [args.ratio_ymin, args.ratio_ymax],
        "xTitle": "m_{ll} [GeV]",
        "yTitle": "Events / bin",
        "logy": args.logy and not args.linear,
        "rebin": args.rebin,
        "ymax": args.ymax,
        "ymin": args.ymin,
        "custom_bin_edges": custom_bin_edges,
        "option2": args.option2
    }
    
    # Create the plot
    hist_name = args.hist_name
    output_name = args.output
    
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