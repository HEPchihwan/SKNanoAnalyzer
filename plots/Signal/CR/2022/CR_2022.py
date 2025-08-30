#!/usr/bin/env python3

data_directory = "/gv0/Users/achihwan/SKNanoRunlog/out/LRSM_TBChannel/2022"

import ROOT
import cmsstyle as CMS
import os
import glob
from array import array
import gc
import argparse
import numpy as np

# Prevent ROOT from owning Python objects
ROOT.SetOwnership(ROOT.gROOT, False)

# Colors for plotting
DATA_COLOR = ROOT.kBlack
SIGNAL_COLOR = ROOT.TColor.GetColor("#00C853")  # Green for TTLJ signal
BACKGROUND_COLORS = [
    ROOT.TColor.GetColor("#5790fc"),  # Blue
    ROOT.TColor.GetColor("#f89c20"),  # Orange  
    ROOT.TColor.GetColor("#964a8b"),  # Purple
    ROOT.TColor.GetColor("#CDDC39"),  # Lime
    ROOT.TColor.GetColor("#009688"),  # Teal
    ROOT.TColor.GetColor("#795548"),  # Brown
    ROOT.TColor.GetColor("#00bfae"),  # Cyan
    ROOT.TColor.GetColor("#9c9ca1")   # Gray for "Others"
]

class SignalBackgroundCanvas():
    def __init__(self, combined_data_hist, signal_hists, background_hists, config, draw_tb_lines=True, target_tb_sample=None):
        super().__init__()
        
        self.data_hist = combined_data_hist
        self.signal_hists = signal_hists
        self.background_hists = background_hists
        self.config = config
        self.draw_tb_lines = draw_tb_lines
        self.target_tb_sample = target_tb_sample  # Specific TB sample to draw as line
        
        # Keep references to prevent garbage collection
        self._objects_to_keep = []
        
        # Separate TB samples for line drawing (only if enabled)
        self.tb_hists = {}
        if self.draw_tb_lines:
            self._separate_tb_samples()
        
        # Build background stack (excludes TB samples if line drawing is enabled)
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
        
        # Style TB histograms for line drawing
        self._style_tb_histograms()
        
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
            "TTZ,TTW" : [],
            "TTH" : [],
            "TTTT" : [],
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
            if name.startswith("TB") and self.draw_tb_lines:
                # Skip ALL TB samples if TB lines are enabled - they will be drawn as lines only
                continue
            elif name.startswith("TTLJ"):
                groups["TTLJ"].append((name, hist))
            elif name.startswith("TTLL"):
                groups["TTLL"].append((name, hist))
            elif name.startswith("ST"):
                groups["ST"].append((name, hist))
            elif name.startswith("DYJets"):
                groups["DYJets"].append((name, hist))
            elif name.startswith("TTZ") or name.startswith("TTW"):
                groups["TTZ,TTW"].append((name, hist))
            elif name.startswith("TTH"):
                groups["TTH"].append((name, hist))
            elif name.startswith("TTTT"):
                groups["TTTT"].append((name, hist))
            elif name.startswith("QCD"):
                groups["QCD"].append((name, hist))
            elif name.startswith("TB"):
                # Include TB samples in stack only if TB lines are disabled
                groups["Others"].append((name, hist))
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
            "TTZ,TTW": BACKGROUND_COLORS[4],
            "TTH": BACKGROUND_COLORS[5],
            "TTTT": BACKGROUND_COLORS[6],
            "Others": BACKGROUND_COLORS[7]  
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
    
    def _separate_tb_samples(self):
        """Separate TB samples from background histograms for line drawing"""
        # TB samples are now loaded separately and passed via signal_hists
        # Extract TB samples from signal_hists for line drawing
        if self.signal_hists:
            for name, hist in self.signal_hists.items():
                if name.startswith("TB"):
                    self.tb_hists[name] = hist
                    print(f"Using TB sample for line drawing: {name}")
        
        if not self.tb_hists:
            print("No TB samples found for line drawing")
    
    def _style_tb_histograms(self):
        """Style TB histograms as red lines"""
        for name, hist in self.tb_hists.items():
            if not hist.GetSumw2N():
                hist.Sumw2()
            hist.SetLineColor(ROOT.kRed)
            hist.SetLineWidth(2)
            hist.SetMarkerSize(0)
            hist.SetFillStyle(0)  # No fill
            hist.SetStats(0)
            self._objects_to_keep.append(hist)
            print(f"Styled TB histogram {name}: {hist.Integral():.1f} events")
    
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
                ymax = 1e8  # Default fallback
        
        ymin = self.config.get("ymin", 0.1 if self.config.get("logy", False) else 0.)
        
        # Apply log scale defaults
        if self.config.get('logy', False):
            if "ymin" not in self.config.keys():
                ymin = 1e-1
            if ymax < 1000:
                ymax = 1e8
        
        # Debug information
        print(f"Y-axis range: [{ymin}, {ymax}]")
        if self.data_hist:
            print(f"Combined Data histogram max: {self.data_hist.GetMaximum()}")
        if self.total_background:
            print(f"Total stack histogram max: {self.total_background.GetMaximum()}")
        
        # Set CMS style
        CMS.SetEnergy(13.6)
        CMS.SetLumi("7.9804 fb^{-1}")
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
        
        # Store common bin range for both pads - calculate once here
        if self.total_background:
            self.bin_min = self.total_background.FindBin(self.xmin)
            self.bin_max = self.total_background.FindBin(self.xmax)
        elif self.data_hist:
            self.bin_min = self.data_hist.FindBin(self.xmin)
            self.bin_max = self.data_hist.FindBin(self.xmax)
        else:
            # Fallback - create a dummy histogram to get bin numbers
            dummy = ROOT.TH1F("dummy_for_bins", "", 200, 0, 2000)
            self.bin_min = dummy.FindBin(self.xmin)
            self.bin_max = dummy.FindBin(self.xmax)
            dummy.Delete()
        
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
                # Use the pre-calculated common bin range
                stack_hist.GetXaxis().SetRange(self.bin_min, self.bin_max)
                stack_hist.GetXaxis().SetRangeUser(self.xmin, self.xmax)
                stack_hist.GetYaxis().SetTitle(self.config["yTitle"])
                stack_hist.GetYaxis().SetTitleSize(0.05)
                stack_hist.GetYaxis().SetTitleOffset(1.2)
                stack_hist.GetYaxis().SetLabelSize(0.04)
                stack_hist.GetXaxis().SetLabelSize(0)  # Hide x-axis labels on upper pad
                stack_hist.GetXaxis().SetTickLength(0.03)
                ROOT.SetOwnership(stack_hist, False)
                self._objects_to_keep.append(stack_hist)
        
        # Draw MC error band for main plot
        if self.total_background:
            mc_error_band = self.total_background.Clone(f"mc_error_band_{id(self)}")
            mc_error_band.SetFillColor(ROOT.kBlack)
            mc_error_band.SetFillStyle(3013)  # Crosshatch pattern
            mc_error_band.SetMarkerSize(0)
            mc_error_band.GetXaxis().SetRange(self.bin_min, self.bin_max)
            mc_error_band.GetXaxis().SetRangeUser(self.xmin, self.xmax)
            ROOT.SetOwnership(mc_error_band, False)
            self._objects_to_keep.append(mc_error_band)
            mc_error_band.Draw("E2 SAME")
        
        # Draw combined data histogram (as points with error bars)
        if self.data_hist:
            # Use the same pre-calculated common bin range
            self.data_hist.GetXaxis().SetRange(self.bin_min, self.bin_max)
            self.data_hist.GetXaxis().SetRangeUser(self.xmin, self.xmax)
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
        
        # Draw TB histograms as red lines on pad1 only (if enabled)
        if self.draw_tb_lines:
            for name, hist in self.tb_hists.items():
                # Apply same bin range as other histograms
                hist.GetXaxis().SetRange(self.bin_min, self.bin_max)
                hist.GetXaxis().SetRangeUser(self.xmin, self.xmax)
                hist.Draw("HIST SAME")
        
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
        
        # Add TB samples to legend (if drawn as lines)
        if self.draw_tb_lines:
            for name, hist in self.tb_hists.items():
                # Extract legend name starting from WR for TB samples
                if name.startswith("TB") and "WR" in name:
                    # Find the position of WR and take substring from there
                    wr_pos = name.find("WR")
                    legend_name = name[wr_pos:] if wr_pos != -1 else name
                else:
                    legend_name = name
                self.leg.AddEntry(hist, legend_name, "L")
        
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
        latex.DrawLatex(0.73, 0.93, "7.9804 fb^{-1} (13.6 TeV)")
        
        ROOT.SetOwnership(latex, False)
        self._objects_to_keep.append(latex)
        
        # Redraw axis
        self.pad1.RedrawAxis()
        
        # =====================================
        # Lower pad - Ratio plot
        # =====================================
        self.pad2.cd()
        
        if self.ratio:
            # Create dummy histogram for ratio pad axis using exactly the same binning as pad1
            if self.background_stack and self.background_stack.GetHistogram():
                ratio_dummy = self.background_stack.GetHistogram().Clone(f"ratio_dummy_{id(self)}")
                ratio_dummy.Reset()  # Clear content but keep binning
            elif self.total_background:
                ratio_dummy = self.total_background.Clone(f"ratio_dummy_{id(self)}")
                ratio_dummy.Reset()  # Clear content but keep binning
            else:
                ratio_dummy = ROOT.TH1F(f"ratio_dummy_{id(self)}", "", 100, self.xmin, self.xmax)
            
            # Use exactly the same pre-calculated bin range as pad1
            ratio_dummy.GetXaxis().SetRange(self.bin_min, self.bin_max)
            ratio_dummy.GetXaxis().SetRangeUser(self.xmin, self.xmax)
            
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
            ratio_dummy.GetXaxis().SetTickLength(0.03)  # Match pad1 tick length
            ratio_dummy.SetDirectory(0)
            
            ROOT.SetOwnership(ratio_dummy, False)
            self._objects_to_keep.append(ratio_dummy)
            ratio_dummy.Draw()
            
            # Create MC error band for ratio (total stack uncertainty)
            if self.total_background:
                mc_ratio_error = self.total_background.Clone(f"mc_ratio_error_{id(self)}")
                mc_ratio_error.Reset()  # Clear content but keep binning
                # Apply same pre-calculated bin range as pad1
                mc_ratio_error.GetXaxis().SetRange(self.bin_min, self.bin_max)
                mc_ratio_error.GetXaxis().SetRangeUser(self.xmin, self.xmax)
                
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
            
            # Draw combined data ratio histogram with same pre-calculated bin range
            self.ratio.GetXaxis().SetRange(self.bin_min, self.bin_max)
            self.ratio.GetXaxis().SetRangeUser(self.xmin, self.xmax)
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

def rebin_to_common_binning(hist, xmin=0, xmax=2000, nbins=200, custom_bins=None):
    """Rebin histogram to common binning to avoid merge issues"""
    if not hist:
        return None
    
    # Create new histogram with common binning
    if custom_bins is not None:
        # Use custom bin edges, but handle overflow bins intelligently
        bin_edges = np.array(custom_bins, dtype=float)
        
        # If the last bin edge is much larger than xmax, treat it as overflow
        # and truncate at xmax for display purposes
        if len(bin_edges) >= 2 and bin_edges[-1] > xmax * 2:
            # Find bins up to xmax and add xmax as final edge
            valid_bins = bin_edges[bin_edges <= xmax]
            if len(valid_bins) == 0 or valid_bins[-1] < xmax:
                bin_edges = np.append(valid_bins, xmax)
            else:
                bin_edges = valid_bins
        
        new_hist = ROOT.TH1F(f"{hist.GetName()}_rebinned", hist.GetTitle(), len(bin_edges)-1, bin_edges)
        nbins = len(bin_edges) - 1
    else:
        # Use uniform binning
        new_hist = ROOT.TH1F(f"{hist.GetName()}_rebinned", hist.GetTitle(), nbins, xmin, xmax)
    
    new_hist.SetDirectory(0)
    
    # Fill new histogram by integrating over bin ranges
    for i in range(1, nbins + 1):
        bin_low = new_hist.GetBinLowEdge(i)
        bin_high = new_hist.GetBinLowEdge(i + 1)
        
        # Get content and error from original histogram in this range
        content = hist.Integral(hist.FindBin(bin_low), hist.FindBin(bin_high))
        
        # Calculate error properly
        error_sq = 0
        for j in range(hist.FindBin(bin_low), hist.FindBin(bin_high) + 1):
            if j >= 1 and j <= hist.GetNbinsX():
                error_sq += hist.GetBinError(j) ** 2
        
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
    
    # Apply common binning to avoid merge issues - but skip for discrete variables
    if ("jetnum" in hist_name.lower() or "num" in hist_name.lower() or 
        hist_clone.GetNbinsX() <= 20):  # Assume discrete if <= 20 bins
        return hist_clone
    else:
        # Check if custom_bins is available in the global scope (passed from main)
        custom_bins = getattr(load_histogram, 'custom_bins', None)
        # Get xmax from current xRange if available
        xmax_for_rebin = getattr(load_histogram, 'xmax_for_rebin', 2000)
        rebinned_hist = rebin_to_common_binning(hist_clone, xmin=0, xmax=xmax_for_rebin, nbins=200, custom_bins=custom_bins)
        return rebinned_hist if rebinned_hist else hist_clone

def combine_Muon_Fata(data_dir, hist_name, systematic="Central"):
    """Load and combine Muon C, D and SingleMuon data histograms"""
    print("Loading data from Muon_E, Muon_F, and SingleMuon files")
    
    # Get all ROOT files in directory
    root_files = glob.glob(os.path.join(data_dir, "*.root"))
    
    data_hists = []
    total_events = 0
    
    for file_path in root_files:
        filename = os.path.basename(file_path)
        
        # Include only Muon_E, Muon_F, and SingleMuon files for data
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
        print("Warning: No data histograms found (Muon_E, Muon_F, or Muon_G)")
        return None
    
    # Create combined histogram starting with the first one
    combined_hist = data_hists[0][1].Clone("Combined_Muon_Fata")
    combined_hist.SetDirectory(0)
    
    # Add the rest
    for i in range(1, len(data_hists)):
        combined_hist.Add(data_hists[i][1])
    
    ROOT.SetOwnership(combined_hist, False)
    print(f"Combined Data Total (Muon_E + Muon_F + Muon_G): {total_events:.1f} events")
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

def load_background_histograms(data_dir, hist_name, systematic="Central", draw_tb_lines=True):
    """Load background histograms from all ROOT files except data files"""
    background_hists = {}
    
    # Get all ROOT files in directory
    root_files = glob.glob(os.path.join(data_dir, "*.root"))
    
    for file_path in root_files:
        filename = os.path.basename(file_path)
        sample_name = filename.replace(".root", "")
        
        # Skip data files (Muon_E, Muon_F, Muon_G) - everything else becomes background
        if ((filename.startswith("Muon_C.root") or 
            filename.startswith("Muon_D.root") or
            filename.startswith("SingleMuon") or
            filename.startswith("Muon_G.root"))or
            filename.startswith("Muon_") or
            filename.startswith("EGamma_") or
            filename.startswith("MuonEG_") or
            filename.startswith("SingleMuon") or
            filename.startswith("ZZTo") or
            filename.startswith("WZTo") or
            filename.startswith("DYJets_MG") or
            filename.startswith("DYG") or
            filename.startswith("TTG") or
            filename.startswith("DYJets_MG") or
            filename.startswith("DYJets10to50_MG") or
            filename.endswith("bcToE.root") ):
            continue
        
        # Skip TB samples if TB lines are disabled (no-tb-lines option)
        if sample_name.startswith("TB") and not draw_tb_lines:
            print(f"Skipping TB sample (TB lines disabled): {sample_name}")
            continue
        
        hist = load_histogram(file_path, hist_name, systematic, silent=True)
        if hist:
            background_hists[sample_name] = hist
            print(f"Loaded background {sample_name}: {hist.Integral():.1f} events")
    
    print(f"Total background samples loaded: {len(background_hists)}")
    return background_hists

def get_tb_samples(data_dir):
    """Get list of TB sample names from the data directory"""
    tb_samples = []
    root_files = glob.glob(os.path.join(data_dir, "*.root"))
    
    for file_path in root_files:
        filename = os.path.basename(file_path)
        sample_name = filename.replace(".root", "")
        if sample_name.startswith("TB"):
            tb_samples.append(sample_name)
    
    return sorted(tb_samples)

def format_tb_legend_name(tb_sample_name):
    """Format TB sample name for legend display (show from WR onwards)"""
    if tb_sample_name.startswith("TB") and "WR" in tb_sample_name:
        wr_pos = tb_sample_name.find("WR")
        return tb_sample_name[wr_pos:] if wr_pos != -1 else tb_sample_name
    return tb_sample_name

def load_tb_sample_directly(data_dir, tb_sample_name, hist_name, systematic="Central"):
    """Load a specific TB sample histogram directly"""
    file_path = os.path.join(data_dir, f"{tb_sample_name}.root")
    if os.path.exists(file_path):
        hist = load_histogram(file_path, hist_name, systematic, silent=True)
        if hist:
            print(f"Loaded TB sample directly: {tb_sample_name}: {hist.Integral():.1f} events")
            return hist
    print(f"TB sample file not found: {file_path}")
    return None

def plot_signal_background_comparison(data_dir, hist_name, config, output_name="signal_background_comparison", systematic="Central", draw_tb_lines=True, target_tb_sample=None):
    """Create Signal+Background vs Data comparison plot"""
    
    print(f"Creating Signal+Background vs Data comparison for histogram: {hist_name}")
    print(f"Data directory: {data_dir}")
    print(f"TB lines drawing: {'Enabled' if draw_tb_lines else 'Disabled'}")
    if target_tb_sample:
        print(f"Target TB sample: {target_tb_sample}")
    
    # Load and combine muon data (Muon_E, Muon_F, Muon_G)
    combined_data_hist = combine_Muon_Fata(data_dir, hist_name, systematic)
    
    # Load signal histograms (empty - all MC becomes background)
    signal_hists = {}
    
    # Load background histograms (all MC files, excluding TB samples)
    background_hists = load_background_histograms(data_dir, hist_name, systematic, draw_tb_lines=draw_tb_lines)
    
    # Load TB sample separately if TB lines are enabled and target is specified
    if draw_tb_lines and target_tb_sample:
        tb_hist = load_tb_sample_directly(data_dir, target_tb_sample, hist_name, systematic)
        if tb_hist:
            # Add to signal_hists so it gets passed to the canvas but won't be in background stack
            signal_hists[target_tb_sample] = tb_hist
    
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
            print(f"Rebinned data histogram: {combined_data_hist.Integral()} events")
        for hist in signal_hists.values():
            hist.Rebin(config["rebin"])
        for hist in background_hists.values():
            hist.Rebin(config["rebin"])
    
    # Create and draw the plot
    canvas = None
    try:
        # Store additional info for TB sample loading
        canvas = SignalBackgroundCanvas(combined_data_hist, signal_hists, background_hists, config, draw_tb_lines=draw_tb_lines, target_tb_sample=target_tb_sample)
        canvas.data_dir = data_dir
        canvas.hist_name = hist_name
        canvas.systematic = systematic
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

def parse_arguments():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(description="TTbar Stacked MC vs Data Comparison Plotter")
    
    # Basic options
    parser.add_argument("--data-dir", type=str, default=data_directory,
                      help=f"Data directory path (default: {data_directory})")
    parser.add_argument("--hist-name", type=str, default="Topjetnum",
                      help="Histogram name to plot (default: Topjetnum)")
    parser.add_argument("--output-name", type=str, default="Topjetnum_CR_2022",
                      help="Output file name prefix (default: Topjetnum_CR_2022)")
    
    # TB sample options
    parser.add_argument("--no-tb-lines", action="store_true",
                      help="Disable TB sample line drawing (include TB in stack instead)")
    
    # Custom binning options
    parser.add_argument("--custom-bins", type=str, default=None,
                      help="Custom bin edges as comma-separated values (e.g., '0,120,300')")
    
    # Plot configuration
    parser.add_argument("--x-range", type=str, default="0,5",
                      help="X-axis range as 'min,max' (default: '0,5')")
    parser.add_argument("--y-range", type=str, default="0.5,1.5",
                      help="Ratio plot Y-axis range as 'min,max' (default: '0.5,1.5')")
    parser.add_argument("--x-title", type=str, default="m_{ll} [GeV]",
                      help="X-axis title (default: 'Topjetnum [GeV]')")
    parser.add_argument("--y-title", type=str, default="Events / bin",
                      help="Y-axis title (default: 'Events / bin')")
    parser.add_argument("--logy", action="store_true",
                      help="Use logarithmic Y-axis")
    parser.add_argument("--rebin", type=int, default=1,
                      help="Rebin factor (default: 1)")
    parser.add_argument("--ymax", type=float, default=1e6,
                      help="Y-axis maximum (auto if not specified)")
    parser.add_argument("--ymin", type=float, default=1e-3,
                      help="Y-axis minimum (auto if not specified)")
    
    return parser.parse_args()

def main():
    """Main function to create TTbar stacked MC vs Data plots"""
    
    print("TTbar Stacked MC vs Data Comparison Plotter")
    print("=" * 50)
    
    # Parse command line arguments
    args = parse_arguments()
    
    # Check if data directory exists
    if not os.path.exists(args.data_dir):
        print(f"Error: Data directory {args.data_dir} does not exist!")
        return
    
    # Parse custom bins if provided
    custom_bins = None
    if args.custom_bins:
        try:
            custom_bins = [float(x.strip()) for x in args.custom_bins.split(',')]
            print(f"Using custom bins: {custom_bins}")
            # Store custom bins globally for use in load_histogram
            load_histogram.custom_bins = custom_bins
        except ValueError:
            print(f"Error: Invalid custom bins format '{args.custom_bins}'. Use comma-separated numbers.")
            return
    
    # Parse ranges
    try:
        x_range = [float(x.strip()) for x in args.x_range.split(',')]
        y_range = [float(x.strip()) for x in args.y_range.split(',')]
        # Store xmax for rebinning
        load_histogram.xmax_for_rebin = x_range[-1]
    except ValueError:
        print(f"Error: Invalid range format. Use 'min,max' format.")
        return
    
    # Configuration for the plot
    plot_config = {
        "xRange": x_range,
        "yRange": y_range,
        "xTitle": args.x_title,
        "yTitle": args.y_title,
        "logy": args.logy,
        "rebin": args.rebin,
    }
    
    # Add ymax and ymin if specified
    if args.ymax is not None:
        plot_config["ymax"] = args.ymax
    if args.ymin is not None:
        plot_config["ymin"] = args.ymin
    elif args.logy:
        plot_config["ymin"] = 1e-2  # Default for log scale
    
    # Determine TB line drawing option
    draw_tb_lines = not args.no_tb_lines
    
    print(f"Configuration:")
    print(f"  Data directory: {args.data_dir}")
    print(f"  Histogram: {args.hist_name}")
    print(f"  Output: {args.output_name}")
    print(f"  TB lines: {'Enabled' if draw_tb_lines else 'Disabled'}")
    print(f"  Custom bins: {custom_bins if custom_bins else 'None'}")
    print(f"  X-range: {x_range}")
    print(f"  Y-range: {y_range}")
    print(f"  Log Y: {args.logy}")
    
    try:
        if draw_tb_lines:
            # Get list of TB samples
            tb_samples = get_tb_samples(args.data_dir)
            
            if tb_samples:
                print(f"Found {len(tb_samples)} TB samples: {tb_samples}")
                
                # Create separate plot for each TB sample
                canvases = []
                for tb_sample in tb_samples:
                    print(f"\n{'='*50}")
                    print(f"Creating plot for TB sample: {tb_sample}")
                    print(f"{'='*50}")
                    
                    # Create output name with TB sample
                    tb_output_name = f"{args.output_name}_{tb_sample}"
                    
                    canvas = plot_signal_background_comparison(
                        args.data_dir, 
                        args.hist_name, 
                        plot_config, 
                        tb_output_name,
                        draw_tb_lines=draw_tb_lines,
                        target_tb_sample=tb_sample
                    )
                    
                    if canvas:
                        canvases.append((tb_sample, canvas))
                        print(f"Plot for {tb_sample} created successfully!")
                    else:
                        print(f"Plot creation failed for {tb_sample}!")
                
                if canvases:
                    print(f"\n{'='*50}")
                    print(f"Successfully created {len(canvases)} plots:")
                    for tb_sample, canvas in canvases:
                        print(f"  - {args.output_name}_{tb_sample}.png/pdf")
                    
                    
                    
                    # Clean up canvases
                    for tb_sample, canvas in canvases:
                        canvas.close()
                else:
                    print("All plot creations failed!")
            else:
                print("No TB samples found in the data directory.")
                print("Creating plot without TB lines...")
                
                # Fallback to single plot without TB lines
                canvas = plot_signal_background_comparison(
                    args.data_dir, 
                    args.hist_name, 
                    plot_config, 
                    args.output_name,
                    draw_tb_lines=False
                )
                
                if canvas:
                    print("Plot creation successful!")
                    
                    canvas.close()
                else:
                    print("Plot creation failed!")
        else:
            # Single plot without TB lines
            canvas = plot_signal_background_comparison(
                args.data_dir, 
                args.hist_name, 
                plot_config, 
                args.output_name,
                draw_tb_lines=draw_tb_lines
            )
            
            if canvas:
                print("Plot creation successful!")
                
                canvas.close()
            else:
                print("Plot creation failed!")
            
    except Exception as e:
        print(f"Error in main: {e}")

if __name__ == "__main__":
    main()