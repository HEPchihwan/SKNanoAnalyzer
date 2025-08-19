directory = "/gv0/Users/achihwan/SKNanoRunlog/out/LRSM_TBChannel/2022EE"





import ROOT
import cmsstyle as CMS
import os
import glob
from array import array
import gc

from math import gcd
from functools import reduce

def _gcd_list(vals):
    return reduce(gcd, vals)

def harmonize_uniform_binning(signal_hists, background_hists):
    """
    모든 히스토그램의 nbins, xmin/xmax를 조사해서
    - 범위가 동일하면 nbins의 gcd로 공통 bin 수 결정
    - 각 히스토그램을 그 bin 수로 리빈
    반환: (signal_hists, background_hists, common_bins)
    """
    all_h = {**signal_hists, **background_hists}
    if not all_h:
        return signal_hists, background_hists, None

    nbins_list = [h.GetNbinsX() for h in all_h.values()]
    xmin_list  = [h.GetXaxis().GetXmin() for h in all_h.values()]
    xmax_list  = [h.GetXaxis().GetXmax() for h in all_h.values()]

    if len(set(xmin_list)) != 1 or len(set(xmax_list)) != 1:
        # 범위가 다르면 여기서 변수 binning으로 가거나 그냥 패스
        print("[WARN] X ranges differ. Use variable binning instead.")
        return signal_hists, background_hists, None

    xmin = xmin_list[0]; xmax = xmax_list[0]
    common_bins = _gcd_list(nbins_list)  # e.g. gcd(100,4,...)=4

    if common_bins <= 0:
        print("[WARN] gcd returned non-positive value. Skip.")
        return signal_hists, background_hists, None

    # 각 히스토그램 리빈
    for d in (signal_hists, background_hists):
        for name, hist in d.items():
            nb = hist.GetNbinsX()
            if nb == common_bins:  # 이미 같음
                continue
            if nb % common_bins != 0:
                print(f"[FAIL] {name}: {nb} bins can't be reduced to {common_bins}.")
                return signal_hists, background_hists, None
            factor = nb // common_bins
            new_hist = hist.Rebin(factor, f"{hist.GetName()}_rb")
            new_hist.SetDirectory(0)
            d[name] = new_hist

    print(f"[INFO] Unified to {common_bins} bins, width={(xmax-xmin)/common_bins:.1f} GeV")
    return signal_hists, background_hists, common_bins

def force_variable_edges(signal_hists, background_hists, edges):
    """변수 bin edge 배열로 강제 리빈"""
    nb = len(edges) - 1
    def _rb(h):
        new_h = h.Rebin(nb, f"{h.GetName()}_vb", array('d', edges))
        new_h.SetDirectory(0)
        return new_h
    for d in (signal_hists, background_hists):
        for k, h in list(d.items()):
            d[k] = _rb(h)
    return signal_hists, background_hists

def clamp_negative_bins(hist):
    """빈 내용이 음수면 0으로 클램프 (시각화용)"""
    for i in range(1, hist.GetNbinsX()+1):
        if hist.GetBinContent(i) < 0:
            hist.SetBinContent(i, 0.)
            hist.SetBinError(i, 0.)

# Prevent ROOT from owning Python objects
ROOT.SetOwnership(ROOT.gROOT, False)

SIGNAL_COLOR = ROOT.TColor.GetColor("#e42536")  # Red for signal
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
    def __init__(self, signal_hists, background_hists, config):
        super().__init__()
        
        self.signal_hists = signal_hists
        self.background_hists = background_hists
        self.config = config
        
        # Keep references to prevent garbage collection
        self._objects_to_keep = []
        
        # 빈 정보 자세히 확인
        self._check_histogram_consistency()
        
        # Initialize variables (stack will be built after rebinning)
        self.background_stack = None
        self.total_background = None
        
        # Sort backgrounds by integral (highest first)
        bg_integrals = [(name, hist.Integral()) for name, hist in background_hists.items()]
        bg_integrals.sort(key=lambda x: x[1], reverse=True)
        
        # Separate top 2 and others
        self.top_backgrounds = []
        self.other_backgrounds = []
        
        for i, (name, integral) in enumerate(bg_integrals):
            if i < 5:  # Top 2 backgrounds
                self.top_backgrounds.append((name, background_hists[name]))
            else:  # All other backgrounds
                self.other_backgrounds.append((name, background_hists[name]))
        
        
        
        
        # 리빈 후에도 빈 정보 다시 확인
        if "rebin" in self.config.keys() or len(self.config.get("xRange", [])) > 2:
            print("\n=== After Rebinning ===")
            self._check_histogram_consistency()
        
        # Create "Others" combined histogram if there are more than 2 backgrounds
        self.others_hist = None
        if len(self.other_backgrounds) > 0:
            for i, (name, hist) in enumerate(self.other_backgrounds):
                if self.others_hist is None:
                    self.others_hist = hist.Clone("Others")
                    self.others_hist.SetDirectory(0)  # Detach from ROOT file system
                    self._objects_to_keep.append(self.others_hist)
                else:
                    self.others_hist.Add(hist)
        
        # Build background stack after rebinning
        self._build_background_stack()
        
        # Process signal histograms
        self.total_signal = None
        for name, hist in self.signal_hists.items():
            if not hist.GetSumw2N():  # Sumw2가 아직 호출되지 않았다면
                hist.Sumw2()
            hist.SetLineColor(SIGNAL_COLOR)
            hist.SetLineWidth(3)
            hist.SetFillStyle(0)  # No fill for signal
            hist.SetStats(0)  # Disable statistics box for this histogram
            
            if self.total_signal is None:
                self.total_signal = hist.Clone("total_signal")
                self.total_signal.SetDirectory(0)
                self.total_signal.SetStats(0)  # Disable stats for cloned histogram
                self._objects_to_keep.append(self.total_signal)
            else:
                self.total_signal.Add(hist)
        
        # Create ratio histogram
        if self.total_background is not None and self.total_signal is not None:
            self.ratio = self.total_signal.Clone("signal_ratio")
            self.ratio.SetDirectory(0)
            self.ratio.SetStats(0)  # Disable statistics box for ratio
            self.ratio.Divide(self.total_background)
            self._objects_to_keep.append(self.ratio)
        
        # Set up canvas
        self._setup_canvas()
    
    def _check_histogram_consistency(self):
        """Signal과 Background 히스토그램의 빈 정보를 자세히 비교"""
        print("\n" + "="*60)
        print("HISTOGRAM BINNING CONSISTENCY CHECK")
        print("="*60)
        
        all_hists = {}
        
        # Signal 히스토그램들 수집
        for name, hist in self.signal_hists.items():
            all_hists[f"SIGNAL_{name}"] = hist
        
        # Background 히스토그램들 수집
        for name, hist in self.background_hists.items():
            all_hists[f"BACKGROUND_{name}"] = hist
        
        # 각 히스토그램의 상세 정보 출력
        hist_info = {}
        for name, hist in all_hists.items():
            nbins = hist.GetNbinsX()
            xmin = hist.GetXaxis().GetXmin()
            xmax = hist.GetXaxis().GetXmax()
            bin_width = hist.GetBinWidth(1)
            
            # 처음 3개와 마지막 3개 빈 경계 확인
            bin_edges = []
            for i in range(1, min(4, nbins+1)):  # 처음 3개 빈
                bin_edges.append(hist.GetXaxis().GetBinLowEdge(i))
            if nbins > 6:
                bin_edges.append("...")
                for i in range(max(1, nbins-2), nbins+2):  # 마지막 3개 빈
                    bin_edges.append(hist.GetXaxis().GetBinLowEdge(i))
            
            hist_info[name] = {
                'nbins': nbins,
                'xmin': xmin,
                'xmax': xmax,
                'bin_width': bin_width,
                'entries': hist.GetEntries(),
                'integral': hist.Integral(),
                'bin_edges': bin_edges[:7]  # 처음 몇 개만 저장
            }
            
            print(f"{name:30} | Bins: {nbins:4d} | Range: [{xmin:7.1f}, {xmax:7.1f}] | Width: {bin_width:6.2f} | Entries: {hist.GetEntries():8.0f}")
        
        # 일관성 확인
        print("\n" + "-"*60)
        print("CONSISTENCY ANALYSIS")
        print("-"*60)
        
        # 빈 수 일관성 확인
        signal_bins = [info['nbins'] for name, info in hist_info.items() if name.startswith('SIGNAL_')]
        bg_bins = [info['nbins'] for name, info in hist_info.items() if name.startswith('BACKGROUND_')]
        
        if len(set(signal_bins)) > 1:
            print(f"WARNING: Signal histograms have different bin counts: {set(signal_bins)}")
        if len(set(bg_bins)) > 1:
            print(f"WARNING: Background histograms have different bin counts: {set(bg_bins)}")
        if len(set(signal_bins + bg_bins)) > 1:
            print(f"WARNING: Signal and Background have different bin counts:")
            print(f"  Signal bins: {set(signal_bins)}")
            print(f"  Background bins: {set(bg_bins)}")
        else:
            print(f"✓ All histograms have consistent bin count: {signal_bins[0] if signal_bins else 'N/A'}")
        
        # X 범위 일관성 확인
        all_xmin = [info['xmin'] for info in hist_info.values()]
        all_xmax = [info['xmax'] for info in hist_info.values()]
        
        if len(set(all_xmin)) > 1:
            print(f"WARNING: Histograms have different xmin values: {set(all_xmin)}")
        if len(set(all_xmax)) > 1:
            print(f"WARNING: Histograms have different xmax values: {set(all_xmax)}")
        if len(set(all_xmin)) == 1 and len(set(all_xmax)) == 1:
            print(f"✓ All histograms have consistent X range: [{all_xmin[0]:.1f}, {all_xmax[0]:.1f}]")
        
        # 빈 너비 일관성 확인
        all_bin_widths = [info['bin_width'] for info in hist_info.values()]
        if len(set(all_bin_widths)) > 1:
            print(f"WARNING: Histograms have different bin widths: {set(all_bin_widths)}")
        else:
            print(f"✓ All histograms have consistent bin width: {all_bin_widths[0]:.2f}")
        
        # 빈 경계 세부 확인 (처음 3개 빈)
        print(f"\nDetailed bin edges (first few bins):")
        for name, info in list(hist_info.items())[:3]:  # 처음 3개 히스토그램만
            edges_str = ", ".join([f"{edge:.2f}" if isinstance(edge, (int, float)) else str(edge) for edge in info['bin_edges']])
            print(f"  {name:30}: [{edges_str}]")
        
        # 리빈 가능성 분석
        if signal_bins:
            reference_bins = signal_bins[0]
            possible_rebins = [i for i in range(1, min(reference_bins+1, 21)) if reference_bins % i == 0]
            common_rebins = [r for r in [1, 2, 4, 5, 8, 10, 16, 20] if r in possible_rebins]
            print(f"\nPossible rebin factors for {reference_bins} bins: {possible_rebins}")
            if common_rebins:
                print(f"Recommended rebin factors: {common_rebins}")
        
        print("="*60)
        
    
    
    def _build_background_stack(self):
        """Build the background stack after rebinning is complete"""
        # Create fresh background stack
        self.background_stack = ROOT.THStack("bg_stack", "Background")
        self._objects_to_keep.append(self.background_stack)
        self.total_background = None
        
        # Add "Others" first (bottom of stack)
        if self.others_hist is not None:
            if not self.others_hist.GetSumw2N():  # Sumw2가 아직 호출되지 않았다면
                self.others_hist.Sumw2()
            self.others_hist.SetFillColor(BACKGROUND_COLORS[-1])  # Gray
            self.others_hist.SetLineColor(BACKGROUND_COLORS[-1])
            self.others_hist.SetStats(0)  # Disable statistics box
            
            # 음수 값 제거 (옵션)
            if self.config.get("fix_negative_bins", True):
                for i in range(1, self.others_hist.GetNbinsX() + 1):
                    if self.others_hist.GetBinContent(i) < 0:
                        self.others_hist.SetBinContent(i, 0)
                        self.others_hist.SetBinError(i, 0)
            
            self.background_stack.Add(self.others_hist)
            
            self.total_background = self.others_hist.Clone("total_background")
            self.total_background.SetDirectory(0)
            self.total_background.SetStats(0)  # Disable stats for total background
            self._objects_to_keep.append(self.total_background)
        
        # Add top 2 backgrounds (from lowest to highest for proper stacking)
        for i, (name, hist) in enumerate(reversed(self.top_backgrounds)):
            if not hist.GetSumw2N():  # Sumw2가 아직 호출되지 않았다면
                hist.Sumw2()
            hist.SetFillColor(BACKGROUND_COLORS[len(self.top_backgrounds)-1-i])
            hist.SetLineColor(BACKGROUND_COLORS[len(self.top_backgrounds)-1-i])
            hist.SetStats(0)  # Disable statistics box for each background
            
            # 음수 값 제거 (옵션)
            if self.config.get("fix_negative_bins", True):
                for j in range(1, hist.GetNbinsX() + 1):
                    if hist.GetBinContent(j) < 0:
                        hist.SetBinContent(j, 0)
                        hist.SetBinError(j, 0)
            
            self.background_stack.Add(hist)
            
            if self.total_background is None:
                self.total_background = hist.Clone("total_background")
                self.total_background.SetDirectory(0)
                self.total_background.SetStats(0)  # Disable stats
                self._objects_to_keep.append(self.total_background)
            else:
                self.total_background.Add(hist)
    
    def _setup_canvas(self):
        """Setup the canvas with proper styling"""
        # Disable statistics box globally
        ROOT.gStyle.SetOptStat(0)  # Turn off statistics box
        ROOT.gStyle.SetOptTitle(0)  # Turn off title box
        
        # Set y-axis range
        if "ymax" in self.config.keys():
            ymax = self.config["ymax"]
            print(f"Using explicit ymax from config: {ymax}")
        else:
            bg_max = self.total_background.GetMaximum() if self.total_background else 0
            sig_max = self.total_signal.GetMaximum() if self.total_signal else 0
            ymax = max(bg_max, sig_max) * 1.5
            print(f"Auto-calculated ymax: {ymax} (bg_max: {bg_max}, sig_max: {sig_max})")
        
        ymin = self.config.get("ymin", 0.)
        
        # Apply log scale defaults if not specified by user
        if "logy" in self.config.keys() and self.config['logy']:
            if "ymin" not in self.config.keys():
                ymin = 1e-3
            if "ymax" not in self.config.keys():
                ymax = max(ymax, 1000)
        
        print(f"Final y-range: [{ymin}, {ymax}]")
        
        # Set CMS style
        CMS.SetEnergy(13.6)
        CMS.SetLumi("26.6717 fb^{-1}")
        CMS.SetExtraText("Simulation Preliminary")
        
        # Create canvas with unique name to avoid conflicts
        canvas_name = f"canvas_{id(self)}"
        self.canv = ROOT.TCanvas(canvas_name, canvas_name, 800, 800)
        ROOT.SetOwnership(self.canv, False)  # Python keeps ownership
        self._objects_to_keep.append(self.canv)
        
        # Create two pads manually with better control
        pad1_name = f"pad1_{id(self)}"
        pad2_name = f"pad2_{id(self)}"
        
        pad1 = ROOT.TPad(pad1_name, pad1_name, 0, 0.3, 1, 1.0)
        pad1.SetBottomMargin(0.02)  # Small bottom margin for upper pad
        pad1.SetTopMargin(0.08)     # Top margin for CMS label
        pad1.SetLeftMargin(0.12)    # Left margin for y-axis label
        pad1.SetRightMargin(0.05)   # Right margin
        ROOT.SetOwnership(pad1, False)
        self._objects_to_keep.append(pad1)
        
        pad2 = ROOT.TPad(pad2_name, pad2_name, 0, 0.0, 1, 0.3)
        pad2.SetTopMargin(0.02)     # Small top margin for lower pad
        pad2.SetBottomMargin(0.35)  # Large bottom margin for x-axis label
        pad2.SetLeftMargin(0.12)    # Match upper pad
        pad2.SetRightMargin(0.05)   # Match upper pad
        ROOT.SetOwnership(pad2, False)
        self._objects_to_keep.append(pad2)
        
        # Apply log scale if needed
        if "logy" in self.config.keys() and self.config['logy']:
            pad1.SetLogy()
            pad2.SetLogy()
        pad1.Draw()
        pad2.Draw()
        
        self.pad1 = pad1
        self.pad2 = pad2
        
        # Set up axis ranges manually
        self.xmin = self.config["xRange"][0]
        self.xmax = self.config["xRange"][-1]
        self.ymin = ymin
        self.ymax = ymax
        self.ratio_ymin = self.config["yRange"][0]
        self.ratio_ymax = self.config["yRange"][1]
        
        # Adjust legend size based on number of entries
        n_entries = len(self.top_backgrounds) + (1 if self.others_hist else 0) + (1 if self.total_signal else 0)
        leg_height = 0.05 * (n_entries + 1)
        self.leg = ROOT.TLegend(0.75, 0.9 - leg_height, 0.92, 0.9)
        self.leg.SetTextSize(0.02)
        self.leg.SetBorderSize(0)
        self.leg.SetFillStyle(0)
        ROOT.SetOwnership(self.leg, False)
        self._objects_to_keep.append(self.leg)
    
    def draw(self):
        """Draw the signal and background comparison"""
        
        # =====================================
        # 상단 플롯 (Upper pad) - Signal/Background
        # =====================================
        self.pad1.cd()
        
        # First draw the background stack to establish the plot
        if self.background_stack:
            self.background_stack.Draw("HIST")
            
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
        else:
            # If no background stack, create dummy histogram
            dummy = ROOT.TH1F(f"dummy_{id(self)}", "", 100, self.xmin, self.xmax)
            dummy.SetMinimum(self.ymin)
            dummy.SetMaximum(self.ymax)
            dummy.GetYaxis().SetTitle(self.config["yTitle"])
            dummy.GetYaxis().SetTitleSize(0.05)
            dummy.GetYaxis().SetTitleOffset(1.2)
            dummy.GetYaxis().SetLabelSize(0.04)
            dummy.GetXaxis().SetLabelSize(0)
            dummy.SetDirectory(0)
            ROOT.SetOwnership(dummy, False)
            self._objects_to_keep.append(dummy)
            dummy.Draw()
        
        # Draw signal on top
        if self.total_signal:
            self.total_signal.Draw("HIST SAME")
        
        # Add reference line at y=1 for upper plot (optional)
        if self.config.get("show_unity_line", False):
            unity_line_upper = ROOT.TLine(self.xmin, 1., self.xmax, 1.)
            unity_line_upper.SetLineStyle(2)  # Dashed line
            unity_line_upper.SetLineColor(ROOT.kGray+2)
            unity_line_upper.SetLineWidth(1)
            ROOT.SetOwnership(unity_line_upper, False)
            self._objects_to_keep.append(unity_line_upper)
            unity_line_upper.Draw()
        
        # Add legend entries in the correct order (for stacked histograms)
        for name, hist in self.top_backgrounds:
            clean_name = name.replace("_powheg", "").replace("_pythia", "")
            self.leg.AddEntry(hist, clean_name, "F")
        
        if self.others_hist:
            self.leg.AddEntry(self.others_hist, "Others", "F")
        
        if self.total_signal:
            signal_label = getattr(self, 'signal_mass_point', 'Signal')
            self.leg.AddEntry(self.total_signal, signal_label, "L")
        
        self.leg.Draw()
        
        # Add CMS labels manually - MOVED HIGHER
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextSize(0.045)
        latex.SetTextFont(62)
        latex.DrawLatex(0.12, 0.93, "CMS")  # Moved from 0.91 to 0.95
        
        latex.SetTextSize(0.035)
        latex.SetTextFont(52)
        latex.DrawLatex(0.20, 0.93, "Simulation Preliminary")  # Moved from 0.91 to 0.95
        
        latex.SetTextFont(42)
        latex.DrawLatex(0.73, 0.93, "26.6717 fb^{-1} (13.6 TeV)")  # Moved from 0.70, 0.91 to 0.65, 0.95
        
        ROOT.SetOwnership(latex, False)
        self._objects_to_keep.append(latex)
        
        # Redraw axis to make sure it's on top
        self.pad1.RedrawAxis()
        
        # =====================================
        # 하단 플롯 (Lower pad) - Ratio plot
        # =====================================
        self.pad2.cd()
        
        if hasattr(self, 'ratio') and self.ratio:
            # Create dummy histogram for ratio pad axis
            nb = self.total_background.GetNbinsX() if self.total_background else 100
            ratio_dummy = ROOT.TH1F(f"ratio_dummy_{id(self)}", "", nb, self.xmin, self.xmax)
            ratio_dummy.SetMinimum(self.ratio_ymin)
            ratio_dummy.SetMaximum(self.ratio_ymax)
            ratio_dummy.GetXaxis().SetTitle(self.config["xTitle"])
            ratio_dummy.GetYaxis().SetTitle("Signal / Background")
            
            # Adjust font sizes for ratio pad
            ratio_dummy.GetXaxis().SetTitleSize(0.12)
            ratio_dummy.GetXaxis().SetTitleOffset(1.0)
            ratio_dummy.GetXaxis().SetLabelSize(0.10)
            ratio_dummy.GetYaxis().SetTitleSize(0.08)
            ratio_dummy.GetYaxis().SetTitleOffset(0.5)
            ratio_dummy.GetYaxis().SetLabelSize(0.05)
            ratio_dummy.GetYaxis().SetNdivisions(505)
            ratio_dummy.SetDirectory(0)
            
            ROOT.SetOwnership(ratio_dummy, False)
            self._objects_to_keep.append(ratio_dummy)
            ratio_dummy.Draw()
            
            # Draw horizontal line at y=1
            line = ROOT.TLine(self.xmin, 1., self.xmax, 1.)
            line.SetLineStyle(2)  # Dashed line
            line.SetLineColor(ROOT.kBlack)
            ROOT.SetOwnership(line, False)
            self._objects_to_keep.append(line)
            line.Draw()
            
            # Draw ratio histogram
            self.ratio.SetLineColor(SIGNAL_COLOR)
            self.ratio.SetLineWidth(2)
            self.ratio.SetMarkerStyle(20)
            self.ratio.SetMarkerSize(0.8)
            self.ratio.Draw("HIST SAME")
            
            # Redraw axis for ratio pad too
            self.pad2.RedrawAxis()
        
        # Update canvas
        self.canv.Update()
    
    def save_as(self, filename):
        """Safe save method"""
        self.canv.SaveAs(filename)
    
    def close(self):
        """Safely close canvas and clean up"""
        if hasattr(self, 'canv') and self.canv:
            self.canv.Close()
        
        # Clear all references
        self._objects_to_keep.clear()
        
        # Force garbage collection
        gc.collect()
    
    def __del__(self):
        """Destructor to clean up resources"""
        try:
            self.close()
        except:
            pass

def load_histograms(file_path, hist_name, systematic="Central", silent=False):
    """Load histogram from ROOT file with proper memory management"""
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
    hist_clone.SetDirectory(0)  # Detach from file
    ROOT.SetOwnership(hist_clone, False)  # Python owns the object
    root_file.Close()
    
    return hist_clone

def check_histogram_binning(file_path, hist_name, systematic="Central"):
    """
    히스토그램의 비닝 정보를 확인하는 유틸리티 함수
    """
    print(f"Checking histogram: {hist_name} in {file_path}")
    
    root_file = ROOT.TFile.Open(file_path)
    if not root_file or root_file.IsZombie():
        print(f"Error: Cannot open file {file_path}")
        return None
    
    directory = root_file.Get(systematic)
    if not directory:
        print(f"Error: Cannot find directory {systematic}")
        root_file.Close()
        return None
    
    hist = directory.Get(hist_name)
    if not hist:
        print(f"Error: Cannot find histogram {hist_name}")
        root_file.Close()
        return None
    
    # 히스토그램 정보 출력
    nbins = hist.GetNbinsX()
    xmin = hist.GetXaxis().GetXmin()
    xmax = hist.GetXaxis().GetXmax()
    bin_width = hist.GetBinWidth(1)
    
    print(f"  Total number of bins: {nbins}")  # 총 빈수 명시적 출력
    print(f"  X-axis range: [{xmin:.1f}, {xmax:.1f}]")
    print(f"  Bin width: {bin_width:.2f}")
    print(f"  Entries: {hist.GetEntries():.0f}")
    print(f"  Integral: {hist.Integral():.2f}")
    
    # 가능한 리빈 팩터들 계산
    possible_rebins = [i for i in range(1, min(nbins+1, 21)) if nbins % i == 0]
    print(f"  Possible rebin factors (up to 20): {possible_rebins}")
    
    # 추가: 일반적으로 많이 사용되는 리빈 팩터들 강조
    common_rebins = [r for r in [1, 2, 4, 5, 8, 10, 16, 20] if r in possible_rebins]
    if common_rebins:
        print(f"  Commonly used rebin factors available: {common_rebins}")
    
    root_file.Close()
    return {
        'nbins': nbins,
        'xmin': xmin,
        'xmax': xmax,
        'bin_width': bin_width,
        'entries': hist.GetEntries(),
        'integral': hist.Integral(),
        'possible_rebins': possible_rebins
    }

def compare_histogram_details(data_dir, hist_name, systematic="Central"):
    """
    Signal과 Background 히스토그램들의 상세 비교 분석
    """
    root_files = glob.glob(os.path.join(data_dir, "*.root"))
    
    print("="*80)
    print(f"DETAILED HISTOGRAM COMPARISON: {hist_name}")
    print("="*80)
    
    signal_hists = {}
    background_hists = {}
    
    # 파일들을 분류하고 히스토그램 로드
    for file_path in root_files:
        filename = os.path.basename(file_path)
        is_signal = filename.startswith("TB")
        hist = load_histograms(file_path, hist_name, systematic, silent=True)
        
        if hist is None:
            continue
            
        if is_signal:
            signal_hists[filename] = hist
        else:
            background_hists[filename] = hist
    
    print(f"Found {len(signal_hists)} signal files and {len(background_hists)} background files")
    print()
    
    # Signal 히스토그램들 상세 분석
    if signal_hists:
        print("SIGNAL HISTOGRAMS:")
        print("-" * 80)
        for i, (filename, hist) in enumerate(signal_hists.items()):
            print(f"{i+1:2d}. {filename}")
            _print_histogram_details(hist, "    ")
            if i == 2:  # 처음 3개만 상세히 출력
                remaining = len(signal_hists) - 3
                if remaining > 0:
                    print(f"    ... and {remaining} more signal files")
                break
        print()
    
    # Background 히스토그램들 상세 분석
    if background_hists:
        print("BACKGROUND HISTOGRAMS:")
        print("-" * 80)
        for i, (filename, hist) in enumerate(background_hists.items()):
            print(f"{i+1:2d}. {filename}")
            _print_histogram_details(hist, "    ")
            if i == 2:  # 처음 3개만 상세히 출력
                remaining = len(background_hists) - 3
                if remaining > 0:
                    print(f"    ... and {remaining} more background files")
                break
        print()
    
    # 비교 분석
    _analyze_histogram_differences(signal_hists, background_hists)
    
    return signal_hists, background_hists

def _print_histogram_details(hist, indent=""):
    """히스토그램의 상세 정보를 출력"""
    nbins = hist.GetNbinsX()
    xmin = hist.GetXaxis().GetXmin()
    xmax = hist.GetXaxis().GetXmax()
    bin_width = hist.GetBinWidth(1)
    entries = hist.GetEntries()
    integral = hist.Integral()
    
    print(f"{indent}Bins: {nbins:4d} | Range: [{xmin:7.1f}, {xmax:7.1f}] | Width: {bin_width:6.2f}")
    print(f"{indent}Entries: {entries:8.0f} | Integral: {integral:10.2f}")
    
    # 처음 5개 빈의 경계와 내용 확인
    print(f"{indent}First 5 bin edges: ", end="")
    for i in range(1, min(6, nbins+1)):
        edge = hist.GetXaxis().GetBinLowEdge(i)
        print(f"{edge:6.1f}", end=" ")
    print()
    
    # 비어있지 않은 빈의 개수
    non_empty_bins = sum(1 for i in range(1, nbins+1) if hist.GetBinContent(i) > 0)
    print(f"{indent}Non-empty bins: {non_empty_bins}/{nbins} ({100*non_empty_bins/nbins:.1f}%)")
    print()

def _analyze_histogram_differences(signal_hists, background_hists):
    """Signal과 Background 히스토그램 간의 차이점 분석"""
    print("COMPARATIVE ANALYSIS:")
    print("-" * 80)
    
    all_hists = {}
    all_hists.update({f"SIGNAL_{k}": v for k, v in signal_hists.items()})
    all_hists.update({f"BACKGROUND_{k}": v for k, v in background_hists.items()})
    
    if not all_hists:
        print("No histograms to analyze!")
        return
    
    # 기본 통계
    stats = {}
    for name, hist in all_hists.items():
        stats[name] = {
            'nbins': hist.GetNbinsX(),
            'xmin': hist.GetXaxis().GetXmin(),
            'xmax': hist.GetXaxis().GetXmax(),
            'bin_width': hist.GetBinWidth(1),
            'entries': hist.GetEntries(),
            'integral': hist.Integral()
        }
    
    # 빈 수 분석
    signal_bins = [stats[k]['nbins'] for k in stats if k.startswith('SIGNAL_')]
    bg_bins = [stats[k]['nbins'] for k in stats if k.startswith('BACKGROUND_')]
    
    print(f"Signal bin counts: {set(signal_bins)} (total: {len(signal_bins)} files)")
    print(f"Background bin counts: {set(bg_bins)} (total: {len(bg_bins)} files)")
    
    if len(set(signal_bins + bg_bins)) == 1:
        print("✓ All histograms have the same number of bins")
        common_bins = signal_bins[0] if signal_bins else bg_bins[0]
        possible_rebins = [i for i in range(1, min(common_bins+1, 21)) if common_bins % i == 0]
        recommended = [r for r in [1, 2, 4, 5, 8, 10, 16, 20] if r in possible_rebins]
        print(f"  Possible rebin factors: {possible_rebins}")
        print(f"  Recommended: {recommended}")
    else:
        print("⚠ WARNING: Histograms have different bin counts!")
        print("  This may cause issues when combining or comparing histograms")
    
    # X 범위 분석
    all_xmin = [stats[k]['xmin'] for k in stats]
    all_xmax = [stats[k]['xmax'] for k in stats]
    
    if len(set(all_xmin)) == 1 and len(set(all_xmax)) == 1:
        print(f"✓ All histograms have the same X range: [{all_xmin[0]:.1f}, {all_xmax[0]:.1f}]")
    else:
        print("⚠ WARNING: Histograms have different X ranges!")
        print(f"  X-min values: {set(all_xmin)}")
        print(f"  X-max values: {set(all_xmax)}")
    
    # 빈 너비 분석
    all_widths = [stats[k]['bin_width'] for k in stats]
    if len(set(all_widths)) == 1:
        print(f"✓ All histograms have the same bin width: {all_widths[0]:.2f}")
    else:
        print("⚠ WARNING: Histograms have different bin widths!")
        print(f"  Bin widths: {set(all_widths)}")
    
    # 통계 요약
    print(f"\nEntries range: {min(stats[k]['entries'] for k in stats):.0f} - {max(stats[k]['entries'] for k in stats):.0f}")
    print(f"Integral range: {min(stats[k]['integral'] for k in stats):.2f} - {max(stats[k]['integral'] for k in stats):.2f}")
    
    print("="*80)
    """
    디렉토리 내의 모든 ROOT 파일에서 히스토그램 정보 확인
    """
    root_files = glob.glob(os.path.join(data_dir, "*.root"))
    
    print(f"Checking histogram '{hist_name}' in {len(root_files)} files:")
    print("=" * 60)
    
    all_info = {}
    for file_path in root_files[:5]:  # 처음 5개 파일만 확인
        filename = os.path.basename(file_path)
        info = check_histogram_binning(file_path, hist_name, systematic)
        if info:
            all_info[filename] = info
        print("-" * 40)
    
    # 공통 가능한 리빈 팩터 찾기
    if all_info:
        common_rebins = None
        for filename, info in all_info.items():
            if common_rebins is None:
                common_rebins = set(info['possible_rebins'])
            else:
                common_rebins = common_rebins.intersection(set(info['possible_rebins']))
        
        print(f"Common possible rebin factors across all files: {sorted(list(common_rebins))}")
    
    return all_info
    """Load histogram from ROOT file with proper memory management"""
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
    hist_clone.SetDirectory(0)  # Detach from file
    ROOT.SetOwnership(hist_clone, False)  # Python owns the object
    root_file.Close()
    
    return hist_clone

def plot_signal_background(data_dir, hist_name, config, output_name="signal_background_plot", systematic="Central"):
    """
    Main function to create signal vs background plots with safe memory management
    """
    
    # Get all ROOT files in directory
    root_files = glob.glob(os.path.join(data_dir, "*.root"))
    
    signal_hists = {}
    background_hists = {}
    
    # Separate signal and background files
    for file_path in root_files:
        filename = os.path.basename(file_path)
        is_signal = filename.startswith("TB")
        hist = load_histograms(file_path, hist_name, systematic, silent=not is_signal)
        
        if hist is None:
            continue
            
        if is_signal:
            signal_hists[filename.replace(".root", "")] = hist
        else:
            background_hists[filename.replace(".root", "")] = hist
    
    if not signal_hists:
        print("Warning: No signal histograms found!")
    if not background_hists:
        print("Warning: No background histograms found!")
        # (중략) signal_hists/background_hists 채운 뒤

    # ① config["rebin"]이 있더라도 여기서 일단 지움 (충돌 방지)
    if "rebin" in config:
        del config["rebin"]

    # ② 모든 히스토그램을 공통 bin으로 강제
    
    edges = [0, 1000, 2000,2500,3000,3500, 4000,4500,5000, 5500,6000,6500,7000,7500, 8000]
    signal_hists, background_hists = force_variable_edges(signal_hists, background_hists, edges)

    # ③ 음수 빈 클램프
    for h in signal_hists.values():
        clamp_negative_bins(h)
    for h in background_hists.values():
        clamp_negative_bins(h)

    
    # Create and draw the plot
    canvas = None
    try:
        canvas = SignalBackgroundCanvas(signal_hists, background_hists, config)
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

def plot_individual_signals(data_dir, hist_name, config, systematic="Central"):
    root_files = glob.glob(os.path.join(data_dir, "*.root"))

    signal_files = []
    background_hists = {}

    for file_path in root_files:
        filename = os.path.basename(file_path)
        is_signal = filename.startswith("TB")
        hist = load_histograms(file_path, hist_name, systematic, silent=not is_signal)
        if hist is None: 
            continue
        if is_signal:
            signal_files.append((filename, hist))
        else:
            background_hists[filename.replace(".root", "")] = hist

    if not signal_files:
        print("Warning: No signal files found!")
        return []
    if not background_hists:
        print("Warning: No background histograms found!")
        return []

    # 0) rebin 키 제거
    config.pop("rebin", None)

    # 1) 공통 edges 정의 (또는 gcd 방식)
    edges = [0, 1000, 2000,2500,3000,3500, 4000,4500,5000, 5500,6000,6500,7000,7500, 8000]
    nb = len(edges) - 1
    arr_edges = array('d', edges)

    # 2) 배경은 한 번만 리빈 -> 복사 dict
    rb_bkg = {}
    for n, h in background_hists.items():
        new_h = h.Rebin(nb, f"{n}_vb", arr_edges)
        new_h.SetDirectory(0)
        clamp_negative_bins(new_h)
        rb_bkg[n] = new_h

    canvases = []

    # 3) 각 시그널마다 개별 리빈 + 캔버스
    for filename, sig_hist in signal_files:
        mass_point = filename.replace("TBChannel_", "").replace(".root", "")

        sig_dict = {}
        new_sig = sig_hist.Rebin(nb, f"{filename}_vb", arr_edges)
        new_sig.SetDirectory(0)
        clamp_negative_bins(new_sig)
        sig_dict[filename.replace(".root","")] = new_sig

        try:
            print(f"Creating individual plot for mass point: {mass_point}")
            canvas = SignalBackgroundCanvas(sig_dict, rb_bkg, config)
            canvas.signal_mass_point = mass_point
            canvas.draw()

            out = f"{hist_name}_{mass_point}_vs_backgrounds"
            canvas.save_as(f"{out}.pdf")
            canvases.append(canvas)
        except Exception as e:
            print(f"Error creating plot for {mass_point}: {e}")
            if 'canvas' in locals():
                canvas.close()

    return canvases

import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))



def plot_individual_mass_points():
    """Create individual plots for each signal mass point"""
    
    # Path to your ROOT files
    data_directory = "/gv0/Users/achihwan/SKNanoRunlog/out/LRSM_TBChannel/2022EE"
    
    # Configuration for individual mass point plots
    plot_configs = {
        "WRMass_Central": {
            "xRange": [0, 8000],
            "yRange": [0.001, 1e+3],
            "xTitle": "M_{W_{R}} [GeV]",
            "yTitle": "Events / bin",
            "logy": True,
            "rebin": 10,
            "canvas_width": 700,   # Even smaller size
            "canvas_height": 700
            ,"show_unity_line": True
        }
    }
    
    # Create individual plots for each histogram
    for hist_name, config in plot_configs.items():
        print(f"Creating individual signal plots for: {hist_name}")
        try:
            canvases = plot_individual_signals(data_directory, hist_name, config)
            print(f"Created {len(canvases)} individual plots for {hist_name}")
        except Exception as e:
            print(f"Error creating individual plots for {hist_name}: {e}")
    
    print("Finished creating individual signal mass point plots!")


if __name__ == "__main__":
    print("LRSM TBChannel Signal vs Background Plotter Example")
    print("=" * 60)
    
    # Check if data directory exists
    data_dir = directory
    if not os.path.exists(data_dir):
        print(f"Error: Data directory {data_dir} does not exist!")
        sys.exit(1)
    
    # Ask user which type of plots to create
    print("Choose plotting option:")
    print("1. Combined signal vs background plots (all signals together)")
    print("2. Individual signal mass point plots (separate plot for each mass)")
    print("3. Both")
    plot_individual_mass_points()
    
    