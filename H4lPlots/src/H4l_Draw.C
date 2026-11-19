#include <ROOT/RDF/HistoModels.hxx>
#include <ROOT/RDF/InterfaceUtils.hxx>
#include <Rtypes.h>
#include <RtypesCore.h>
#include <TAxis.h>
#include <TBox.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraphAsymmErrors.h>
#include <TH1.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMathBase.h>
#include <TPad.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#pragma cling add_include_path(                                                \
    "/eos/user/m/mtokarev/work_area/currentCMSSW/src")

#include "PlotScripts/H4lPlots/interface/Logger.hxx"
#include "PlotScripts/H4lPlots/src/ChannelNames.hxx"
#include "PlotScripts/H4lPlots/src/H4l_Common.C"
#include "PlotScripts/H4lPlots/src/Paths.hxx"
#include "PlotScripts/H4lPlots/src/Style.hxx"
#include "PlotScripts/H4lPlots/src/Variables.C"

#define N_SR_HISTOS 5
#define N_CR_HISTOS 5

const std::map<std::string, float> LumiPerPeriod = {
    {"2024", 109.95},
    {"2025", 110.63},
}; // Place-holder value

const std::map<std::string, std::string> latexVariableNames = {
    {"ZZMass", "m_{#it{4l}}"},
    {"ZZMass_4e", "m_{#it{4e}}"},
    {"ZZMass_2e2mu", "m_{#it{2e2#mu}}"},
    {"ZZMass_4mu", "m_{#it{4#mu}}"},
    {"Z1Mass", "m_{Z_{1}}"},
    {"Z2Mass", "m_{Z_{2}}"},
};

const std::unordered_map<std::string, std::vector<std::string>>
    dataSetPerPeriod = {{"Run3",
                         {"2022", "2022EE", "2023preBPix", "2023postBPix",
                          "2024", "2025", "2026"}}};

struct blindInterval {
  Double_t low, up;
  bool isIn(Double_t value) { return this->low <= value && value <= this->up; }
};

blindInterval blindIntervals[] = {
    blindInterval{.low = 105, .up = 140}, blindInterval{.low = 300, .up = 500},
    blindInterval{.low = 500, .up = TMath::Infinity()}};

TH1D *getZX(TH1D *model, FinalStateFlavor fsFlav) {

  int nEntries(10000);
  float binDown(70), binUp(3000), lumiRun3(62 * 1000);
  TF1 *distro4e =
      new TF1("f_4e_comb", "TMath::Landau(x, [0], [1])", binDown, binUp);
  TF1 *distro4mu =
      new TF1("f_4mu_comb", "TMath::Landau(x, [0], [1])", binDown, binUp);
  TF1 *distro2e2mu = new TF1(
      "f_2e2mu_comb",
      "[0] * TMath::Landau(x, [1], [2]) + [3] * TMath::Landau(x, [4], [5])",
      binDown, binUp);

  distro4e->SetParameters(141.9, 21.3);
  distro4mu->SetParameters(130.4, 15.6);
  distro2e2mu->SetParameters(0.45, 131.1, 18.1, 0.55, 133.8, 18.9);

  float yield4eRun3 = 24.256 / lumiRun3;
  float yield4muRun3 = 73.199 / lumiRun3;
  float yield2e2muRun3 = 86.14 / lumiRun3;

  TH1D *histo4e = (TH1D *)model->Clone("ZX_4e");
  histo4e->Clear();
  TH1D *histo4mu = (TH1D *)histo4e->Clone("ZX_4mu");
  TH1D *histo2e2mu = (TH1D *)histo4e->Clone("ZX_2e2mu");

  histo4e->FillRandom(distro4e, nEntries);
  histo4mu->FillRandom(distro4mu, nEntries);
  histo2e2mu->FillRandom(distro2e2mu, nEntries);

  histo4e->Scale(yield4eRun3 / histo4e->Integral());
  histo4mu->Scale(yield4muRun3 / histo4mu->Integral());
  histo2e2mu->Scale(yield2e2muRun3 / histo2e2mu->Integral());

  TH1D *outHisto;
  switch (fsFlav) {
  case FinalStateFlavor::L4:
    outHisto = (TH1D *)histo4e->Clone("ZX_tot");
    outHisto->Add(histo4mu);
    outHisto->Add(histo2e2mu);
    break;
  case FinalStateFlavor::E4:
    outHisto = histo4e;
    break;
  case FinalStateFlavor::Mu4:
    outHisto = histo4mu;
    break;
  case FinalStateFlavor::E2Mu2:
    outHisto = histo2e2mu;
    break;
  }

  return outHisto;
}

template <size_t N>
TH1D *agregateHistogram(TFile *file, const char *variable, int binning,
                        unitIdx unit,
                        const std::array<const char *, N> &samples,
                        const std::string &sampleBundleName,
                        FinalStateFlavor fsFlav) {

  std::string fsString("");

  if (fsFlav != FinalStateFlavor::L4) {
    fsString = makeFSString(fsFlav);
  }

  TH1D *histogram = nullptr;

  for (auto sample : samples) {
    TString inHistName = makeHistoName(variable, binning, unit, sample, fsFlav);
    LOG("[agregateHistogram]\tExtracting the histogram "
            << inHistName << " for the variable " << variable,
        3);

    TH1D *currHist = (TH1D *)file->Get(inHistName);

    if (!currHist)
      throw std::runtime_error("Histogram not found: " + inHistName +
                               " in file " + file->GetPath());

    if (!histogram) { // Initialize first histogram
      histogram = (TH1D *)currHist->Clone(sampleBundleName.c_str());
    } else {
      histogram->Add(currHist);
    }
  }
  return histogram;
}

std::array<TH1D *, N_SR_HISTOS> getMCHistosSR(TFile *inFile, Float_t lumi,
                                              const char *variable,
                                              FinalStateFlavor fsFlav,
                                              int binning, unitIdx unit) {

  auto mkBundleName = [variable, fsFlav, binning, unit](
                          const std::string &sampleName) -> const std::string {
    return TString::Format("%s_%i%s_%s%s", variable, binning,
                           getUnitNameForTitle(unit).c_str(),
                           sampleName.c_str(), makeFSString(fsFlav).c_str())
        .Data();
  };

  auto mkHistoSum = [inFile, variable, binning, unit, fsFlav, mkBundleName](
                        const std::string &sampleName, const auto &dataSet) {
    return agregateHistogram(inFile, variable, binning, unit, dataSet,
                             mkBundleName(sampleName), fsFlav);
  };

  TH1D *signal = mkHistoSum("signal", MC_SignalH125);
  signal->Scale(lumi * 1000.);
  signal->SetLineColor(TColor::GetColor(204, 00, 00));
  signal->SetFillColor(TColor::GetColor(255, 155, 155));

  TH1D *EW = mkHistoSum("EW", MC_BkgEW);
  EW->Scale(lumi * 1000.);
  EW->SetLineColor(TColor::GetColor(0, 0, 153));
  EW->SetFillColor(TColor::GetColor(3, 49, 185));

  TH1D *ZZTo4l = mkHistoSum("qqZZ", MC_BkgQQZZ);
  ZZTo4l->Scale(lumi * 1000.);
  ZZTo4l->SetLineColor(TColor::GetColor(0, 0, 153));
  ZZTo4l->SetFillColor(TColor::GetColor(155, 204, 255));

  TH1D *ggTo = mkHistoSum("ggZZ", MC_BkgGGTo);
  ggTo->Scale(lumi * 1000.);
  ggTo->SetLineColor(TColor::GetColor(0, 0, 153));
  ggTo->SetFillColor(TColor::GetColor(75, 120, 255));

  TH1D *DY = mkHistoSum("DY", MC_BkgDY);
  DY->Scale(lumi * 1000.);
  DY->SetLineColor(TColor::GetColor(0, 0, 153));
  DY->SetFillColor(kOrange);

  TH1D *TT = mkHistoSum("ttTo", MC_BkgTT);
  TT->Scale(lumi * 1000.);

  TH1D *ZX = getZX(ZZTo4l, fsFlav);
  ZX->Scale(lumi * 1000.);
  ZX->SetLineColor(TColor::GetColor(0, 51, 0));
  ZX->SetFillColor(TColor::GetColor(102, 153, 102));

  return {ZX, /*DY, TT, */ EW, ggTo, ZZTo4l, signal};
}

namespace MCSRHistoIdx {
enum MCHistoIdx {
  ZXIdx,
  EWIdx,
  ggToIdx,
  ZZTo4lIdx,
  /*DYIdx, TTIdx,*/
  signalIdx
};
}

std::array<TH1D *, N_CR_HISTOS> getMCHistosCR(TFile *inFile, Float_t lumi,
                                              const char *variable,
                                              FinalStateFlavor fsFlav,
                                              int binning, unitIdx unit) {
  auto mkBundleName = [variable, fsFlav, binning, unit](
                          const std::string &sampleName) -> const std::string {
    return TString::Format("%s_%i%s_%s%s", variable, binning,
                           getUnitNameForTitle(unit).c_str(),
                           sampleName.c_str(), makeFSString(fsFlav).c_str())
        .Data();
  };

  auto mkHistoSum = [inFile, variable, binning, unit, fsFlav, mkBundleName](
                        const std::string &sampleName, const auto &dataSet) {
    return agregateHistogram(inFile, variable, binning, unit, dataSet,
                             mkBundleName(sampleName), fsFlav);
  };
  TH1D *DY = mkHistoSum("DY_CR", MC_BkgDY);
  DY->Scale(lumi * 1000.);
  DY->SetLineColor(TColor::GetColor(0, 0, 153));
  DY->SetFillColor(kOrange);

  TH1D *ggTo = mkHistoSum("ggTo_CR", MC_BkgGGTo);
  ggTo->Scale(lumi * 1000.);
  ggTo->SetLineColor(TColor::GetColor(0, 0, 153));
  ggTo->SetFillColor(TColor::GetColor(102, 255, 102));

  TH1D *TTTo = mkHistoSum("TTTo_CR", MC_BkgTTTo);
  TTTo->Scale(lumi * 1000.);
  TTTo->SetLineColor(TColor::GetColor(0, 0, 153));
  TTTo->SetFillColor(TColor::GetColor(255, 153, 102));

  TH1D *WZTo = mkHistoSum("WZTo_CR", std::array<const char *, 1>{"WZto3LNu"});
  WZTo->Scale(lumi * 1000.);
  WZTo->SetLineColor(TColor::GetColor(0, 0, 153));
  WZTo->SetFillColor(TColor::GetColor(175, 204, 255));

  TH1D *ZZTo = mkHistoSum("ZZTo_CR", MC_BkgQQZZ);
  ZZTo->Scale(lumi * 1000.);
  ZZTo->SetLineColor(TColor::GetColor(0, 0, 153));
  ZZTo->SetFillColor(TColor::GetColor(155, 204, 255));

  return {ZZTo, WZTo, TTTo, ggTo, DY};
}

namespace MCCRHistoIdx {
enum MCHistoIdx { ZZToIdx, WZToIdx, TTToIdx, ggToIdx, DYIdx };
}

template <bool isCR>
std::array<TH1D *, isCR ? N_CR_HISTOS : N_SR_HISTOS>
getMCHistos(TFile *inFile, Float_t lumi, const char *variable,
            FinalStateFlavor fsFlav, int binning, unitIdx unit) {
  LOG("[getMCHistos]\tGetting the histograms for the variable " << variable, 3);
  if constexpr (isCR) {
    return getMCHistosCR(inFile, lumi, variable, fsFlav, binning, unit);
  } else {
    return getMCHistosSR(inFile, lumi, variable, fsFlav, binning, unit);
  }
}

template <size_t N>
THStack *stackHistograms(const std::array<TH1D *, N> &histograms,
                         const char *variable, std::string xlabel,
                         FinalStateFlavor fsFlav, int binning, unitIdx unit) {

  TString stackName = makeCanvasName(variable, binning, unit, fsFlav);

  TString stackTitle;
  if (getUnitName(unit) == "")
    stackTitle = TString::Format(";%s;Events", xlabel.c_str());
  else
    stackTitle = TString::Format(";%s; Events / %i %s", xlabel.c_str(), binning,
                                 getUnitName(unit).c_str());

  LOG("[stackHistograms]\tCreating stack", 2);
  THStack *stack = new THStack(stackName, stackTitle);

  LOG("[stackHistograms]\tStacking histograms", 2);
  for (TH1D *histo : histograms) {
    if (histo->GetMinimum() < 0)
      LOGERR("[stackHistograms][WARNING]\tHistogram "
                 << histo->GetName() << " has negative minimum\n",
             1);
    stack->Add(histo, "HISTO");
  }

  return stack;
}

TGraphAsymmErrors *getDataGraph(TH1D *data, const bool isBlinded = true) {
  Int_t nBins = data->GetNbinsX();
  Double_t x, y, UpErr, LowErr, xErr(0);

  TGraphAsymmErrors *outGraph = new TGraphAsymmErrors();
  for (int i = 1; i <= data->GetNbinsX(); ++i) {
    x = data->GetBinCenter(i);
    bool inBlindRegion(false);
    for (blindInterval interval : blindIntervals) {
      inBlindRegion = inBlindRegion || interval.isIn(x);
    }
    if (isBlinded && inBlindRegion)
      continue;

    y = data->GetBinContent(i);

    if (y == 0.)
      continue;

    UpErr = data->GetBinErrorUp(i);
    LowErr = data->GetBinErrorLow(i);

    outGraph->AddPointError(x, y, 0, 0, LowErr, UpErr);
  }
  return outGraph;
}

TH1D *TGraphAsymmErrorsToTH1(TGraphAsymmErrors *graph, TH1D *model,
                             const std::string &name = "dataHisto") {
  Int_t nBins = model->GetNbinsX();
  TH1D *outHisto =
      new TH1D(name.c_str(), "", nBins, model->GetXaxis()->GetXmin(),
               model->GetXaxis()->GetXmax());

  for (size_t i = 0; i < graph->GetN(); i++) {
    Double_t x, y;
    graph->GetPoint(i, x, y);
    Int_t binIdx = outHisto->FindBin(x);
    outHisto->SetBinContent(binIdx, y);
    outHisto->SetBinError(
        binIdx, TMath::Max(graph->GetErrorYlow(i), graph->GetErrorYhigh(i)));
  }

  return outHisto;
}

template <bool isCR>
void addLegendEntries(
    TLegend *legend, TGraphAsymmErrors *data,
    const std::array<TH1D *, isCR ? N_CR_HISTOS : N_SR_HISTOS> &mcHistos,
    const char *variable) {
  legend->AddEntry(data, "Data", "p");
  if constexpr (isCR) {
    legend->AddEntry(mcHistos[MCCRHistoIdx::DYIdx], "DY", "f");
    legend->AddEntry(mcHistos[MCCRHistoIdx::ZZToIdx], "ZZ->4l", "f");
    legend->AddEntry(mcHistos[MCCRHistoIdx::TTToIdx], "TT->2l 2#nu", "f");
    legend->AddEntry(mcHistos[MCCRHistoIdx::WZToIdx], "WZ->3l 1#nu", "f");
    legend->AddEntry(mcHistos[MCCRHistoIdx::ggToIdx], "gg->4l", "f");
  } else {
    legend->AddEntry(mcHistos[MCSRHistoIdx::signalIdx], "H(125)", "f");
    legend->AddEntry(mcHistos[MCSRHistoIdx::ZZTo4lIdx],
                     "q#bar{q}#rightarrow ZZ", "f");
    legend->AddEntry(mcHistos[MCSRHistoIdx::ggToIdx], "gg#rightarrow ZZ", "f");
    legend->AddEntry(mcHistos[MCSRHistoIdx::EWIdx], "EW", "f");
    if (strcmp(variable, "ZZmass"))
      legend->AddEntry(mcHistos[MCSRHistoIdx::ZXIdx], "Z+X", "f");
  }
}

template <bool isCR>
void drawPlots(
    TCanvas *canvas, THStack *mcStack,
    const std::array<TH1D *, isCR ? N_CR_HISTOS : N_SR_HISTOS> &mcHistos,
    TGraphAsymmErrors *data, const char *variable, const bool blinded = true) {

  canvas->SetTicks();

  TPad *padStack = new TPad("stackPad", "stackPad", 0, 0.3, 1, 1);
  padStack->SetTopMargin(0.09);
  padStack->SetBottomMargin(0.01);
  padStack->Draw();

  TPad *padRatio = new TPad("ratioPad", "rationPad", 0, 0, 1, 0.3);
  padRatio->SetTopMargin(0.05);
  padRatio->SetBottomMargin(0.4);
  padRatio->Draw();

  Double_t xmin, xmax, ymin, ymax;
  data->ComputeRange(xmin, ymin, xmax, ymax);
  mcStack->SetMaximum(
      TMath::Ceil(1.1 * TMath::Max(mcStack->GetMaximum(), ymax)));

  padStack->cd();

  mcStack->Draw("HISTO");

  data->SetMarkerStyle(20);
  data->SetLineColor(kBlack);
  data->SetMarkerSize(0.9);

  data->Draw("SAMEPE1");

  TH1 *axisHist = mcStack->GetHistogram();
  axisHist->GetYaxis()->SetTitleSize(0.07);
  axisHist->GetYaxis()->SetLabelSize(0.05);
  axisHist->GetYaxis()->SetTitleOffset(0.8);

  mcStack->GetXaxis()->SetLabelSize(0);

  if (blinded && strcmp(variable, "ZZMass") == 0) {
    Float_t epsilon(0.1);
    for (blindInterval interval : blindIntervals) {
      TBox *bblind =
          new TBox(TMath::Max(interval.low, mcStack->GetXaxis()->GetXmin()), 0,
                   TMath::Min(interval.up, mcStack->GetXaxis()->GetXmax()),
                   mcStack->GetMaximum() - epsilon);
      bblind->SetFillColor(kGray);
      bblind->SetFillStyle(3002);
      bblind->Draw("SAME");
    }
  }

  gPad->RedrawAxis();

  TLegend *legend = new TLegend(0.65, 0.40, 0.95, 0.78);
  addLegendEntries<isCR>(legend, data, mcHistos, variable);
  legend->SetFillColor(kWhite);
  legend->SetLineColor(kWhite);
  legend->SetTextFont(43);
  legend->SetTextSize(32);

  legend->Draw();

  padRatio->cd();
  bool firstHisto(true);
  TH1D *hsum;
  for (TH1D *hist : mcHistos) {
    if (firstHisto) {
      hsum = (TH1D *)hist->Clone("hsum");
      firstHisto = false;
      continue;
    }
    hsum->Add(hist);
  }

  TH1D *dataHisto = TGraphAsymmErrorsToTH1(data, hsum);
  TH1D *ratio = (TH1D *)dataHisto->Clone("Ratio");
  delete dataHisto;

  ratio->Divide(hsum);

  ratio->SetMarkerStyle(20);
  ratio->SetTitle("");

  ratio->GetYaxis()->SetTitle("Data / MC");
  ratio->GetYaxis()->SetNdivisions(505);
  ratio->GetYaxis()->SetTitleSize(0.15);
  ratio->GetYaxis()->SetLabelSize(0.12);
  ratio->GetYaxis()->SetTitleOffset(0.35);

  ratio->GetXaxis()->SetTitleSize(0.15);
  ratio->GetXaxis()->SetLabelSize(0.12);
  ratio->GetXaxis()->SetTitle(mcStack->GetXaxis()->GetTitle());
  ratio->SetMinimum(0.5);
  ratio->SetMaximum(1.5);
  ratio->SetStats(false);
  ratio->Draw("PE1");

  TLine *line = new TLine(ratio->GetXaxis()->GetXmin(), 1.,
                          ratio->GetXaxis()->GetXmax(), 1.);
  line->SetLineColor(kBlack);
  line->SetLineStyle(2);
  line->SetLineWidth(2);
  line->Draw("SAME");

  canvas->Update();
}

void assertPeriod(const std::string &period) {
  assert(LumiPerPeriod.find(period) != LumiPerPeriod.end() &&
         "Luminosity not configured for current period");
  assert(MC_FILL_FILES.find(period) != MC_FILL_FILES.end() &&
         "H4lFill_MC file path not configured for current period");
  assert(DATA_FILL_FILES.find(period) != DATA_FILL_FILES.end() &&
         "H4lFill_Data fill path not configured for current period");
}

template <bool isMC> constexpr auto selectFillFiles() {
  if constexpr (isMC) {
    return MC_FILL_FILES;
  } else {
    return DATA_FILL_FILES;
  }
}

template <bool isMC>
void openFillFile(const std::string &period,
                  std::map<std::string, TFile *> &inFiles) {

  auto FILL_FILES = selectFillFiles<isMC>();

  if (inFiles.find(period) == inFiles.end()) {
    LOG("[openFillFile]\tOpening file " << FILL_FILES.at(period), 2);
    inFiles[period] = TFile::Open(FILL_FILES.at(period).c_str(), "READ");
    if (!inFiles[period] || inFiles[period]->IsZombie())
      throw std::runtime_error("Could not read file " +
                               FILL_FILES.at(period).string() + "\n");
  }
}

template <bool isCR> auto selectIdx() {
  if constexpr (isCR) {
    return std::array<MCCRHistoIdx::MCHistoIdx, 5>{
        MCCRHistoIdx::DYIdx, MCCRHistoIdx::TTToIdx, MCCRHistoIdx::WZToIdx,
        MCCRHistoIdx::ZZToIdx, MCCRHistoIdx::ggToIdx};
  } else {
    return std::array<MCSRHistoIdx::MCHistoIdx, 7>{
        MCSRHistoIdx::signalIdx, MCSRHistoIdx::EWIdx, MCSRHistoIdx::ZZTo4lIdx,
        MCSRHistoIdx::ggToIdx, /*MCSRHistoIdx::DYIdx, MCSRHistoIdx::TTIdx,*/
        MCSRHistoIdx::ZXIdx};
  }
}

std::string makeXlabel(const char *variable, bool isCR,
                       FinalStateFlavor fsFlav = FinalStateFlavor::L4,
                       CRIdx crIdx = P2F2) {
  std::string xlabel;
  LOG("[makeXlabel]\tCreating the xlabel of the canvas", 3);
  if (isCR && crIdx == ZL) {
    if (labelPerVariableZL.find(variable) == labelPerVariableZL.end()) {
      LOGERR("[makeXlabel]\tCould not find latex variable name for "
                 << variable << " using " << variable << " instead",
             0);
      return variable;
    }
    return labelPerVariableZL.at(variable);
  }

  if (labelPerVariable.find(variable) == labelPerVariable.end()) {
    LOGERR("[makeXlabel]\tCould not find latex variable name for "
               << variable << " using " << variable << " instead",
           0);
    return variable;
  }
  if (fsFlav != FinalStateFlavor::L4) {
    TString xlabelPreformat = labelPerVariable.at(variable);
    TString fsString =
        TString(makeFSString(fsFlav).substr(1)).ReplaceAll("mu", "#mu");
    return TString(xlabelPreformat.ReplaceAll("4l", fsString)).Data();
  }
  return labelPerVariable.at(variable);
}

template <typename T, size_t N>
std::array<T, N - 1> drop_elem(std::array<T, N> &arr, size_t Nth) {
  static_assert(N > 0, "Can't drop from empty array");
  assert(N > Nth);
  std::array<T, N - 1> out;
  std::copy(arr.begin(), arr.begin() + Nth, out.begin());
  std::copy(arr.begin() + Nth + 1, arr.end(), out.begin() + Nth);
  return out;
}

template <bool isCR, CRIdx crIdx = P2F2>
void writePlot(TFile *outFile, std::map<std::string, TFile *> &inMCFiles,
               std::map<std::string, TFile *> &inDataFiles,
               const std::vector<std::string> &periods, const char *variable,
               FinalStateFlavor fsFlav, int binning, unitIdx unit,
               bool isBlinded) {
  std::array<TH1D *, isCR ? N_CR_HISTOS : N_SR_HISTOS> MCHistograms;
  TH1D *dataHisto;
  std::string fullVariableName = getVariablePrefix(isCR, crIdx) + variable;

  TString canvasName =
      makeCanvasName(fullVariableName.c_str(), binning, unit, fsFlav);

  Int_t canvasXSize(1200), canvasYSize(900);

  std::string outDir;
  if constexpr (isCR) {
    outDir = "/CR" + CRNames.at(crIdx);
  } else {
    outDir = "/SR";
  }

  outFile->cd(outDir.c_str());

  TCanvas *canvas =
      new TCanvas(canvasName, canvasName, canvasXSize, canvasYSize);

  canvas->cd(); // To associate the histogram ptrs to the canvas

  bool firstPeriod(true);
  for (const std::string &period : periods) {

    float lumi = LumiPerPeriod.at(period);

    LOG("[writePlot]\tExtracting histograms for fullVariableName="
            << fullVariableName << " in the final state "
            << static_cast<Int_t>(fsFlav) << " at " << binning << " "
            << getUnitName(unit) << " binning",
        2);
    if (firstPeriod) {
      MCHistograms =
          getMCHistos<isCR>(inMCFiles[period], lumi, fullVariableName.c_str(),
                            fsFlav, binning, unit);

      dataHisto = agregateHistogram(
          inDataFiles[period], fullVariableName.c_str(), binning, unit,
          Data_dataSets2025,
          TString::Format("Data_histogram%s", period.c_str()).Data(), fsFlav);
      firstPeriod = !firstPeriod;
    } else {
      const auto currMCHistos = getMCHistos</*isCR=*/false>(
          inMCFiles[period], lumi, fullVariableName.c_str(), fsFlav, binning,
          unit);
      for (auto idx : selectIdx<isCR>())
        MCHistograms[idx]->Add(currMCHistos[idx]);

      dataHisto->Add(agregateHistogram(
          inDataFiles[period], fullVariableName.c_str(), binning, unit,
          Data_dataSets2025,
          TString::Format("Data_histogram%s", period.c_str()).Data(), fsFlav));
    }
  }

  LOG("[writePlot]\tStacking histograms", 1);
  THStack *stack;

  if (!isCR && fullVariableName == "ZZmass")
    stack = stackHistograms(MCHistograms, fullVariableName.c_str(),
                            makeXlabel(variable, isCR, fsFlav, crIdx), fsFlav,
                            binning, unit);
  else
    // Drop ZX
    stack = stackHistograms(
        drop_elem(MCHistograms, MCSRHistoIdx::ZXIdx), fullVariableName.c_str(),
        makeXlabel(variable, isCR, fsFlav, crIdx), fsFlav, binning, unit);

  LOG("[writePlot]\tTransforming data histograms to a graph", 1);
  TGraphAsymmErrors *data =
      getDataGraph(dataHisto, isBlinded && "ZZmass" == fullVariableName);

  LOG("[writePlot]\tLaying down the canvas", 1);
  drawPlots<isCR>(canvas, stack, MCHistograms, data, fullVariableName.c_str(),
                  isBlinded);

  outFile->cd(outDir.c_str());
  canvas->Write();
  delete canvas;
}

void H4l_Draw(const std::string &period = "2025", const bool &isBlinded = true,
              char verboseLevel = 0) {
  logging::verbose = verboseLevel;

  std::vector<std::string> periods;
  std::map<std::string, TFile *> inMCFiles, inDataFiles;

  if (dataSetPerPeriod.find(period) == dataSetPerPeriod.end()) {
    periods = {period};
    LOG("[H4l_Draw]\tPeriod "
            << period << " not found in the list, proceding with \"" << period
            << "\" "
               "as explicit period\n",
        0);
  } else {
    periods = dataSetPerPeriod.at(period);
    std::stringstream ss;
    ss << "Period used: " << period
       << "\n\t\tWill load the following data sets:";
    for (const auto &p : periods) {
      ss << "\n\t\t\t" << p;
    }
    LOG(ss.str(), 0);
  }

  LOG("[H4l_Draw]\tOpening files corresponding to " << period << " data set",
      1);
  for (const auto &period : periods) {
    assertPeriod(period);
    openFillFile</*isMC=*/true>(period, inMCFiles);
    openFillFile</*isMC=*/false>(period, inDataFiles);
  }

  setStyle();

  std::filesystem::path outPath = REPO_PATH / "outputs" / "H4lDraw.root";

  LOG("[H4l_Draw]\tOutputs will be written to " << outPath, 1);

  TFile *outFile = TFile::Open(outPath.c_str(), "RECREATE");

  outFile->cd(); // To associate the canvas ptrs here
  outFile->mkdir("SR");
  outFile->mkdir("CR2P2F");
  outFile->mkdir("CR3P1F");
  outFile->mkdir("CRSS");
  outFile->mkdir("CRZL");

  for (const auto &[unit, variables] : variablesPerUnit) {
    std::string unitName = getUnitName(unit);
    for (const auto &variable : variables) {
      for (HistogramBinning binning : binGranularites.at(unit))
        for (FinalStateFlavor fsFlav : fsFlavors) {
          writePlot</*isCR=*/false>(outFile, inMCFiles, inDataFiles, periods,
                                    variable.c_str(), fsFlav,
                                    int(binning.getBinSize()), unit, isBlinded);
          writePlot</*isCR=*/true, P2F2>(
              outFile, inMCFiles, inDataFiles, periods, variable.c_str(),
              fsFlav, int(binning.getBinSize()), unit, isBlinded);
          writePlot</*isCR=*/true, P3F1>(
              outFile, inMCFiles, inDataFiles, periods, variable.c_str(),
              fsFlav, int(binning.getBinSize()), unit, isBlinded);
          writePlot</*isCR=*/true, SS>(
              outFile, inMCFiles, inDataFiles, periods, variable.c_str(),
              fsFlav, int(binning.getBinSize()), unit, isBlinded);
        }
    }
  }

  for (const auto &[unit, variables] : variablesPerUnitZL)
    for (const auto &variable : variables) {
      std::string xlabel = labelPerVariableZL.at(variable);
      for (HistogramBinning binning : binGranularites.at(unit)) {
        writePlot</*isCR=*/true, ZL>(
            outFile, inMCFiles, inDataFiles, periods, variable.c_str(),
            FinalStateFlavor::L4, int(binning.getBinSize()), unit, isBlinded);
      }
    }

  LOG("[H4l_Draw]\tClosing file " << outPath, 0);
  outFile->cd();
  outFile->Close();
  return;
}
