#include <ROOT/RDF/HistoModels.hxx>
#include <ROOT/RDF/InterfaceUtils.hxx>
#include <Rtypes.h>
#include <RtypesCore.h>
#include <TCanvas.h>
#include <TColor.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <array>
#include <cstddef>
#include <filesystem>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#pragma cling add_include_path(                                                \
    "/eos/user/m/mtokarev/work_area/currentCMSSW/src")

#include "PlotScripts/H4lPlots/interface/ChannelNames.hxx"
#include "PlotScripts/H4lPlots/interface/Paths.hxx"

const std::map<std::string, float> LumiPerPeriod = {
    {"2025", 109.}}; // Place-holder value

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

inline std::string makeFSString(FinalStateFlavor fsFlav) {

  std::string fsString("");
  switch (fsFlav) {
  case FinalStateFlavor::E4:
    fsString = "_4e";
    break;
  case FinalStateFlavor::E2Mu2:
    fsString = "_2e2mu";
    break;
  case FinalStateFlavor::Mu4:
    fsString = "_4mu";
    break;
  case FinalStateFlavor::L4:
  default:
    break;
  }
  return fsString;
}

template <size_t N>
TH1D *agregateHistogram(
    TFile *file, const char *variable, int energyBinning,
    const std::array<const char *, N> &samples,
    const std::string &sampleBundleName,
    const std::optional<FinalStateFlavor> &fsFlav = std::nullopt) {

  std::string fsString("");

  if (fsFlav.has_value()) {
    fsString = makeFSString(fsFlav.value());
  }

  TH1D *histogram;

  for (auto sample : samples) {
    TString inHistName = TString::Format(
        "%s_%iGeV_%s%s", variable, energyBinning, sample, fsString.c_str());
    TH1D *currHist = (TH1D *)file->Get(inHistName);

    if (!currHist)
      throw std::runtime_error("Histogram not found: " + inHistName +
                               " in file " + file->GetPath());

    if (!histogram) {
      histogram = (TH1D *)currHist->Clone(sampleBundleName.c_str());
    } else {
      histogram->Add(currHist);
    }
  }
  return histogram;
}

enum MCHistoIdx { signalIdx, EWIdx, ZZTo4lIdx, ggToIdx, DYIdx, TTIdx, ZXIdx };
std::array<TH1D *, 7> getMCHistos(TFile *inFile, Float_t lumi,
                                  const char *variable, FinalStateFlavor fsFlav,
                                  int energyBinning, const char *period = "") {

  std::string sampleBundleName;
  auto mkBundleName = [variable, fsFlav, energyBinning, period](
                          const std::string &sampleName) -> const std::string {
    return TString::Format("%s_%s%s%s", variable, sampleName.c_str(),
                           makeFSString(fsFlav).c_str(), period)
        .Data();
  };

  sampleBundleName = mkBundleName("signal");
  TH1D *signal = agregateHistogram(inFile, variable, energyBinning /*GeV*/,
                                   MC_SignalH125, sampleBundleName);
  signal->Scale(lumi * 1000.);

  sampleBundleName = mkBundleName("EW");
  TH1D *EW = agregateHistogram(inFile, variable, energyBinning /*GeV*/,
                               MC_BkgEW, sampleBundleName);
  EW->Scale(lumi * 1000.);

  sampleBundleName = mkBundleName("qqZZ");
  TH1D *ZZTo4l = agregateHistogram(inFile, variable, energyBinning /*GeV*/,
                                   MC_BkgQQZZ, sampleBundleName);
  ZZTo4l->Scale(lumi * 1000.);

  sampleBundleName = mkBundleName("ggZZ");
  TH1D *ggTo = agregateHistogram(inFile, variable, energyBinning /*GeV*/,
                                 MC_BkgGGTo, sampleBundleName);
  ggTo->Scale(lumi * 1000.);

  sampleBundleName = mkBundleName("DY");
  TH1D *DY = agregateHistogram(inFile, variable, energyBinning /*GeV*/,
                               MC_BkgDY, sampleBundleName);
  DY->Scale(lumi * 1000.);

  sampleBundleName = mkBundleName("ttTo");
  TH1D *TT = agregateHistogram(inFile, variable, energyBinning /*GeV*/,
                               MC_BkgTT, sampleBundleName);
  TT->Scale(lumi * 1000.);

  TH1D *ZX = getZX(signal, fsFlav);
  ZX->Scale(lumi * 1000.);

  return {signal, EW, ZZTo4l, ggTo, DY, TT, ZX};
}
const std::map<const char *, std::string> latexVariableNames = {
    {"ZZMass", "m_{#it{4l}}"},
    {"Z1Mass", "m_{Z_{1}}"},
    {"Z2Mass", "m_{Z_{2}}"},
};

THStack stackHistograms(const std::array<TH1D *, 7> &histograms,
                        const char *variable, FinalStateFlavor fsFlav,
                        int energyBinning) {
  auto &[signal, EW, ZZTo4l, ggToZZ, DY, TT, ZX] = histograms;
  std::string fsString = makeFSString(fsFlav);

  TString stackName =
      TString::Format("%s_%iGeV%s", variable, energyBinning, fsString.c_str());

  std::string latexVariable;
  auto it = latexVariableNames.find(variable);
  if (it == latexVariableNames.end())
    latexVariable = variable;
  else
    latexVariable = latexVariableNames.at(variable);

  TString stackTitle = TString::Format(";%s; Events / %i GeV",
                                       latexVariable.c_str(), energyBinning);
  THStack stack = THStack(stackName, stackTitle);

  Int_t lineColor;
  signal->SetLineColor(TColor::GetColor(204, 00, 00));
  signal->SetFillColor(TColor::GetColor(255, 155, 155));

  EW->SetFillColor(TColor::GetColor(lineColor, 0, 0, 153));
  EW->SetFillColor(TColor::GetColor(3, 49, 185));

  ZZTo4l->SetLineColor(lineColor);
  ZZTo4l->SetFillColor(TColor::GetColor(155, 204, 255));

  ggToZZ->SetLineColor(lineColor);
  ggToZZ->SetFillColor(TColor::GetColor(75, 120, 255));

  DY->SetLineColor(lineColor);
  DY->SetFillColor(kOrange);

  ZX->SetLineColor(TColor::GetColor(0, 51, 0));
  ZX->SetFillColor(TColor::GetColor(102, 153, 102));

  stack.Add(ZX, "HISTO");
  stack.Add(EW, "HISTO");
  stack.Add(ggToZZ, "HISTO");
  stack.Add(ZZTo4l, "HISTO");
  // stack.Add(DY, "HISTO");
  // stack.Add(TT, "HISTO");

  stack.Add(signal, "HISTO");

  return stack;
}

void setStyle() {
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetLabelOffset(0.008, "XYZ");
  gStyle->SetLabelSize(0.04, "XYZ");
  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(true);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetTitleOffset(1.00, "X");
  gStyle->SetTitleOffset(1.25, "Y");
  gStyle->SetLabelOffset(0.008, "XYZ");
  gStyle->SetLabelSize(0.04, "XYZ");
}

const std::unordered_map<std::string, std::vector<std::string>>
    dataSetPerPeriod = {{"Run3",
                         {"2022", "2022EE", "2023preBPix", "2023postBPix",
                          "2024", "2025", "2026"}}};

void H4l_Draw(const std::string &period = "2025",
              const bool &isBlinded = true) {

  std::vector<std::string> periods;

  if (dataSetPerPeriod.find(period) == dataSetPerPeriod.end())
    periods = {period};
  else
    periods = dataSetPerPeriod.at(period);

  setStyle();

  std::filesystem::path outPath = REPO_PATH / "outputs" / "H4lDraw.root";

  TFile *outFile = TFile::Open(outPath.c_str(), "RECREATE");
  bool firstFile(true);
  std::map<std::string, TFile *> inFiles;

  std::array<TH1D *, 7> MCHistograms;

  const char *variable = "ZZMass";
  int energyBinning(4);
  float lumi(100);

  for (FinalStateFlavor fsFlav : fsFlavors) {
    for (const std::string &period : periods) {

      if (inFiles.find(period) == inFiles.end()) {
        inFiles[period] = TFile::Open(MC_FILL_FILES.at(period).c_str(), "READ");
        if (!inFiles[period] || inFiles[period]->IsZombie())
          throw std::runtime_error("Could not read file " +
                                   MC_FILL_FILES.at(period).string() + "\n");
      }

      if (firstFile) {
        MCHistograms =
            getMCHistos(inFiles[period], lumi, variable, fsFlav, energyBinning);
      } else {
        const auto currMCHistos =
            getMCHistos(inFiles[period], lumi, variable, fsFlav, energyBinning);
        for (MCHistoIdx idx :
             {signalIdx, EWIdx, ZZTo4lIdx, ggToIdx, DYIdx, TTIdx, ZXIdx}) {
          MCHistograms[idx]->Add(currMCHistos[idx]);
        }
      }
    }

    THStack stack =
        stackHistograms(MCHistograms, variable, fsFlav, energyBinning);
    TCanvas canvas = TCanvas();
    outFile->cd();
    stack.Write();
  }

  outFile->Close();
  return;
}
