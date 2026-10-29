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
#include <TMath.h>
#include <TMathBase.h>
#include <TPad.h>
#include <TROOT.h>
#include <TString.h>
#include <TStyle.h>
#include <TVirtualPad.h>
#include <array>
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#pragma cling add_include_path(                                                \
    "/eos/user/m/mtokarev/work_area/currentCMSSW/src")

#include "PlotScripts/H4lPlots/src/ChannelNames.hxx"
#include "PlotScripts/H4lPlots/src/H4l_Common.C"
#include "PlotScripts/H4lPlots/src/Paths.hxx"
#include "PlotScripts/H4lPlots/src/Style.hxx"

const std::map<std::string, float> LumiPerPeriod = {
    {"2025", 109.}}; // Place-holder value

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
TH1D *agregateHistogram(
    TFile *file, const char *variable, int energyBinning,
    const std::array<const char *, N> &samples,
    const std::string &sampleBundleName,
    const std::optional<FinalStateFlavor> &fsFlav = std::nullopt) {

  std::string fsString("");

  if (fsFlav.has_value()) {
    fsString = makeFSString(fsFlav.value());
  }

  TH1D *histogram = nullptr;

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
    return TString::Format("%s_%iGeV_%s%s%s", variable, energyBinning,
                           sampleName.c_str(), makeFSString(fsFlav).c_str(),
                           period)
        .Data();
  };

  auto mkHisto = [inFile, variable, energyBinning, fsFlav, mkBundleName](
                     const std::string &sampleName, const auto &dataSet) {
    return agregateHistogram(inFile, variable, energyBinning /*GeV*/, dataSet,
                             mkBundleName(sampleName), fsFlav);
  };

  TH1D *signal = mkHisto("signal", MC_SignalH125);
  signal->Scale(lumi * 1000.);

  TH1D *EW = mkHisto("EW", MC_BkgEW);
  EW->Scale(lumi * 1000.);

  TH1D *ZZTo4l = mkHisto("qqZZ", MC_BkgQQZZ);
  ZZTo4l->Scale(lumi * 1000.);

  TH1D *ggTo = mkHisto("ggZZ", MC_BkgGGTo);
  ggTo->Scale(lumi * 1000.);

  TH1D *DY = mkHisto("DY", MC_BkgDY);
  DY->Scale(lumi * 1000.);

  TH1D *TT = mkHisto("ttTo", MC_BkgTT);
  TT->Scale(lumi * 1000.);

  TH1D *ZX = getZX(ZZTo4l, fsFlav);
  ZX->Scale(lumi * 1000.);

  return {signal, EW, ZZTo4l, ggTo, DY, TT, ZX};
}

THStack *stackHistograms(const std::array<TH1D *, 7> &histograms,
                         const char *variable, FinalStateFlavor fsFlav,
                         int energyBinning) {
  auto &[signal, EW, ZZTo4l, ggToZZ, DY, TT, ZX] = histograms;
  std::string fsString = makeFSString(fsFlav);

  TString stackName =
      TString::Format("%s_%iGeV%s", variable, energyBinning, fsString.c_str());

  std::string latexVariable;
  auto it = latexVariableNames.find((variable + makeFSString(fsFlav)));
  if (it == latexVariableNames.end()) {
    printf("Could not find latex variable for %s\n",
           (variable + makeFSString(fsFlav)).c_str());
    latexVariable = variable;
  } else
    latexVariable =
        latexVariableNames.at((variable + makeFSString(fsFlav)).c_str());

  TString stackTitle = TString::Format(";%s; Events / %i GeV",
                                       latexVariable.c_str(), energyBinning);
  THStack *stack = new THStack(stackName, stackTitle);

  signal->SetLineColor(TColor::GetColor(204, 00, 00));
  signal->SetFillColor(TColor::GetColor(255, 155, 155));

  EW->SetLineColor(TColor::GetColor(0, 0, 153));
  EW->SetFillColor(TColor::GetColor(3, 49, 185));

  ZZTo4l->SetLineColor(TColor::GetColor(0, 0, 153));
  ZZTo4l->SetFillColor(TColor::GetColor(155, 204, 255));

  ggToZZ->SetLineColor(TColor::GetColor(0, 0, 153));
  ggToZZ->SetFillColor(TColor::GetColor(75, 120, 255));

  DY->SetLineColor(TColor::GetColor(0, 0, 153));
  DY->SetFillColor(kOrange);

  ZX->SetLineColor(TColor::GetColor(0, 51, 0));
  ZX->SetFillColor(TColor::GetColor(102, 153, 102));

  stack->Add(ZX, "HISTO");
  stack->Add(EW, "HISTO");
  stack->Add(ggToZZ, "HISTO");
  stack->Add(ZZTo4l, "HISTO");
  // stack->Add(DY, "HISTO");
  // stack->Add(TT, "HISTO");

  for (TH1D *histo : histograms)
    if (histo->GetMinimum() < 0)
      printf("Histogram %s has negative minimum\n", histo->GetName());

  stack->Add(signal, "HISTO");

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

void drawPlots(TCanvas *canvas, THStack *mcStack,
               const std::array<TH1D *, 7> &mcHistos, TGraphAsymmErrors *data,
               const bool blinded = true) {

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
  data->Draw("SAMEPE1");

  TH1 *axisHist = mcStack->GetHistogram();
  axisHist->GetYaxis()->SetTitleSize(0.07);
  axisHist->GetYaxis()->SetLabelSize(0.05);
  axisHist->GetYaxis()->SetTitleOffset(0.8);

  mcStack->GetXaxis()->SetLabelSize(0);

  if (blinded) {
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
  legend->AddEntry(data, "Data", "p");
  legend->AddEntry(mcHistos[signalIdx], "H(125)", "f");
  legend->AddEntry(mcHistos[ZZTo4lIdx], "q#bar{q}#rightarrow ZZ", "f");
  legend->AddEntry(mcHistos[ggToIdx], "gg#rightarrow ZZ", "f");
  legend->AddEntry(mcHistos[EWIdx], "EW", "f");
  legend->AddEntry(mcHistos[ZXIdx], "Z+X", "f");

  legend->SetFillColor(kWhite);
  legend->SetLineColor(kWhite);
  legend->SetTextFont(43);
  legend->SetTextSize(32);

  legend->Draw();

  padRatio->cd();

  TH1D *hsum = (TH1D *)mcHistos[signalIdx]->Clone("hsum");
  for (MCHistoIdx idx : {EWIdx, ZZTo4lIdx, ggToIdx, /*DYIdx, TTIdx,*/ ZXIdx})
    hsum->Add(mcHistos[idx]);

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

  canvas->Update();
}

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
  std::map<std::string, TFile *> inMCFiles, inDataFiles;

  const char *variable = "ZZMass";
  int energyBinning(4);
  float lumi(109);

  outFile->cd(); // To associate the canvas ptrs here

  for (FinalStateFlavor fsFlav : fsFlavors) {
    for (const auto &energyBinning : {2, 4, 10}) {
      bool firstFile(true);
      std::array<TH1D *, 7> MCHistograms;
      TH1D *dataHisto;

      TString canvasName = TString::Format(
          "%s_%iGeV%s", variable, energyBinning, makeFSString(fsFlav).c_str());
      Int_t canvasXSize(1200), canvasYSize(900);

      TCanvas *canvas =
          new TCanvas(canvasName, canvasName, canvasXSize, canvasYSize);

      canvas->cd(); // To associate the histogram ptrs to the canvas

      for (const std::string &period : periods) {

        if (inMCFiles.find(period) == inMCFiles.end()) {
          inMCFiles[period] =
              TFile::Open(MC_FILL_FILES.at(period).c_str(), "READ");
          if (!inMCFiles[period] || inMCFiles[period]->IsZombie())
            throw std::runtime_error("Could not read file " +
                                     MC_FILL_FILES.at(period).string() + "\n");
        }

        if (inDataFiles.find(period) == inDataFiles.end()) {
          inDataFiles[period] =
              TFile::Open(DATA_FILL_FILES.at(period).c_str(), "READ");
          if (!inDataFiles[period] || inDataFiles[period]->IsZombie())
            throw std::runtime_error("Could not read file " +
                                     DATA_FILL_FILES.at(period).string() +
                                     "\n");
        }

        if (firstFile) {
          MCHistograms = getMCHistos(inMCFiles[period], lumi, variable, fsFlav,
                                     energyBinning);
          dataHisto = agregateHistogram(
              inDataFiles[period], variable, energyBinning, Data_dataSets,
              TString::Format("Data_histogram%s", period.c_str()).Data(),
              fsFlav);
          firstFile = !firstFile;
        } else {
          const auto currMCHistos = getMCHistos(
              inMCFiles[period], lumi, variable, fsFlav, energyBinning);
          for (MCHistoIdx idx :
               {signalIdx, EWIdx, ZZTo4lIdx, ggToIdx, DYIdx, TTIdx, ZXIdx})
            MCHistograms[idx]->Add(currMCHistos[idx]);

          dataHisto->Add(agregateHistogram(
              inDataFiles[period], variable, energyBinning, Data_dataSets,
              TString::Format("Data_histogram%s", period.c_str()).Data(),
              fsFlav));
        }
      }

      THStack *stack =
          stackHistograms(MCHistograms, variable, fsFlav, energyBinning);

      TGraphAsymmErrors *data = getDataGraph(dataHisto, isBlinded);

      drawPlots(canvas, stack, MCHistograms, data, isBlinded);

      outFile->cd();
      canvas->Write();
    }
  }
  outFile->Close();
  return;
}
