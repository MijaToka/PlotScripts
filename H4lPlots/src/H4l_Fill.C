#include <ROOT/RDF/HistoModels.hxx>
#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <ROOT/RVec.hxx>
#include <RtypesCore.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <filesystem>
#include <map>
#include <memory>
#include <string>
#include <vector>

#pragma cling add_include_path(                                                \
    "/eos/user/m/mtokarev/work_area/currentCMSSW/src")

#include "PlotScripts/H4lPlots/interface/Logger.hxx"
#include "PlotScripts/H4lPlots/src/ChannelNames.hxx"
#include "PlotScripts/H4lPlots/src/HistogramModel.hxx"
#include "PlotScripts/H4lPlots/src/Paths.hxx"

Double_t get_genEventSumw(const std::string &file_path) {
  TFile *file = TFile::Open(file_path.c_str(), "READ");

  TTree *runs = (TTree *)file->Get("Runs");
  // TTree *events = (TTree *)file->Get("Events");

  Double_t genEventSumw(0), eventSumw;
  Long64_t genEventCount(0), eventCount;

  runs->SetBranchAddress("genEventSumw", &eventSumw);
  runs->SetBranchAddress("genEventCount", &eventCount);

  for (int iRun = 0; iRun < runs->GetEntries(); iRun++) {
    runs->GetEntry(iRun);
    genEventSumw += eventSumw;
    genEventCount += eventCount;
  }

  file->Close();
  return genEventSumw;
};

struct DataSet {
  std::string name;
  std::filesystem::path path;
};

std::vector<std::string> getFileNames(const DataSet &sample) {
  std::vector<std::string> files;
  for (auto const &file : std::filesystem::directory_iterator{sample.path})
    if (file.path().extension() == ".root" && file.is_regular_file()) {
      files.push_back(file.path().string());
      LOG("[getFileNames]\tAdding file " << file.path(), 1);
    }
  return files;
}

Float_t get_genEventSumwSample(const std::vector<std::string> &files) {
  Float_t genEventSumw(0);
  for (auto const &file : files) {
    Float_t curr;
    genEventSumw += get_genEventSumw(file);
    LOG("[get_genEventSumwSample]\tFile " << file
                                          << "'s eventSumw=" << genEventSumw,
        2);
  }
  return genEventSumw;
}

ROOT::RDF::RNode initializeRDF(const std::vector<std::string> &files) {
  ROOT::RDF::RNode df = ROOT::RDataFrame("Events", files);
  df = df.Filter(
             [](Short_t bestIdx, Bool_t HLT_passZZ4l) {
               return bestIdx >= 0 && HLT_passZZ4l;
             },
             {"bestCandIdx", "HLT_passZZ4l"})
           .Define("eventWeight", []() { return 1.; }, {});
  return df;
}

ROOT::RDF::RNode calculateEventWeight(Float_t genEventSumw,
                                      ROOT::RDF::RNode df) {
  df = df.Redefine(
      "eventWeight",
      [genEventSumw](Float_t overallWeight, ROOT::RVecF dataMCWeight,
                     Short_t bestIdx) -> Float_t {
        if (genEventSumw == 0 || bestIdx == -1)
          return 0.;
        else
          return overallWeight * dataMCWeight[bestIdx] / genEventSumw;
      },
      {"overallEventWeight", "ZZCand_dataMCWeight", "bestCandIdx"});

  return df;
}

#define MK_FLATTEN_VARIABLE(columnTag, containerClass)                         \
  "ZZ" columnTag,                                                              \
      [](const containerClass &ZZCand, Short_t bestIdx) {                      \
        return ZZCand[bestIdx];                                                \
      },                                                                       \
      {"ZZCand_" columnTag, "bestCandIdx"}

ROOT::RDF::RNode initializeVariables(ROOT::RDF::RNode df) {
  df = df.Define(MK_FLATTEN_VARIABLE("mass", ROOT::RVecF))
           .Define("Z1mass",
                   [](ROOT::RVecF Z1Mass, Short_t bestIdx) {
                     return Z1Mass[bestIdx];
                   },
                   {"ZZCand_Z1mass", "bestCandIdx"})
           .Define("Z2mass",
                   [](ROOT::RVecF Z2Mass, Short_t bestIdx) {
                     return Z2Mass[bestIdx];
                   },
                   {"ZZCand_Z2mass", "bestCandIdx"})
           .Define("FSflav",
                   [](ROOT::RVecI Z1Flav, ROOT::RVecI Z2Flav, Short_t bestIdx) {
                     return Z1Flav[bestIdx] * Z2Flav[bestIdx];
                   },
                   {"ZZCand_Z1flav", "ZZCand_Z2flav", "bestCandIdx"})
           .Define(MK_FLATTEN_VARIABLE("pt", ROOT::RVecF))
           .Define(MK_FLATTEN_VARIABLE("eta", ROOT::RVecF))
           .Define(MK_FLATTEN_VARIABLE("phi", ROOT::RVecF))
           .Define(MK_FLATTEN_VARIABLE("costheta1", ROOT::RVecF))
           .Define(MK_FLATTEN_VARIABLE("costheta2", ROOT::RVecF))
           .Define(MK_FLATTEN_VARIABLE("costhetastar", ROOT::RVecF))
           .Define(MK_FLATTEN_VARIABLE("Phi", ROOT::RVecF))
           .Define(MK_FLATTEN_VARIABLE("Phi1", ROOT::RVecF));

  return df;
}

const std::vector<HistogramBinning> massBinGranularities = {
    HistogramBinning{.nBins = 65, .xlow = 70., .xhigh = 200.},
    HistogramBinning{.nBins = 233, .xlow = 70., .xhigh = 1002.},
    HistogramBinning{.nBins = 93, .xlow = 70., .xhigh = 1000.},
};
const std::vector<HistogramBinning> ptBinGranularities = {
    HistogramBinning{.nBins = 100, .xlow = 0., .xhigh = 500},
};
const std::vector<HistogramBinning> cosBinGranularities = {
    HistogramBinning{.nBins = 30, .xlow = -1, .xhigh = 1},
};
const std::vector<HistogramBinning> etaBinGranularities = {
    HistogramBinning{.nBins = 100, .xlow = -4, .xhigh = 4},
};
const std::vector<HistogramBinning> angleBinGranularities = {
    HistogramBinning{.nBins = 100, .xlow = -TMath::Pi(), .xhigh = TMath::Pi()},
};

const std::map<std::string, std::string> labelPerVariable = {
    {"ZZmass", "m_{#it{4l}}"},
    {"Z1mass", "m_{Z_1}"},
    {"Z2mass", "m_{Z_2}"},
    {"ZZpt", "p_{T}^{#it{4l}}"},
    {"ZZeta", "#eta_{#it{4l}}"},
    {"ZZphi", "#phi_{#it{4l}}"},
    {"ZZcostheta1", "cos#theta_{1}"},
    {"ZZcostheta2", "cos#theta_{2}"},
    {"ZZcosthetastar", "cos#theta^{*}"},
    {"ZZPhi", "#Phi_{#it{4l}}"},
    {"ZZPhi1", "#Phi1_{#it{4l}}"},
};

enum unitIdx { GeV, GeV_c, angle, cosAngle, eta };
const std::map<unitIdx, std::vector<HistogramBinning>> binGranularites = {
    {GeV, massBinGranularities},    {GeV_c, ptBinGranularities},
    {angle, angleBinGranularities}, {cosAngle, cosBinGranularities},
    {eta, etaBinGranularities},
};
const std::map<unitIdx, std::vector<std::string>> variablesPerUnit = {
    {GeV, {"ZZmass", "Z1mass", "Z2mass"}},
    {GeV_c, {"ZZpt"}},
    {angle, {"ZZphi", "ZZPhi", "ZZPhi1"}},
    {cosAngle, {"ZZcostheta1", "ZZcostheta2", "ZZcosthetastar"}},
    {eta, {"ZZeta"}},
};

std::string getUnitName(unitIdx unit) {
  switch (unit) {
  case GeV:
    return "GeV";
  case GeV_c:
    return "GeV/c";
  case angle:
  case cosAngle:
  case eta:
  default:
    return "";
  }
}

std::vector<HistogramModel>
initializeHistogramModels(const std::string &sampleName) {
  std::vector<HistogramModel> models;

  for (const auto &[unit, variables] : variablesPerUnit) {
    std::string unitName = getUnitName(unit);
    for (const auto &variable : variables) {
      std::string xlabel = labelPerVariable.at(variable);
      for (HistogramBinning binning : binGranularites.at(unit)) {
        for (FinalStateFlavor fsFlav : fsFlavors) {
          auto currModel = HistogramModel(variable, sampleName, xlabel, binning,
                                          fsFlav, unitName);
          LOG("[initializeHistogramModels]\t Created model: "
                  << currModel.getModel().fTitle << " for the variable "
                  << currModel.getVariable(),
              3);
          models.push_back(currModel);
        }
      }
    }
  }

  HistogramModel flavorPlot = HistogramModel(
      /*variable=*/"FSflav", sampleName, /*xlabel=*/"Event flavor pdgID",
      /*binning=*/
      HistogramBinning{.nBins = static_cast<Int_t>(FinalStateFlavor::Mu4) -
                                static_cast<Int_t>(FinalStateFlavor::E4) + 1,
                       .xlow = static_cast<Int_t>(FinalStateFlavor::E4),
                       .xhigh = static_cast<Int_t>(FinalStateFlavor::Mu4)},
      FinalStateFlavor::L4, /*unitName=*/"");

  //  HistogramModel flavorPlot =
  //      HistogramModel{/*.model=*/
  //                     ROOT::RDF::TH1DModel(
  //                         TString::Format("%s_%s", "FSflav",
  //                         sampleName.c_str()), TString::Format(";Event flavor
  //                         pdgID"), static_cast<Int_t>(FinalStateFlavor::Mu4)
  //                         -
  //                             static_cast<Int_t>(FinalStateFlavor::E4) + 1,
  //                         static_cast<Int_t>(FinalStateFlavor::E4),
  //                         static_cast<Int_t>(FinalStateFlavor::Mu4)),
  //                     /*.variable=*/"FSflav"};

  models.push_back(flavorPlot);

  return models;
}

template <bool isMC> void fillHistos(const DataSet &sample, TFile &outFile) {
  std::vector<std::string> files = getFileNames(sample);
  ROOT::RDF::RNode df = initializeRDF(files);

  if constexpr (isMC) {
    Float_t genEventSumw = get_genEventSumwSample(files);
    df = calculateEventWeight(genEventSumw, df);
  }
  LOG("[fillHistos]\tInitializing variables", 1);
  df = initializeVariables(df);

  std::vector<ROOT::RDF::RResultPtr<TH1D>> histograms;

  LOG("[fillHistos]\tCreating TH1D models", 1);
  std::vector<HistogramModel> histoModels =
      initializeHistogramModels(sample.name);

  for (auto model : histoModels) {
    auto selection = df;
    LOG("[fillHistos]\tChecking if model " << model.getModel().fTitle
                                           << " has filter",
        2);
    if (model.hasFilter()) {
      LOG("[fillHistos]\tApplying filter: " << model.getFilter(), 2);
      selection = selection.Filter(model.getFilter());
    }
    LOG("[fillHistos]\tCreating histogram "
            << model.getModel().fName << " with title "
            << model.getModel().fTitle << " for the variable "
            << model.getVariable(),
        2);
    auto hist =
        selection.Histo1D(model.getModel(), model.getVariable(), "eventWeight");
    histograms.push_back(hist);
  }

  outFile.cd();

  LOG("[fillHistos]\tWriting histograms to output file", 1);
  for (auto &hist : histograms)
    hist->Write();
}

template <bool isMC> constexpr auto selectDataSet() {
  if constexpr (isMC)
    return MC_dataSets;
  else
    return Data_dataSets2025;
}

template <bool isMC> void run(const std::string &period = "2025") {

  TString outName = TString::Format("H4lFill_%s%s.root", isMC ? "MC" : "Data",
                                    period.c_str());
  auto dataSets = selectDataSet<isMC>();
  auto INPUT_PATH = isMC ? MC_PATHS : DATA_PATHS;

  auto outFile = std::unique_ptr<TFile>(TFile::Open(outName, "RECREATE"));

  LOG("[run]\tOutput will be saved to " << outName, 0);

  for (std::string channel : dataSets) {
    DataSet sample =
        DataSet{.name = channel,
                .path = std::filesystem::path(INPUT_PATH[period] / channel)};

    LOG("[run]\tFilling histograms from ROOT file(s) in " << sample.path, 0);

    fillHistos<isMC>(sample, *outFile);
  }
  LOG("[run]\tClosing file " << outName, 0);
  outFile->Close();
}

void runMC(const std::string &period = "2025") { run<true>(period); }

void runData(const std::string &period = "2025") { run<false>(period); }

void H4l_Fill(const std::vector<std::string> &years = {"2025"},
              int nEvents = -1) {
  logging::verbose = 1;
  for (std::string year : years) {
    runMC(year);
    runData(year);
  }
  return;
}
