#include <ROOT/RDF/HistoModels.hxx>
#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultPtr.hxx>
#include <ROOT/RVec.hxx>
#include <RtypesCore.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <cstdio>
#include <filesystem>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#pragma cling add_include_path(                                                \
    "/eos/user/m/mtokarev/work_area/currentCMSSW/src")

#include "PlotScripts/H4lPlots/src/ChannelNames.hxx"
#include "PlotScripts/H4lPlots/src/H4l_Common.C"
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
      printf("\tAdding file %s\n", file.path().c_str());
    }
  return files;
}

Float_t get_genEventSumwSample(const std::vector<std::string> &files) {
  Float_t genEventSumw(0);
  for (auto const &file : files) {
    Float_t curr;
    genEventSumw += get_genEventSumw(file);
    printf("\tFile %s's eventSumw=%f\n", file.c_str(), genEventSumw);
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

ROOT::RDF::RNode initializeVariables(ROOT::RDF::RNode df) {
  df = df.Define("ZZMass",
                 [](ROOT::RVecF ZZMass, Short_t bestIdx) -> Float_t {
                   return ZZMass[bestIdx];
                 },
                 {"ZZCand_mass", "bestCandIdx"})
           .Define("Z1Mass",
                   [](ROOT::RVecF Z1Mass, Short_t bestIdx) {
                     return Z1Mass[bestIdx];
                   },
                   {"ZZCand_Z1mass", "bestCandIdx"})
           .Define("Z2Mass",
                   [](ROOT::RVecF Z2Mass, Short_t bestIdx) {
                     return Z2Mass[bestIdx];
                   },
                   {"ZZCand_Z2mass", "bestCandIdx"})
           .Define("FSflav",
                   [](ROOT::RVecI Z1Flav, ROOT::RVecI Z2Flav, Short_t bestIdx) {
                     return Z1Flav[bestIdx] * Z2Flav[bestIdx];
                   },
                   {"ZZCand_Z1flav", "ZZCand_Z2flav", "bestCandIdx"});

  return df;
}

struct HistogramBinning {
  int GeVVal, nBins;
  float xlow, xhigh;
};

const HistogramBinning binGranularities[] = {
    HistogramBinning{.GeVVal = 2, .nBins = 65, .xlow = 70., .xhigh = 200.},
    HistogramBinning{.GeVVal = 4, .nBins = 233, .xlow = 70., .xhigh = 1002.},
    HistogramBinning{.GeVVal = 10, .nBins = 93, .xlow = 70., .xhigh = 1000.}};

const std::map<std::string, std::string> labelPerVariable = {
    {"ZZMass", "m_{#it{4l}}"}, {"Z1Mass", "m_{Z_1}"}, {"Z2Mass", "m_{Z_2}"}};

struct HistogramModel {
  ROOT::RDF::TH1DModel model;
  std::string variable;
  std::optional<TString> filter;
};

HistogramModel
makeModel(const std::string &variable, const std::string &sampleName,
          const std::string &xlabel, const HistogramBinning &binning,
          const std::optional<FinalStateFlavor> fsFlav = std::nullopt) {
  TString filterString;

  TString modelName = TString::Format("%s_%iGeV_%s", variable.c_str(),
                                      binning.GeVVal, sampleName.c_str());

  if (fsFlav.has_value()) {
    modelName += makeFSString(fsFlav.value());
  }

  HistogramModel histogram{
      .model = ROOT::RDF::TH1DModel(
          modelName,
          TString::Format("%s;%s;Events / %i GeV", modelName.Data(),
                          xlabel.c_str(), binning.GeVVal),
          binning.nBins, binning.xlow, binning.xhigh),
      .variable = variable,
  };

  if (fsFlav.has_value() && fsFlav.value() != FinalStateFlavor::L4)
    histogram.filter =
        TString::Format("FSflav==static_cast<Int_t>(%i)", fsFlav.value());

  return histogram;
}

std::vector<HistogramModel>
initializeHistogramModels(const std::string &sampleName) {
  std::vector<HistogramModel> models;

  for (const auto &[variable, xlabel] : labelPerVariable)
    for (HistogramBinning binning : binGranularities) {
      for (FinalStateFlavor fsFlav : fsFlavors) {
        auto currModel =
            makeModel(variable, sampleName, xlabel, binning, fsFlav);
        models.push_back(currModel);
      }
    }

  HistogramModel flavorPlot =
      HistogramModel{.model = ROOT::RDF::TH1DModel(
                         TString::Format("%s_%s", "FSflav", sampleName.c_str()),
                         TString::Format(";Event flavor pdgID"),
                         static_cast<Int_t>(FinalStateFlavor::Mu4) -
                             static_cast<Int_t>(FinalStateFlavor::E4) + 1,
                         static_cast<Int_t>(FinalStateFlavor::E4),
                         static_cast<Int_t>(FinalStateFlavor::Mu4)),
                     .variable = "FSflav"};

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

  df = initializeVariables(df);

  std::vector<ROOT::RDF::RResultPtr<TH1D>> histograms;

  std::vector<HistogramModel> histoModels =
      initializeHistogramModels(sample.name);

  for (auto model : histoModels) {
    auto selection = df;
    if (model.filter.has_value())
      selection = selection.Filter(model.filter.value().Data());

    auto hist = selection.Histo1D(model.model, model.variable, "eventWeight");
    histograms.push_back(hist);
  }

  outFile.cd();

  for (auto &hist : histograms)
    hist->Write();
}

template <bool isMC> constexpr auto selectDataSet() {
  if constexpr (isMC)
    return MC_dataSets;
  else
    return Data_dataSets;
}

template <bool isMC> void run(const std::string &period = "2025") {

  TString outName = TString::Format("H4lFill_%s%s.root", isMC ? "MC" : "Data",
                                    period.c_str());
  auto dataSets = selectDataSet<isMC>();
  auto INPUT_PATH = isMC ? MC_PATHS : DATA_PATHS;
  auto outFile = std::unique_ptr<TFile>(TFile::Open(outName, "RECREATE"));

  printf("Output will be saved to %s\n", outName.Data());

  for (std::string channel : dataSets) {
    DataSet sample =
        DataSet{.name = channel,
                .path = std::filesystem::path(INPUT_PATH[period] / channel)};

    printf("Filling histograms from ROOT file(s) in %s \n",
           sample.path.c_str());

    fillHistos<isMC>(sample, *outFile);
  }
  printf("Closing file %s\n", outName.Data());
  outFile->Close();
}

void runMC(const std::string &period = "2025") { run<true>(period); }

void runData(const std::string &period = "2025") { run<false>(period); }

void H4l_Fill(const std::vector<std::string> &years = {"2025"},
              int nEvents = -1) {
  for (std::string year : years) {
    runMC(year);
    runData(year);
  }
  return;
}
