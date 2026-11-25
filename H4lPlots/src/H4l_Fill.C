#include <ROOT/RDF/HistoModels.hxx>
#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RResultHandle.hxx>
#include <ROOT/RResultPtr.hxx>
#include <ROOT/RVec.hxx>
#include <Rtypes.h>
#include <RtypesCore.h>
#include <TFile.h>
#include <TH1.h>
#include <TH1F.h>
#include <TMath.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <algorithm>
#include <array>
#include <cstring>
#include <filesystem>
#include <string>
#include <sys/types.h>
#include <vector>

#pragma cling add_include_path(                                                \
    "/eos/user/m/mtokarev/work_area/currentCMSSW/src")

#include "PlotScripts/H4lPlots/interface/Logger.hxx"
#include "PlotScripts/H4lPlots/src/ChannelNames.hxx"
#include "PlotScripts/H4lPlots/src/HistogramModel.hxx"
#include "PlotScripts/H4lPlots/src/Paths.hxx"
#include "PlotScripts/H4lPlots/src/Variables.C"

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

template <bool hasCR>
ROOT::RDF::RNode initializeRDF(const std::vector<std::string> &files) {
  ROOT::RDF::RNode df = ROOT::RDataFrame("Events", files);
  df = df.Define("eventWeight", []() { return 1.; }, {});

  if constexpr (hasCR) {
    df = df.Define("eventWeightCR2P2F", []() { return 1.; }, {})
             .Define("eventWeightCR3P1F", []() { return 1.; }, {})
             .Define("eventWeightCRSS", []() { return 1.; }, {})
             .Define("eventWeightCRZL", []() { return 1.; }, {});
  }
  return df;
}

template <bool hasCR>
ROOT::RDF::RNode calculateEventWeight(Float_t genEventSumw,
                                      ROOT::RDF::RNode df) {

  auto calculateWeight = [genEventSumw](Float_t overallWeight,
                                        ROOT::RVecF dataMCWeight,
                                        Short_t bestIdx) -> Float_t {
    if (genEventSumw == 0 || bestIdx == -1)
      return 0.;
    else
      return overallWeight * dataMCWeight[bestIdx] / genEventSumw;
  };

  df =
      df.Redefine("eventWeight", calculateWeight,
                  {"overallEventWeight", "ZZCand_dataMCWeight", "bestCandIdx"});

  if constexpr (hasCR) {
    df = df.Redefine(
               "eventWeightCR2P2F", calculateWeight,
               {"overallEventWeight", "ZLLCand_dataMCWeight", "ZLLbest2P2FIdx"})
             .Redefine("eventWeightCR3P1F", calculateWeight,
                       {"overallEventWeight", "ZLLCand_dataMCWeight",
                        "ZLLbest3P1FIdx"})
             .Redefine(
                 "eventWeightCRSS", calculateWeight,
                 {"overallEventWeight", "ZLLCand_dataMCWeight", "ZLLbestSSIdx"})
             .Redefine(
                 "eventWeightCRZL",
                 [genEventSumw](
                     Float_t overallWeight, ROOT::RVecF ElectronWeights,
                     ROOT::RVecF MuonWeights, int nElectron, Short_t bestZIdx,
                     Short_t bestZLIdx, ROOT::VecOps::RVec<Short_t> l1Idxs,
                     ROOT::VecOps::RVec<Short_t> l2Idxs) -> Float_t {
                   if (genEventSumw == 0 || bestZLIdx == -1) {
                     return 0;
                   } else {
                     Float_t dataMCWeight(1);
                     for (Short_t lIdx : {l1Idxs[bestZIdx], l2Idxs[bestZIdx]}) {
                       if (lIdx < 0)
                         continue;

                       if (lIdx < nElectron) {
                         dataMCWeight *= ElectronWeights[lIdx];
                       } else {
                         dataMCWeight *= MuonWeights[lIdx - nElectron];
                       }
                     }

                     return overallWeight * dataMCWeight / genEventSumw;
                   }
                 },
                 {"overallEventWeight", "Electron_dataMC", "Muon_dataMC",
                  "nElectron", "bestZIdx", "ZLCand_lepIdx", "ZCand_l1Idx",
                  "ZCand_l2Idx"});
  }

  return df;
}

template <bool hasCR>
std::vector<HistogramModel>
initializeHistogramModels(const std::string &sampleName) {
  std::vector<HistogramModel> models;

  // Z+ll plots
  for (const auto &[unit, variables] : variablesPerUnit)
    for (const auto &variable : variables) {
      std::string xlabel = labelPerVariable.at(variable);
      for (HistogramBinning binning : binGranularites.at(unit)) {
        for (FinalStateFlavor fsFlav : fsFlavors) {
          auto currModel = HistogramModel(/*isCR=*/false, variable, sampleName,
                                          xlabel, binning, unit, fsFlav);
          LOG("[initializeHistogramModels]\t Created model: "
                  << currModel.getModel().fTitle,
              3);
          models.push_back(currModel);

          // Z+X control region plots
          if constexpr (hasCR) {
            for (CRIdx crIdx : {P2F2, P3F1, SS}) {
              auto currModel =
                  HistogramModel(/*isCR=*/true, variable, sampleName, xlabel,
                                 binning, unit, fsFlav, crIdx);

              LOG("[initializeHistogramModels]\t Created model: "
                      << currModel.getModel().fTitle << " for the variable "
                      << currModel.getVariable(),
                  3);
              models.push_back(currModel);
            }
          }
        }
      }
    }

  // Z+l control region plots
  if constexpr (hasCR) {
    for (const auto &[unit, variables] : variablesPerUnitZL)
      for (const auto &variable : variables) {
        std::string xlabel = labelPerVariableZL.at(variable);
        for (HistogramBinning binning : binGranularites.at(unit)) {
          for (FinalStateFlavor fsFlav : fsFlavorsZL) {
            auto currModel = HistogramModel(
                /*isCR=*/true, variable, sampleName, xlabel, binning, unit,
                fsFlav, /*crIdx=*/ZL);
            LOG("[initializeHistogramModels]\t Created model: "
                    << currModel.getModel().fTitle << " for the variable "
                    << currModel.getVariable(),
                3);
            models.push_back(currModel);
          }
        }
      }
  }

  // Final state flavor plots as a sanity check
  HistogramModel flavorPlot = HistogramModel(
      /*isCR=*/false,
      /*variable=*/"FSflav", sampleName, /*xlabel=*/"Event flavor pdgID",
      /*binning=*/
      HistogramBinning{.nBins = static_cast<Int_t>(FinalStateFlavor::Mu4) -
                                static_cast<Int_t>(FinalStateFlavor::E4) + 1,
                       .xlow = static_cast<Int_t>(FinalStateFlavor::E4),
                       .xhigh = static_cast<Int_t>(FinalStateFlavor::Mu4)},
      /*unit=*/Idx, FinalStateFlavor::L4);

  models.push_back(flavorPlot);
  if constexpr (hasCR) {
    for (CRIdx crIdx : {P2F2, P3F1, SS}) {
      HistogramModel flavorPlot = HistogramModel(
          /*isCR=*/true,
          /*variable=*/"FSflav", sampleName, /*xlabel=*/"Event flavor pdgID",
          /*binning=*/
          HistogramBinning{
              .nBins = 2 * static_cast<Int_t>(FinalStateFlavor::Mu4) + 1,
              .xlow = -static_cast<Int_t>(FinalStateFlavor::Mu4),
              .xhigh = static_cast<Int_t>(FinalStateFlavor::Mu4)},
          /*unit=*/Idx, FinalStateFlavor::L4, crIdx);

      models.push_back(flavorPlot);
    }
  }
  return models;
}

template <bool isMC, bool hasCR>
void fillHistos(const DataSet &sample, TFile &outFile) {
  std::vector<std::string> files = getFileNames(sample);
  ROOT::RDF::RNode df = initializeRDF<hasCR>(files);

  if constexpr (isMC) {
    Float_t genEventSumw = get_genEventSumwSample(files);
    df = calculateEventWeight<hasCR>(genEventSumw, df);
  }
  LOG("[fillHistos]\tInitializing variables", 1);
  df = initializeVariables<hasCR>(df);

  std::vector<ROOT::RDF::RResultPtr<TH1D>> histograms;
  std::vector<ROOT::RDF::RResultHandle> sumWs;

  std::vector<ROOT::RDF::RResultHandle> handles;

  LOG("[fillHistos]\tCreating TH1D models", 1);
  std::vector<HistogramModel> histoModels =
      initializeHistogramModels<hasCR>(sample.name);

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
            << model.getVariable() << " with weightColmn " << model.getWeight(),
        2);
    auto hist = selection.Histo1D(model.getModel(), model.getVariable(),
                                  model.getWeight());
    histograms.push_back(hist);
    handles.push_back(hist);

    auto sumW = selection.Sum(model.getWeight());
    sumWs.push_back(sumW);
    handles.push_back(sumW);
  }

  LOG("[fillHistos]\tRunning event loop for " << sample.name, 1);
  ROOT::RDF::RunGraphs(handles);

  outFile.cd();

  LOG("[fillHistos]\tWriting histograms to output file", 1);

  for (size_t i = 0; i < histograms.size(); ++i) {
    LOG("[fillHistos]\tWriting hitogram "
            << histograms[i]->GetName() << " with "
            << sumWs[i].GetValue<Double_t>() << " events",
        2);
    histograms[i]->Write();
  }
}

template <bool isMC> constexpr auto selectDataSet(const std::string &period) {
  if constexpr (isMC)
    return MC_dataSets;
  else {
    assert(periodIndex.find(period) != periodIndex.end() &&
           "Data period not configured");
    return Data_dataSets[periodIndex.at(period)];
  }
}

template <bool isMC> void run(const std::string &period = "2025") {

  TString outName = TString::Format("H4lFill_%s%s.root", isMC ? "MC" : "Data",
                                    period.c_str());
  View<const char *> dataSets = selectDataSet<isMC>(period);
  auto INPUT_PATH = isMC ? MC_PATHS : DATA_PATHS;

  auto outFile = std::unique_ptr<TFile>(TFile::Open(outName, "RECREATE"));

  LOG("[run]\tOutput will be saved to " << outName, 0);

  for (std::string channel : dataSets) {
    DataSet sample =
        DataSet{.name = channel,
                .path = std::filesystem::path(INPUT_PATH[period] / channel)};

    bool hasCR =
        !isMC || std::find_if(MC_ControlRegion.begin(), MC_ControlRegion.end(),
                              [sample](const char *s) {
                                return std::strcmp(s, sample.name.c_str()) == 0;
                              }) != MC_ControlRegion.end();
    LOG("[run]\tFilling histograms from ROOT file(s) in "
            << sample.path << " with isMC: " << isMC << " hasCR: " << hasCR,
        0);

    if (hasCR)
      fillHistos<isMC, true>(sample, *outFile);
    else
      fillHistos<isMC, false>(sample, *outFile);
  }
  LOG("[run]\tClosing file " << outName, 0);
  outFile->Close();
}

void runMC(const std::string &period = "2025") { run<true>(period); }

void runData(const std::string &period = "2025") { run<false>(period); }

void H4l_Fill(const std::vector<std::string> &years = {"2025"},
              char verboseLevel = 0) {
  logging::verbose = verboseLevel;
  ROOT::EnableThreadSafety();
  ROOT::EnableImplicitMT();
  for (std::string year : years) {
    runMC(year);
    runData(year);
  }
  return;
}
