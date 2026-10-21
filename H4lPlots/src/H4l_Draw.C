#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TString.h>
#include <filesystem>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

#pragma cling add_include_path(                                                \
    "/eos/user/m/mtokarev/work_area/currentCMSSW/src")

#include "PlotScripts/H4lPlots/interface/ChannelNames.hxx"
#include "PlotScripts/H4lPlots/interface/Paths.hxx"

const int energyBinning[] = {2, 4, 10};

enum stacksIndex { per2GeV, per4GeV, per10GeV };
const stacksIndex stackIndices[] = {per2GeV, per4GeV, per10GeV};

template <bool isMC = true> struct HistoAggregate {
  std::vector<std::unique_ptr<std::conditional<isMC, THStack, TH1D>>>
      aggregates;
  std::vector<std::unique_ptr<TH1D>> histograms;
};

template <bool isMC> constexpr auto selectDataSets() {
  if constexpr (isMC)
    return MC_dataSets;
  else
    return Data_dataSets;
}

template <bool isMC> HistoAggregate<isMC> agregateHistograms() {

  std::filesystem::path inputFileName =
      isMC ? "H4lPlots_MC2025.root" : "H4lPlots_Data2025.root";

  std::filesystem::path inputPath(REPO_PATH / "output" / inputFileName);

  TFile *inMCFile = TFile::Open(inputPath.c_str(), "READ");

  if (!inMCFile || inMCFile->IsZombie()) {
    std::stringstream ss;
    ss << "File " << inputPath << " is invalid.\n";
    throw std::runtime_error(ss.str());
  }

  HistoAggregate<isMC> out;

  for (int E : energyBinning) {
    TString stackName = TString::Format("ZZmass_%iGeV_MC", E);
    TString stackTitle =
        TString::Format("ZZmass_%iGeV_%s;m_{#it{4l}};Events / %i GeV", E,
                        isMC ? "MC" : "Data", E);

    auto stack = std::make_unique<THStack>(stackName, stackTitle);

    for (auto &sampleName : selectDataSets<isMC>()) {
      TString histoName =
          TString::Format("ZZmass_%iGev_%s", E, sampleName.c_str());
      TH1D *histo = (TH1D *)inMCFile->Get(histoName);

      if (!histo) {
        throw std::runtime_error(
            ("Histogram " + histoName + " not found in " + inputPath.c_str())
                .Data());
      }

      // detatch from TFile
      histo->SetDirectory(nullptr);

      // Take ownership so it survives the closing of the TFile
      out.histograms.emplace_back(histo);

      out.aggregates.back()->Add(histo);
    }
    out.aggregates.push_back(std::move(stack));
  }
  return out;
}

void H4l_Stack() { return; }
