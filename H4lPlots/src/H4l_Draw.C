#include <ROOT/RDF/InterfaceUtils.hxx>
#include <TFile.h>
#include <TH1.h>
#include <THStack.h>
#include <TString.h>
#include <array>
#include <cstddef>
#include <map>
#include <optional>
#include <shared_mutex>
#include <stdexcept>
#include <string>

#pragma cling add_include_path(                                                \
    "/eos/user/m/mtokarev/work_area/currentCMSSW/src")

#include "PlotScripts/H4lPlots/interface/ChannelNames.hxx"
#include "PlotScripts/H4lPlots/interface/Paths.hxx"

const std::map<std::string, float> LumiPerPeriod = {
    {"2025", 109.}}; // Place-holder value

template <size_t N>
TH1D *agregateHistogram(
    TFile *file, const std::string &variable, const int &energyBinning,
    const std::array<char *, N> &samples, const std::string &sampleBundleName,
    const std::optional<FinalStateFlavor> &fsFlav = std::nullopt) {

  std::string fsString("");

  if (fsFlav.has_value()) {
    switch (fsFlav.value()) {
    case FinalStateFlavor::E4:
      fsString = "_4e";
      break;
    case FinalStateFlavor::E2Mu2:
      fsString = "_2e2mu";
      break;
    case FinalStateFlavor::Mu4:
      fsString = "_4mu";
      break;
    }
  }

  TH1D *histogram;

  for (auto sample : samples) {
    TString inHistName =
        TString::Format("%s_%iGeV_%s%s", variable.c_str(), energyBinning,
                        sample, fsString.c_str());
    TH1D *currHist = (TH1D *)file->Get(inHistName);

    if (!currHist)
      throw std::runtime_error("Histogram not found: " + inHistName);

    if (!histogram) {
      histogram = (TH1D *)currHist->Clone(sampleBundleName.c_str());
    } else {
      histogram->Add(currHist);
    }
  }
  return histogram;
}

template <bool isMC> constexpr auto selectDataSets(const char *selection) {
  if constexpr (isMC)
    // Here selection is Hmass
    return MC_dataSetsH125;
  else
    return Data_dataSets;
}

void H4l_Draw(const std::string &period = "2025",
              const bool &isBlinded = true) {

  return;
}
