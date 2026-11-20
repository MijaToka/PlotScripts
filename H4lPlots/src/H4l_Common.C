#ifndef H4L_COMMON
#define H4L_COMMON

#include <TString.h>
#include <string>

#pragma cling add_include_path(                                                \
    "/eos/user/m/mtokarev/work_area/currentCMSSW/src")
#include "PlotScripts/H4lPlots/src/ChannelNames.hxx"
#include "PlotScripts/H4lPlots/src/Variables.C"

inline std::string makeFSString(FinalStateFlavor fsFlav) {
  switch (fsFlav) {
  case FinalStateFlavor::E4:
    return "_4e";
  case FinalStateFlavor::E2Mu2:
    return "_2e2mu";
  case FinalStateFlavor::Mu4:
    return "_4mu";

  // Z+l
  case FinalStateFlavor::E3:
    return "_3e";
  case FinalStateFlavor::E2Mu1:
    return "_2e1mu";
  case FinalStateFlavor::Mu2E1:
    return "_2mu1e";
  case FinalStateFlavor::Mu3:
    return "_3mu";
  // Inclusive
  case FinalStateFlavor::L4:
    return "";
  }
}

inline TString makeHistoName(const char *variable, int energyBinning,
                             unitIdx unit, const char *sample,
                             FinalStateFlavor fsFlav) {
  if (getUnitNameForTitle(unit) == "") {
    return TString::Format("%s_%s%s", variable, sample,
                           makeFSString(fsFlav).c_str());
  } else {
    return TString::Format("%s_%i%s_%s%s", variable, energyBinning,
                           getUnitNameForTitle(unit).c_str(), sample,
                           makeFSString(fsFlav).c_str());
  }
}

inline TString makeCanvasName(const char *variable, int energyBinning,
                              unitIdx unit, FinalStateFlavor fsFlav) {
  if (getUnitNameForTitle(unit) == "") {
    return TString::Format("%s%s", variable, makeFSString(fsFlav).c_str());
  } else {
    return TString::Format("%s_%i%s%s", variable, energyBinning,
                           getUnitNameForTitle(unit).c_str(),
                           makeFSString(fsFlav).c_str());
  }
}

#endif
