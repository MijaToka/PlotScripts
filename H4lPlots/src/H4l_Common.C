#ifndef H4L_COMMON
#define H4L_COMMON

#include <TString.h>
#include <string>

#pragma cling add_include_path(                                                \
    "/eos/user/m/mtokarev/work_area/currentCMSSW/src")
#include "PlotScripts/H4lPlots/src/ChannelNames.hxx"
#include "PlotScripts/H4lPlots/src/Variables.C"

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
