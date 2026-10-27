
#include <string>

#pragma cling add_include_path(                                                \
    "/eos/user/m/mtokarev/work_area/currentCMSSW/src")
#include "PlotScripts/H4lPlots/src/ChannelNames.hxx"

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
