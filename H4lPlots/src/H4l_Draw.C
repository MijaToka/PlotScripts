#pragma cling add_include_path(                                                \
    "/eos/user/m/mtokarev/work_area/currentCMSSW/src")

#include "PlotScripts/H4lPlots/interface/ChannelNames.hxx"

#include <THStack.h>
#include <TString.h>
#include <cstddef>
#include <vector>

const int energyBinning[] = {2, 4, 10};

std::vector<THStack> createMCStack() {

  std::vector<THStack> stacks;

  for (size_t i = 0; i < sizeof(energyBinning) / sizeof(int); i++) {
    int E = energyBinning[i];
    TString stackName = TString::Format("ZZmass_%iGeV_MC", E);
    TString stackTitle =
        TString::Format("ZZmass_%iGeV_MC;m_{#it{4l}};Events / %i GeV", E, E);
    stacks.emplace_back(stackName, stackTitle);

    for (auto &sampleName : MC_dataSets) {
      TString histoName =
          TString::Format("ZZmass_%iGev_%s", E, sampleName.c_str());
    }
  }
  stacks.add(histos);
}

void H4l_Stack() {

  std::vector<THStack> stacks = createMCStack();

  runData();
}
