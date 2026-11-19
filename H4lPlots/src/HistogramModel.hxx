#include <ROOT/RDF/HistoModels.hxx>
#include <TString.h>
#include <string>
#include <string_view>

#pragma cling add_include_path(                                                \
    "/eos/user/m/mtokarev/work_area/currentCMSSW/src")
#include "PlotScripts/H4lPlots/src/ChannelNames.hxx"
#include "PlotScripts/H4lPlots/src/H4l_Common.C"
#include "PlotScripts/H4lPlots/src/Variables.C"

class HistogramModel {
  ROOT::RDF::TH1DModel model;
  std::string variable, weightColumn, xlabel, sampleName, filter, name;
  unitIdx unit;
  HistogramBinning binning;
  FinalStateFlavor fsFlav = FinalStateFlavor::L4;
  bool hasFilter_ = false;
  bool isCR_;
  CRIdx crIdx;

  void update() {

    auto mkName = [this](const std::string &variablePrefix) {
      return makeHistoName((variablePrefix + variable).c_str(),
                           int(binning.getBinSize()), unit, sampleName.c_str(),
                           fsFlav);
    };

    weightColumn = "eventWeight";
    name = mkName(getVariablePrefix(isCR_, crIdx));
    if (isCR_) {
      if (crIdx != ZL) {
        addFilter(
            TString::Format("ZLLbest%sIdx>=0", CRNames.at(crIdx).c_str()));
      } else {
        addFilter(TString("ZLCand_lepIdx>=0"));
      }
      weightColumn += "CR" + CRNames.at(crIdx);
    } else {
      addFilter("bestCandIdx>=0");
    }

    TString binSizeLabel =
        binning.getBinSize() < 1
            ? TString::Format("%02d", int(binning.getBinSize() * 100))
            : TString::Format("%i", int(binning.getBinSize()));

    if (fsFlav != FinalStateFlavor::L4) {
      TString flavorFilter("true");
      if (isCR_) {
        if (crIdx != ZL)
          flavorFilter =
              TString::Format("abs(%s%sFSflav)==%i", isCR_ ? "ZLL" : "ZZ",
                              isCR_ ? CRNames.at(crIdx).c_str() : "",
                              static_cast<Int_t>(fsFlav));
      } else {
        flavorFilter =
            TString::Format("abs(ZZFSflav)==%i", static_cast<Int_t>(fsFlav));
      }
      addFilter(flavorFilter);

      TString fsReplace =
          TString(makeFSString(fsFlav).substr(1)).ReplaceAll("mu", "#mu");
      xlabel = TString(xlabel).ReplaceAll("4l", fsReplace);
      std::string unitName = getUnitName(unit);
      xlabel += (unitName == "") ? "" : " [" + unitName + "]";
    }

    TString title;
    std::string unitName = getUnitName(unit);
    if (unitName == "")
      title = TString::Format("%s;%s;Events", name.c_str(), xlabel.c_str());
    else
      title =
          TString::Format("%s;%s;Events / %i %s", name.c_str(), xlabel.c_str(),
                          int(binning.getBinSize()), unitName.c_str());

    model = ROOT::RDF::TH1DModel(name.c_str(), title, binning.nBins,
                                 binning.xlow, binning.xhigh);
  }

public:
  HistogramModel(bool isCR, std::string variable, std::string sampleName,
                 std::string xlabel, HistogramBinning binning, unitIdx unit,
                 FinalStateFlavor fsFlav = FinalStateFlavor::L4,
                 CRIdx crIdx = P2F2)
      : variable(variable), xlabel(xlabel), sampleName(sampleName),
        binning(binning), fsFlav(fsFlav), unit(unit), isCR_(isCR),
        crIdx(crIdx) {
    update();
  }

  void addFilter(TString newFilter) {
    if (!hasFilter_) {
      setFilter(newFilter);
    } else {
      filter += "&&" + newFilter;
    }
  }

  HistogramBinning getBinning() { return binning; }
  std::string_view getFilter() { return filter; }
  ROOT::RDF::TH1DModel getModel() { return model; }
  std::string getVariable() {
    return getVariablePrefix(isCR_, crIdx) + variable;
  }
  std::string getWeight() { return weightColumn; }

  bool hasFilter() { return hasFilter_; }
  bool isCR() { return isCR_; }

  void setBinning(HistogramBinning newBinning) { binning = newBinning; }
  void setFilter(TString newFilter) {
    hasFilter_ = true;
    filter = newFilter;
  }
};
