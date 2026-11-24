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
  bool hasFilter_ = false;
  bool weightOverride = false;
  bool isCR_;
  unitIdx unit;
  CRIdx crIdx;
  FinalStateFlavor fsFlav = FinalStateFlavor::L4;
  HistogramBinning binning;
  std::string variable, xlabel, sampleName, filter;
  std::string filterNameSuffix = "";
  std::string weightColumn = "eventWeight";

  void addCandidateFilter() {
    if (isCR_) {
      if (crIdx != ZL) {
        addFilter(
            TString::Format("ZLLbest%sIdx>=0", CRNames.at(crIdx).c_str()));
      } else {
        addFilter(TString("ZLCand_lepIdx>=0"));
      }
    } else {
      addFilter("bestCandIdx>=0");
    }
  }
  void applyFSFlavorFilter() {
    if (fsFlav == FinalStateFlavor::L4)
      return;
    if (!isCR_ || crIdx == ZL) {
      addFilter(TString::Format("abs(%sFSflav)==%i", isCR_ ? "ZL" : "ZZ",
                                static_cast<Int_t>(fsFlav)));
      return;
    }
    if (crIdx == ZL) {
      addFilter(
          TString::Format("abs(ZLFSflav)==%i", static_cast<Int_t>(fsFlav)));
      return;
    }
    addFilter(TString::Format("abs(ZLL%sFSflav)==%i", CRNames.at(crIdx).c_str(),
                              static_cast<Int_t>(fsFlav)));
    return;
  }
  void initializeXlabel() {
    if (fsFlav != FinalStateFlavor::L4) {
      TString fsReplace =
          TString(makeFSString(fsFlav).substr(1)).ReplaceAll("mu", "#mu");
      xlabel = TString(xlabel).ReplaceAll("4l", fsReplace);
      std::string unitName = getUnitName(unit);
      xlabel += (unitName == "") ? "" : " [" + unitName + "]";
    }
  }

  void initialize() {
    addCandidateFilter();
    applyFSFlavorFilter();
    initializeXlabel();
  }

public:
  HistogramModel(bool isCR, std::string variable, std::string sampleName,
                 std::string xlabel, HistogramBinning binning, unitIdx unit,
                 FinalStateFlavor fsFlav = FinalStateFlavor::L4,
                 CRIdx crIdx = P2F2)
      : variable(variable), xlabel(xlabel), sampleName(sampleName),
        binning(binning), fsFlav(fsFlav), unit(unit), isCR_(isCR),
        crIdx(crIdx) {
    initialize();
  }

  std::string getName() {
    auto mkName = [this](const std::string &variablePrefix) {
      return makeHistoName((variablePrefix + variable).c_str(),
                           int(binning.getBinSize()), unit, sampleName.c_str(),
                           fsFlav);
    };
    return mkName(getVariablePrefix(isCR_, crIdx)).Data() + filterNameSuffix;
  }

  std::string getXlabel() { return xlabel; }

  void addFilter(TString newFilter) {
    if (!hasFilter_) {
      setFilter(newFilter);
    } else {
      filter += "&&" + newFilter;
    }
  }
  void addFilter(TString filter, std::string filterSuffix) {
    addFilter(filter);
    filterNameSuffix += filterSuffix;
  }

  void resetFilter() {
    filter = "";
    hasFilter_ = false;
  }

  HistogramBinning getBinning() { return binning; }

  std::string_view getFilter() { return filter; }

  ROOT::RDF::TH1DModel getModel() {
    std::string name = getName();
    TString title;
    std::string unitName = getUnitName(unit);
    if (unitName == "")
      title = TString::Format("%s;%s;Events", name.c_str(), xlabel.c_str());
    else
      title =
          TString::Format("%s;%s;Events / %i %s", name.c_str(), xlabel.c_str(),
                          int(binning.getBinSize()), unitName.c_str());
    return ROOT::RDF::TH1DModel(name.c_str(), title, binning.nBins,
                                binning.xlow, binning.xhigh);
  }

  std::string getVariable() {
    return getVariablePrefix(isCR_, crIdx) + variable;
  }

  std::string getWeight() {
    if (weightOverride)
      return weightColumn;

    std::string weightColumnSuffix = "";
    if (isCR_) {
      weightColumnSuffix += "CR" + CRNames.at(crIdx);
      if (crIdx == ZL && variable.rfind("lepton", 0) != std::string::npos)
        weightColumnSuffix += "lepton";
    }
    return weightColumn + weightColumnSuffix;
  }

  bool hasFilter() { return hasFilter_; }
  bool isCR() { return isCR_; }

  void setBinning(HistogramBinning newBinning) { binning = newBinning; }
  void setFilter(TString newFilter) {
    hasFilter_ = true;
    filter = newFilter;
  }

  void setWeight(std::string weight) {
    weightColumn = weight;
    weightOverride = true;
  }

  void resetWeight() {
    weightColumn = "eventWeight";
    weightOverride = false;
  }
};
