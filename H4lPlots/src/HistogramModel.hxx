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
  std::string variable, xlabel, sampleName, filter, name;
  unitIdx unit;
  HistogramBinning binning;
  FinalStateFlavor fsFlav = FinalStateFlavor::L4;
  bool hasFilter_ = false;

  void update() {

    TString binSizeLabel =
        binning.getBinSize() < 1
            ? TString::Format("%02d", int(binning.getBinSize() * 100))
            : TString::Format("%i", int(binning.getBinSize()));

    std::string unitName(getUnitNameForTitle(unit));
    if (unitName == "") {
      name = TString::Format("%s_%s", variable.c_str(), sampleName.c_str());
    } else {
      name =
          TString::Format("%s_%s%s_%s", variable.c_str(), binSizeLabel.Data(),
                          unitName.c_str(), sampleName.c_str());
    }

    if (fsFlav != FinalStateFlavor::L4) {
      filter = TString::Format("FSflav==%i", static_cast<Int_t>(fsFlav)).Data();
      hasFilter_ = true;
      std::string fsString = makeFSString(fsFlav);
      name += fsString;
      TString fsReplace = TString(fsString.substr(1)).ReplaceAll("mu", "#mu");
      xlabel = TString(xlabel).ReplaceAll("4l", fsReplace);
      std::string unitName = getUnitName(unit);
      xlabel += unitName == "" ? "" : " [" + unitName + "]";
    }

    TString title;
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
  HistogramModel(ROOT::RDF::TH1DModel model, std::string variable)
      : model(model), variable(variable) {};
  HistogramModel(ROOT::RDF::TH1DModel model, std::string variable,
                 TString filter)
      : model(model), variable(variable), filter(filter), hasFilter_(true) {};

  HistogramModel(std::string variable, std::string sampleName,
                 std::string xlabel, HistogramBinning binning, unitIdx unit,
                 FinalStateFlavor fsFlav = FinalStateFlavor::L4)
      : variable(variable), xlabel(xlabel), sampleName(sampleName),
        binning(binning), fsFlav(fsFlav), unit(unit) {
    update();
  }

  HistogramBinning getBinning() { return binning; }
  std::string_view getFilter() { return filter; }
  ROOT::RDF::TH1DModel getModel() { return model; }
  std::string getVariable() { return variable; }

  bool hasFilter() { return hasFilter_; }

  void setBinning(HistogramBinning newBinning) { binning = newBinning; }
  void setFilter(TString newFilter) {
    hasFilter_ = true;
    filter = newFilter;
  }
};
