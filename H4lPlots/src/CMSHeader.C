#include <TLatex.h>
#include <TPad.h>
#include <iomanip>
#include <ios>
#include <sstream>
#include <string>

class CMS {
  std::string type = "";
  float luminosity, cmEnergy;

public:
  CMS(std::string type, float luminosity, float cmEnergy)
      : type(type), luminosity(luminosity), cmEnergy(cmEnergy) {};
  ~CMS() = default;

  void formatPad(TPad *pad) {
    pad->cd();
    double left = pad->GetLeftMargin();
    double right = 1 - pad->GetRightMargin();
    double topMargin = pad->GetTopMargin();
    double top = 1 - topMargin;
    double y = top + 0.2 * topMargin;

    TLatex CMSHeader;
    CMSHeader.SetNDC();

    // CMS Private/Preliminary/...
    CMSHeader.SetTextFont(42);
    CMSHeader.SetTextSize(0.75 * pad->GetTopMargin());

    std::ostringstream oss;
    oss << "#bf{CMS} #it{" << type << "}";

    CMSHeader.DrawLatex(left, y, oss.str().c_str());
    oss.str("");
    oss.clear();

    // lumi fb-1 (sqrtS TeV)
    oss << std::fixed << std::setprecision(1) << luminosity << " fb^{-1} "
        << "(" << std::fixed << std::setprecision(1) << cmEnergy << " TeV)";
    CMSHeader.SetTextAlign(31);
    CMSHeader.SetTextSize(0.05);
    CMSHeader.DrawLatex(right, y, oss.str().c_str());
  }
};
