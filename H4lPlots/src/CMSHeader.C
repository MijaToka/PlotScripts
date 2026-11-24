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
    TLatex CMSHeader;
    CMSHeader.SetNDC();
    CMSHeader.SetTextSize(0.07);
    std::ostringstream oss;
    oss << "#font[61]{CMS} #font[52]{" << type << "}";
    CMSHeader.DrawLatex(0.12, 0.92, oss.str().c_str());
    oss.clear();
    oss << "#font[42]{" << std::fixed << std::setprecision(1) << luminosity
        << " fb^{-1} (" << std::fixed << std::setprecision(1) << cmEnergy
        << " TeV})";
    CMSHeader.DrawLatex(0.9, 0.92, oss.str().c_str());
  }
};
