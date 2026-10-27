#include <TH1.h>
#include <TStyle.h>

inline void setStyle() {
  TH1::SetDefaultSumw2();
  gStyle->SetErrorX(0);
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadLeftMargin(0.16);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetLabelOffset(0.008, "XYZ");
  gStyle->SetLabelSize(0.04, "XYZ");
  gStyle->SetAxisColor(1, "XYZ");
  gStyle->SetStripDecimals(true);
  gStyle->SetTickLength(0.03, "XYZ");
  gStyle->SetNdivisions(510, "XYZ");
  gStyle->SetPadTickX(1);
  gStyle->SetPadTickY(1);
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetTitleOffset(1.00, "X");
  gStyle->SetTitleOffset(1.25, "Y");
  gStyle->SetLabelOffset(0.008, "XYZ");
  gStyle->SetLabelSize(0.04, "XYZ");
}
