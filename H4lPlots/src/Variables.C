#ifndef VARIABLES
#define VARIABLES
#include <TMath.h>
#include <map>
#include <string>
#include <vector>

struct HistogramBinning {
  int nBins;
  float xlow, xhigh;
  float getBinSize() { return (xhigh - xlow) / nBins; }
};
const std::vector<HistogramBinning> higgsMassBinGranularities = {
    HistogramBinning{.nBins = 65, .xlow = 70., .xhigh = 200.},
    HistogramBinning{.nBins = 233, .xlow = 70., .xhigh = 1002.},
    HistogramBinning{.nBins = 93, .xlow = 70., .xhigh = 1000.},
};
const std::vector<HistogramBinning> zMassBinGranularities = {
    HistogramBinning{.nBins = 35, .xlow = 50., .xhigh = 120},
};
const std::vector<HistogramBinning> ptBinGranularities = {
    HistogramBinning{.nBins = 70, .xlow = 0., .xhigh = 350},
};
const std::vector<HistogramBinning> cosBinGranularities = {
    HistogramBinning{.nBins = 30, .xlow = -1, .xhigh = 1},
};
const std::vector<HistogramBinning> etaBinGranularities = {
    HistogramBinning{.nBins = 50, .xlow = -4, .xhigh = 4},
};
const std::vector<HistogramBinning> angleBinGranularities = {
    HistogramBinning{.nBins = 50, .xlow = -TMath::Pi(), .xhigh = TMath::Pi()},
};

const std::map<std::string, std::string> labelPerVariable = {
    {"ZZmass", "m_{#it{4l}}"},
    {"Z1mass", "m_{Z_1}"},
    {"Z2mass", "m_{Z_2}"},
    {"ZZpt", "p_{T}^{#it{4l}}"},
    {"ZZeta", "#eta_{#it{4l}}"},
    {"ZZphi", "#phi_{#it{4l}}"},
    {"ZZcostheta1", "cos#theta_{1}"},
    {"ZZcostheta2", "cos#theta_{2}"},
    {"ZZcosthetastar", "cos#theta^{*}"},
    {"ZZPhi", "#Phi_{#it{4l}}"},
    {"ZZPhi1", "#Phi1_{#it{4l}}"},
};

enum unitIdx { m_H, m_Z, GeV_c, angle, cosAngle, eta, Idx };
const std::map<unitIdx, std::vector<HistogramBinning>> binGranularites = {
    {m_H, higgsMassBinGranularities}, {m_Z, zMassBinGranularities},
    {GeV_c, ptBinGranularities},      {angle, angleBinGranularities},
    {cosAngle, cosBinGranularities},  {eta, etaBinGranularities},
};
const std::map<unitIdx, std::vector<std::string>> variablesPerUnit = {
    {m_H, {"ZZmass"}},
    {m_Z, {"Z1mass", "Z2mass"}},
    {GeV_c, {"ZZpt"}},
    {angle, {"ZZphi", "ZZPhi", "ZZPhi1"}},
    {cosAngle, {"ZZcostheta1", "ZZcostheta2", "ZZcosthetastar"}},
    {eta, {"ZZeta"}},
};

std::string getUnitName(unitIdx unit) {
  switch (unit) {
  case m_H:
  case m_Z:
    return "GeV";
  case GeV_c:
    return "GeV/c";
  case angle:
  case cosAngle:
  case eta:
  default:
    return "";
  }
}

std::string getUnitNameForTitle(unitIdx unit) {
  switch (unit) {
  case GeV_c:
    return "GeV_c";
  default:
    return getUnitName(unit);
  }
}

#endif
