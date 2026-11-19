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
const std::vector<HistogramBinning> leptonMassBinGranularities = {
    HistogramBinning{.nBins = 35, .xlow = 0., .xhigh = 2},
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
    {"mass", "m_{#it{4l}}"},
    {"Z1mass", "m_{Z_1}"},
    {"Z2mass", "m_{Z_2}"},
    {"pt", "p_{T}^{#it{4l}}"},
    {"eta", "#eta_{#it{4l}}"},
    {"phi", "#phi_{#it{4l}}"},
    {"costheta1", "cos#theta_{1}"},
    {"costheta2", "cos#theta_{2}"},
    {"costhetastar", "cos#theta^{*}"},
    {"Phi", "#Phi_{#it{4l}}"},
    {"Phi1", "#Phi1_{#it{4l}}"},
};

enum unitIdx { m_H, m_Z, m_l, GeV_c, angle, cosAngle, eta, Idx };
const std::map<unitIdx, std::vector<HistogramBinning>> binGranularites = {
    {m_H, higgsMassBinGranularities},  {m_Z, zMassBinGranularities},
    {GeV_c, ptBinGranularities},       {angle, angleBinGranularities},
    {cosAngle, cosBinGranularities},   {eta, etaBinGranularities},
    {m_l, leptonMassBinGranularities},
};

const std::map<unitIdx, std::vector<std::string>> variablesPerUnit = {
    {m_H, {"mass"}},
    {m_Z, {"Z1mass", "Z2mass"}},
    {GeV_c, {"pt"}},
    {angle, {"phi", "Phi", "Phi1"}},
    {cosAngle, {"costheta1", "costheta2", "costhetastar"}},
    {eta, {"eta"}},
};
const std::map<unitIdx, std::vector<std::string>> variablesPerUnitZL = {
    {m_Z, {"mass"}},
    {m_l, {"leptonmass"}},
    {GeV_c, {"pt", "leptonpt"}},
    {angle, {"phi", "leptonphi"}},
    {eta, {"eta", "leptoneta"}},
};

const std::map<std::string, std::string> labelPerVariableZL = {
    {"mass", "m_{Z}"},
    {"pt", "p_{T}^{Z}"},
    {"eta", "#eta_{Z}"},
    {"phi", "#phi_{Z}"},
    {"leptonmass", "m_{#it{l}}"},
    {"leptonpt", "p_{T}^{#it{l}}"},
    {"leptoneta", "#eta_{#it{l}}"},
    {"leptonphi", "#phi_{#it{l}}"},
};

std::string getUnitName(unitIdx unit) {
  switch (unit) {
  case m_H:
  case m_Z:
  case m_l:
    return "GeV";
  case GeV_c:
    return "GeV/c";
  case angle:
  case cosAngle:
  case eta:
  case Idx:
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

enum CRIdx { P2F2, P3F1, SS, ZL };
const std::map<CRIdx, std::string> CRNames = {
    {P2F2, "2P2F"},
    {P3F1, "3P1F"},
    {SS, "SS"},
    {ZL, "ZL"},
};
std::string getVariablePrefix(bool isCR, CRIdx crIdx) {
  if (isCR) {
    if (crIdx == ZL) {
      return "ZL";
    }
    return "ZLL" + CRNames.at(crIdx);
  }
  return "ZZ";
}
#endif
