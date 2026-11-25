#ifndef VARIABLES
#define VARIABLES

#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <Math/Vector4Dfwd.h>
#include <ROOT/RDF/InterfaceUtils.hxx>
#include <ROOT/RDFHelpers.hxx>
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <RtypesCore.h>
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
    HistogramBinning{.nBins = 35, .xlow = 0., .xhigh = 0.2},
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
const std::vector<HistogramBinning> METBinGranularities = {
    HistogramBinning{.nBins = 50, .xlow = 0, .xhigh = 100},
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

enum unitIdx { m_H, m_Z, m_l, MET, GeV_c, angle, cosAngle, eta, Idx };
const std::map<unitIdx, std::vector<HistogramBinning>> binGranularites = {
    {m_H, higgsMassBinGranularities},  {m_Z, zMassBinGranularities},
    {GeV_c, ptBinGranularities},       {angle, angleBinGranularities},
    {cosAngle, cosBinGranularities},   {eta, etaBinGranularities},
    {m_l, leptonMassBinGranularities}, {MET, METBinGranularities},
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
    {MET, {"ETmiss"}},
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
    {"ETmiss", "E_{T}^{miss}"},
};

template <typename T>
inline T flattenVariable(const ROOT::VecOps::RVec<T> &Cand, Short_t bestIdx) {
  return bestIdx >= 0 ? Cand[bestIdx] : -999.;
}

#define SR_FLATTEN_VARIABLE(columnTag, variableType)                           \
  "ZZ" columnTag, flattenVariable<variableType>,                               \
      {"ZZCand_" columnTag, "bestCandIdx"}

#define CR_FLATTEN_VARIABLE(columnTag, variableType, region)                   \
  "ZLL" region columnTag, flattenVariable<variableType>,                       \
      {"ZLLCand_" columnTag, "ZLLbest" region "Idx"}

#define CR2P2F_FLATTEN_VARIABLE(columnTag, variableType)                       \
  CR_FLATTEN_VARIABLE(columnTag, variableType, "2P2F")

#define CR3P1F_FLATTEN_VARIABLE(columnTag, variableType)                       \
  CR_FLATTEN_VARIABLE(columnTag, variableType, "3P1F")

#define CRSS_FLATTEN_VARIABLE(columnTag, variableType)                         \
  CR_FLATTEN_VARIABLE(columnTag, variableType, "SS")

#define ZL_FLATTEN_VARIABLE(columnTag, variableType)                           \
  "ZL" columnTag, flattenVariable<variableType>,                               \
      {"ZCand_" columnTag, "bestZIdx"}

#define LEP_FLATTEN_VARIABLE(columnTag, variableType)                          \
  "ZLlepton" columnTag, flattenVariable<variableType>,                         \
      {"Lepton_" columnTag, "ZLCand_lepIdx"}

#define DEFINE_SINGLE_VARIABLES(X)                                             \
  .Define(X("mass", Float_t))                                                  \
      .Define(X("pt", Float_t))                                                \
      .Define(X("eta", Float_t))                                               \
      .Define(X("phi", Float_t))

#define DEFINE_ZZ_VARIABLES(X)                                                 \
  .Define(X("Z1mass", Float_t))                                                \
      .Define(X("Z2mass", Float_t))                                            \
      .Define(X("costheta1", Float_t))                                         \
      .Define(X("costheta2", Float_t))                                         \
      .Define(X("costhetastar", Float_t))                                      \
      .Define(X("Phi", Float_t))                                               \
      .Define(X("Phi1", Float_t))

#define DEFINE_VARIABLES(X) DEFINE_SINGLE_VARIABLES(X) DEFINE_ZZ_VARIABLES(X)

template <bool hasCR>
ROOT::RDF::RNode initializeVariables(ROOT::RDF::RNode df) {
  df = df.Define("ZZFSflav",
                 [](ROOT::RVecI Z1Flav, ROOT::RVecI Z2Flav, Short_t bestIdx) {
                   return bestIdx >= 0 ? Z1Flav[bestIdx] * Z2Flav[bestIdx] : 0;
                 },
                 {"ZZCand_Z1flav", "ZZCand_Z2flav", "bestCandIdx"})
           DEFINE_VARIABLES(SR_FLATTEN_VARIABLE);

  if constexpr (hasCR) {
    df =
        df.Define("ZLL2P2FFSflav",
                  [](ROOT::RVecI Z1Flav, ROOT::RVecI Z2Flav, Short_t bestIdx) {
                    return bestIdx >= 0 ? Z1Flav[bestIdx] * Z2Flav[bestIdx] : 0;
                  },
                  {"ZLLCand_Z1flav", "ZLLCand_Z2flav", "ZLLbest2P2FIdx"})
            .Define(
                "ZLL3P1FFSflav",
                [](ROOT::RVecI Z1Flav, ROOT::RVecI Z2Flav, Short_t bestIdx) {
                  return bestIdx >= 0 ? Z1Flav[bestIdx] * Z2Flav[bestIdx] : 0;
                },
                {"ZLLCand_Z1flav", "ZLLCand_Z2flav", "ZLLbest3P1FIdx"})
            .Define(
                "ZLLSSFSflav",
                [](ROOT::RVecI Z1Flav, ROOT::RVecI Z2Flav, Short_t bestIdx) {
                  return bestIdx >= 0 ? Z1Flav[bestIdx] * Z2Flav[bestIdx] : 0;
                },
                {"ZLLCand_Z1flav", "ZLLCand_Z2flav", "ZLLbestSSIdx"})
            .Define("ZLFSflav",
                    [](ROOT::RVecI ZFlav, ROOT::RVecI LepFlav, Short_t bestZIdx,
                       Short_t bestZLIdx) {
                      return bestZLIdx >= 0
                                 ? ZFlav[bestZIdx] * LepFlav[bestZLIdx]
                                 : 0.;
                    },
                    {"ZCand_flav", "Lepton_pdgId", "bestZIdx", "ZLCand_lepIdx"})

                DEFINE_VARIABLES(CR2P2F_FLATTEN_VARIABLE)
                    DEFINE_VARIABLES(CR3P1F_FLATTEN_VARIABLE)
                        DEFINE_VARIABLES(CRSS_FLATTEN_VARIABLE)
                            DEFINE_SINGLE_VARIABLES(ZL_FLATTEN_VARIABLE)
                                DEFINE_SINGLE_VARIABLES(LEP_FLATTEN_VARIABLE)
            .Define("ZLETmiss",
                    [](Float_t Zmass, Float_t Zpt, Float_t Zeta, Float_t Zphi,
                       Float_t Lmass, Float_t Lpt, Float_t Leta, Float_t Lphi,
                       Short_t bestZLIdx) -> Float_t {
                      if (bestZLIdx == -1)
                        return -999.;

                      ROOT::Math::PtEtaPhiMVector Z(Zpt, Zeta, Zphi, Zmass);
                      ROOT::Math::PtEtaPhiMVector L(Lpt, Leta, Lphi, Lmass);

                      ROOT::Math::PtEtaPhiMVector missing = (Z + L);
                      return missing.Pt();
                    },
                    {"ZLmass", "ZLpt", "ZLeta", "ZLphi", "ZLleptonmass",
                     "ZLleptonpt", "ZLleptoneta", "ZLleptonphi",
                     "ZLCand_lepIdx"});
  }

  return df;
}

std::string getUnitName(unitIdx unit) {
  switch (unit) {
  case m_H:
  case m_Z:
  case m_l:
  case MET:
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
