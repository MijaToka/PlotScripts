#ifndef CHANNELNAMES_HXX
#define CHANNELNAMES_HXX

#include <array>
#include <cstddef>

template <typename T, size_t LLen, size_t RLen>
constexpr std::array<T, LLen + RLen> concat(const std::array<T, LLen> &a,
                                            const std::array<T, RLen> &b) {
  std::array<T, LLen + RLen> out{};
  for (size_t i = 0; i < LLen; ++i)
    out[i] = a[i];
  for (size_t i = 0; i < RLen; ++i)
    out[LLen + i] = b[i];
  return out;
}

template <typename A, typename B, typename... Rest>
constexpr auto concatAll(const A &a, const B &b, const Rest &...rest) {
  if constexpr (sizeof...(rest) == 0)
    return concat(a, b);
  else
    return concatAll(concat(a, b), rest...);
}

constexpr std::array<const char *, 6> MC_SignalH124 = {
    "ggH124", "VBFH124", "WplusH124", "WminusH124", "ZH124", "ttH124"};

constexpr std::array<const char *, 6> MC_SignalH124p5 = {
    "ggH124p5",     "VBFH124p5", "WplusH124p5",
    "WminusH124p5", "ZH124p5",   "ttH124p5"};

constexpr std::array<const char *, 6> MC_SignalH125 = {
    "ggH125", "VBFH125", "WplusH125", "WminusH125", "ZH125", "ttH125"};

constexpr std::array<const char *, 6> MC_SignalH125p5 = {
    "ggH125p5",     "VBFH125p5", "WplusH125p5",
    "WminusH125p5", "ZH125p5",   "ttH125p5"};

constexpr std::array<const char *, 6> MC_SignalH126 = {
    "ggH126", "VBFH126", "WplusH126", "WminusH126", "ZH126", "ttH126"};

constexpr std::array<const char *, 4> MC_BkgEW = {"WZto3LNu", "WWZ", "WZZ",
                                                  "ZZZ"};

constexpr std::array<const char *, 1> MC_BkgQQZZ = {"ZZTo4l"};

constexpr std::array<const char *, 6> MC_BkgGGTo = {
    "ggTo4e_Contin_MCFM701",      "ggTo2e2mu_Contin_MCFM701",
    "ggTo2e2tau_Contin_MCFM701",  "ggTo4mu_Contin_MCFM701",
    "ggTo2mu2tau_Contin_MCFM701", "ggTo4tau_Contin_MCFM701"};

constexpr std::array<const char *, 2> MC_BkgDY = {"DYJetsTo2E", "DYJetsTo2Mu"};

constexpr std::array<const char *, 3> MC_BkgTT = {"TTto2L2Nu", "TTWW", "TTZZ"};

constexpr std::array<const char *, 46> MC_dataSets = concatAll(
    MC_SignalH124, MC_SignalH124p5, MC_SignalH125, MC_SignalH125p5,
    MC_SignalH126, MC_BkgEW, MC_BkgQQZZ, MC_BkgGGTo, MC_BkgDY, MC_BkgTT);

constexpr std::array<const char *, 22> MC_dataSetsH124 = concatAll(
    MC_SignalH124, MC_BkgEW, MC_BkgQQZZ, MC_BkgGGTo, MC_BkgDY, MC_BkgTT);

constexpr std::array<const char *, 22> MC_dataSetsH124p5 = concatAll(
    MC_SignalH124p5, MC_BkgEW, MC_BkgQQZZ, MC_BkgGGTo, MC_BkgDY, MC_BkgTT);

constexpr std::array<const char *, 22> MC_dataSetsH125 = concatAll(
    MC_SignalH125, MC_BkgEW, MC_BkgQQZZ, MC_BkgGGTo, MC_BkgDY, MC_BkgTT);

constexpr std::array<const char *, 22> MC_dataSetsH125p5 = concatAll(
    MC_SignalH125p5, MC_BkgEW, MC_BkgQQZZ, MC_BkgGGTo, MC_BkgDY, MC_BkgTT);

constexpr std::array<const char *, 22> MC_dataSetsH126 = concatAll(
    MC_SignalH126, MC_BkgEW, MC_BkgQQZZ, MC_BkgGGTo, MC_BkgDY, MC_BkgTT);

constexpr std::array<const char *, 14> Data_EGamma2025 = {
    "EGamma02025Cv1", "EGamma02025Cv2", "EGamma02025Dv1", "EGamma02025Ev1",
    "EGamma02025Fv1", "EGamma02025Fv2", "EGamma02025Gv1", "EGamma12025Cv1",
    "EGamma12025Cv2", "EGamma12025Dv1", "EGamma12025Ev1", "EGamma12025Fv1",
    "EGamma12025Fv2", "EGamma12025Gv1"};

constexpr std::array<const char *, 14> Data_Muon2025 = {
    "Muon02025Cv1", "Muon02025Cv2", "Muon02025Dv1", "Muon02025Ev1",
    "Muon02025Fv1", "Muon02025Fv2", "Muon02025Gv1", "Muon12025Cv1",
    "Muon12025Cv2", "Muon12025Dv1", "Muon12025Ev1", "Muon12025Fv1",
    "Muon12025Fv2", "Muon12025Gv1"};

constexpr std::array<const char *, 7> Data_MuonEG2025 = {
    "MuonEG2025Cv1", "MuonEG2025Cv2", "MuonEG2025Dv1", "MuonEG2025Ev1",
    "MuonEG2025Fv1", "MuonEG2025Fv2", "MuonEG2025Gv1",
};

const std::array<const char *, 35> Data_dataSets =
    concatAll(Data_EGamma2025, Data_Muon2025, Data_MuonEG2025);

enum class FinalStateFlavor { L4, E4 = 14641, E2Mu2 = 20449, Mu4 = 28561 };

constexpr FinalStateFlavor fsFlavors[] = {
    FinalStateFlavor::L4, FinalStateFlavor::E4, FinalStateFlavor::E2Mu2,
    FinalStateFlavor::Mu4};

#endif
