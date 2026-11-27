#ifndef CHANNELNAMES_HXX
#define CHANNELNAMES_HXX

#include <array>
#include <cstddef>
#include <map>
#include <string>

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

template <typename T> struct View {
  const T *data;
  size_t n;

  constexpr View(const T *data, size_t n) : data(data), n(n) {}
  template <size_t N>
  constexpr View(const std::array<T, N> &arr) : data(arr.data()), n(N) {}
  template <size_t N> constexpr View(const std::array<T, N> &&) = delete;

  constexpr const T *begin() const { return data; }
  constexpr const T *end() const { return data + n; }
  constexpr const T &operator[](size_t i) const { return data[i]; }
  constexpr size_t size() const { return n; }
};

template <typename T, size_t N>
constexpr View<T> arrayToView(const std::array<T, N> &arr) {
  return {arr.data(), N};
}

// MC
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
constexpr std::array<const char *, 1> MC_BkgDYUnited = {"DYJetsToLL"};

constexpr std::array<const char *, 2> MC_BkgTT = {"TTWW", "TTZZ"};
constexpr std::array<const char *, 1> MC_BkgTTTo = {"TTto2L2Nu"};

constexpr std::array<const char *, 46> MC_dataSetsFullH125 =
    concatAll(MC_SignalH124, MC_SignalH124p5, MC_SignalH125, MC_SignalH125p5,
              MC_SignalH126, MC_BkgEW, MC_BkgQQZZ, MC_BkgGGTo, MC_BkgDY,
              MC_BkgTT, MC_BkgTTTo);

constexpr std::array<const char *, 45> MC_dataSetsFullH125UnitedDY =
    concatAll(MC_SignalH124, MC_SignalH124p5, MC_SignalH125, MC_SignalH125p5,
              MC_SignalH126, MC_BkgEW, MC_BkgQQZZ, MC_BkgGGTo, MC_BkgDYUnited,
              MC_BkgTT, MC_BkgTTTo);

constexpr std::array<const char *, 22> MC_dataSetsH124 =
    concatAll(MC_SignalH124, MC_BkgEW, MC_BkgQQZZ, MC_BkgGGTo, MC_BkgDY,
              MC_BkgTT, MC_BkgTTTo);

constexpr std::array<const char *, 22> MC_dataSetsH124p5 =
    concatAll(MC_SignalH124p5, MC_BkgEW, MC_BkgQQZZ, MC_BkgGGTo, MC_BkgDY,
              MC_BkgTT, MC_BkgTTTo);

constexpr std::array<const char *, 22> MC_dataSetsH125 =
    concatAll(MC_SignalH125, MC_BkgEW, MC_BkgQQZZ, MC_BkgGGTo, MC_BkgDY,
              MC_BkgTT, MC_BkgTTTo);

constexpr std::array<const char *, 22> MC_dataSetsH125p5 =
    concatAll(MC_SignalH125p5, MC_BkgEW, MC_BkgQQZZ, MC_BkgGGTo, MC_BkgDY,
              MC_BkgTT, MC_BkgTTTo);

constexpr std::array<const char *, 22> MC_dataSetsH126 =
    concatAll(MC_SignalH126, MC_BkgEW, MC_BkgQQZZ, MC_BkgGGTo, MC_BkgDY,
              MC_BkgTT, MC_BkgTTTo);
constexpr std::array<const char *, 1> MC_WZTo =
    std::array<const char *, 1>{"WZto3LNu"};

constexpr std::array<const char *, 12> MC_ControlRegion = concatAll(
    MC_BkgDY, MC_BkgDYUnited, MC_BkgGGTo, MC_BkgTTTo, MC_WZTo, MC_BkgQQZZ);

// Data
constexpr std::array<const char *, 2> Data_EGamma2022 = {
    "EGamma2022C",
    "EGamma2022D",
};

constexpr std::array<const char *, 3> Data_Muon2022 = {
    "Muon2022C",
    "Muon2022D",
    "DoubleMu2022C",
};

constexpr std::array<const char *, 3> Data_MuonEG2022 = {
    "MuonEG2022C",
    "MuonEG2022D",
    "SingleMu2022C",
};
constexpr std::array<const char *, 8> Data_dataSets2022 =
    concatAll(Data_EGamma2022, Data_Muon2022, Data_MuonEG2022);

constexpr std::array<const char *, 3> Data_EGamma2022EE = {
    "EGamma2022E",
    "EGamma2022F",
    "EGamma2022G",
};

constexpr std::array<const char *, 3> Data_Muon2022EE = {
    "Muon2022E",
    "Muon2022F",
    "Muon2022G",
};

constexpr std::array<const char *, 3> Data_MuonEG2022EE = {
    "MuonEG2022E",
    "MuonEG2022F",
    "MuonEG2022G",
};

constexpr std::array<const char *, 9> Data_dataSets2022EE =
    concatAll(Data_EGamma2022EE, Data_Muon2022EE, Data_MuonEG2022EE);

constexpr std::array<const char *, 8> Data_EGamma2023preBPix = {
    "EGamma02023Cv1", "EGamma02023Cv2", "EGamma02023Cv3", "EGamma02023Cv4",
    "EGamma12023Cv1", "EGamma12023Cv2", "EGamma12023Cv3", "EGamma12023Cv4",
};

constexpr std::array<const char *, 8> Data_Muon2023preBPix = {
    "Muon02023Cv1", "Muon02023Cv2", "Muon02023Cv3", "Muon02023Cv4",
    "Muon12023Cv1", "Muon12023Cv2", "Muon12023Cv3", "Muon12023Cv4",
};

constexpr std::array<const char *, 4> Data_MuonEG2023preBPix = {
    "MuonEG2023Cv1",
    "MuonEG2023Cv2",
    "MuonEG2023Cv3",
    "MuonEG2023Cv4",
};

constexpr std::array<const char *, 20> Data_dataSets2023preBPix = concatAll(
    Data_EGamma2023preBPix, Data_Muon2023preBPix, Data_MuonEG2023preBPix);

constexpr std::array<const char *, 4> Data_EGamma2023postBPix = {
    "EGamma02023Dv1",
    "EGamma02023Dv2",
    "EGamma12023Dv1",
    "EGamma12023Dv2",
};

constexpr std::array<const char *, 4> Data_Muon2023postBPix = {
    "Muon02023Dv1",
    "Muon02023Dv2",
    "Muon12023Dv1",
    "Muon12023Dv2",
};

constexpr std::array<const char *, 2> Data_MuonEG2023postBPix = {
    "MuonEG2023Dv1",
    "MuonEG2023Dv2",
};

constexpr std::array<const char *, 10> Data_dataSets2023postBPix = concatAll(
    Data_EGamma2023postBPix, Data_Muon2023postBPix, Data_MuonEG2023postBPix);

constexpr std::array<const char *, 16> Data_EGamma2024 = {
    "EGamma02024Cv1", "EGamma02024Dv1", "EGamma02024Ev1", "EGamma02024Fv1",
    "EGamma02024Gv2", "EGamma02024Hv2", "EGamma02024Iv1", "EGamma02024Iv2",
    "EGamma12024Cv1", "EGamma12024Dv1", "EGamma12024Ev1", "EGamma12024Fv1",
    "EGamma12024Gv2", "EGamma12024Hv1", "EGamma12024Iv1", "EGamma12024Iv2",
};

constexpr std::array<const char *, 16> Data_Muon2024 = {
    "Muon02024Cv1", "Muon02024Dv1", "Muon02024Ev1", "Muon02024Fv1",
    "Muon02024Gv1", "Muon02024Hv1", "Muon02024Iv1", "Muon02024Iv2",
    "Muon12024Cv1", "Muon12024Dv1", "Muon12024Ev1", "Muon12024Fv1",
    "Muon12024Gv2", "Muon12024Hv2", "Muon12024Iv1", "Muon12024Iv2",
};

constexpr std::array<const char *, 8> Data_MuonEG2024 = {
    "MuonEG2024Cv1", "MuonEG2024Dv1", "MuonEG2024Ev1", "MuonEG2024Fv2",
    "MuonEG2024Gv3", "MuonEG2024Hv2", "MuonEG2024Iv2", "MuonEG2024Iv2v2",
};
constexpr std::array<const char *, 40> Data_dataSets2024 =
    concatAll(Data_EGamma2024, Data_Muon2024, Data_MuonEG2024);

constexpr std::array<const char *, 32> Data_EGamma2025 = {
    "EGamma02025Bv1", "EGamma02025Cv1", "EGamma02025Cv2", "EGamma02025Dv1",
    "EGamma02025Ev1", "EGamma12025Bv1", "EGamma02025Fv1", "EGamma02025Fv2",
    "EGamma02025Gv1", "EGamma12025Cv1", "EGamma12025Cv2", "EGamma12025Dv1",
    "EGamma12025Ev1", "EGamma12025Fv1", "EGamma12025Fv2", "EGamma12025Gv1",
    "EGamma22025Bv1", "EGamma22025Cv1", "EGamma22025Cv2", "EGamma22025Dv1",
    "EGamma22025Ev1", "EGamma22025Fv1", "EGamma22025Fv2", "EGamma22025Gv1",
    "EGamma32025Bv1", "EGamma32025Cv1", "EGamma32025Cv2", "EGamma32025Dv1",
    "EGamma32025Ev1", "EGamma32025Fv1", "EGamma32025Fv2", "EGamma32025Gv1",
};

constexpr std::array<const char *, 16> Data_Muon2025 = {
    "Muon02025Bv1", "Muon02025Cv1", "Muon02025Cv2", "Muon02025Dv1",
    "Muon02025Ev1", "Muon02025Fv1", "Muon02025Fv2", "Muon02025Gv1",
    "Muon12025Bv1", "Muon12025Cv1", "Muon12025Cv2", "Muon12025Dv1",
    "Muon12025Ev1", "Muon12025Fv1", "Muon12025Fv2", "Muon12025Gv1"};

constexpr std::array<const char *, 8> Data_MuonEG2025 = {
    "MuonEG2025Bv1", "MuonEG2025Cv1", "MuonEG2025Cv2", "MuonEG2025Dv1",
    "MuonEG2025Ev1", "MuonEG2025Fv1", "MuonEG2025Fv2", "MuonEG2025Gv1",
};

const std::array<const char *, 56> Data_dataSets2025 =
    concatAll(Data_EGamma2025, Data_Muon2025, Data_MuonEG2025);

constexpr std::array<View<const char *>, 2> MC_dataSets = {
    MC_dataSetsFullH125, MC_dataSetsFullH125UnitedDY};
enum mcDataSet { separatedDY, unitedDY };

constexpr std::array<View<const char *>, 6> Data_dataSets{
    Data_dataSets2022,       Data_dataSets2022EE, Data_EGamma2023preBPix,
    Data_EGamma2023postBPix, Data_dataSets2024,   Data_dataSets2025};
enum dataPeriod { y2022, y2022EE, y2023preBPix, y2023postBPix, y2024, y2025 };

constexpr mcDataSet getMCIdx(dataPeriod year) {
  switch (year) {
  case y2022:
  case y2022EE:
  case y2023preBPix:
  case y2023postBPix:
    return unitedDY;
  case y2024:
  case y2025:
    return separatedDY;
  }
}

const std::map<std::string, dataPeriod> periodIndex = {
    {"2022", y2022},
    {"2022EE", y2022EE},
    {"2023preBPix", y2023preBPix},
    {"2023postBPix", y2023postBPix},
    {"2024", y2024},
    {"2025", y2025}};

enum class FinalStateFlavor {
  L4,
  E3 = 1331,
  E2Mu1 = 1573,
  Mu2E1 = 1859,
  Mu3 = 2197,
  E4 = 14641,
  E2Mu2 = 20449,
  Mu4 = 28561
};

constexpr FinalStateFlavor fsFlavors[] = {
    FinalStateFlavor::L4, FinalStateFlavor::E4, FinalStateFlavor::E2Mu2,
    FinalStateFlavor::Mu4};

constexpr FinalStateFlavor fsFlavorsZL[] = {
    FinalStateFlavor::L4, FinalStateFlavor::E3, FinalStateFlavor::E2Mu1,
    FinalStateFlavor::Mu2E1, FinalStateFlavor::Mu3};

#endif
