#ifndef PATHS_HXX
#define PATHS_HXX

#include <filesystem>
#include <map>
#include <string>

inline const std::filesystem::path CJLST_PATH =
    "/eos/cms/store/group/phys_higgs/cmshzz4l/cjlst";
inline const std::filesystem::path PRIVATE_SAMPLES =
    "/eos/user/m/mtokarev/samplesHZZ";

const std::map<std::string, std::filesystem::path> MC_PATHS = {
    {"2022", CJLST_PATH / "HIG-26-STXS/"
                          "260623/2022/"},
    {"2022EE", CJLST_PATH / "HIG-26-STXS/"
                            "260623/2022EE/"},
    {"2023preBPix", CJLST_PATH / "HIG-26-STXS/"
                                 "260623/2023preBPix/"},
    {"2023postBPix", CJLST_PATH / "HIG-26-STXS/"
                                  "260623/2023postBPix/"},
    {"2024", CJLST_PATH / "HIG-26-STXS/"
                          "260623/2024/"},
    {"2025", PRIVATE_SAMPLES / "MC_2025_d25e5b4/"},
};
const std::map<std::string, std::filesystem::path> DATA_PATHS = {
    {"2022", CJLST_PATH / "HIG-26-STXS/"
                          "260623/Data2022/"},
    {"2022EE", CJLST_PATH / "HIG-26-STXS/"
                            "260623/Data2022/"},
    {"2023preBPix", CJLST_PATH / "HIG-26-STXS/"
                                 "260623/Data2023/"},
    {"2023postBPix", CJLST_PATH / "HIG-26-STXS/"
                                  "260623/Data2023/"},
    {"2024", CJLST_PATH / "HIG-26-STXS/"
                          "260623/Data2024/"},
    {"2025", PRIVATE_SAMPLES / "Data_2025_10ba8e7/"},
};

const std::filesystem::path
    REPO_PATH("/eos/user/m/mtokarev/work_area/currentCMSSW/src/");

const std::map<std::string, std::filesystem::path> MC_FILL_FILES = {
    {"2022", REPO_PATH / "outputs" / "Fill" / "H4lFill_MC2022.root"},
    {"2022EE", REPO_PATH / "outputs" / "Fill" / "H4lFill_MC2022EE.root"},
    {"2023preBPix",
     REPO_PATH / "outputs" / "Fill" / "H4lFill_MC2023preBPix.root"},
    {"2023postBPix",
     REPO_PATH / "outputs" / "Fill" / "H4lFill_MC2023postBPix.root"},
    {"2024", REPO_PATH / "outputs" / "Fill" / "H4lFill_MC2024.root"},
    {"2025", REPO_PATH / "outputs" / "Fill" / "H4lFill_MC2025.root"},
};

const std::map<std::string, std::filesystem::path> DATA_FILL_FILES = {
    {"2022", REPO_PATH / "outputs" / "Fill" / "H4lFill_Data2022.root"},
    {"2022EE", REPO_PATH / "outputs" / "Fill" / "H4lFill_Data2022EE.root"},
    {"2023preBPix",
     REPO_PATH / "outputs" / "Fill" / "H4lFill_Data2023preBPix.root"},
    {"2023postBPix",
     REPO_PATH / "outputs" / "Fill" / "H4lFill_Data2023postBPix.root"},
    {"2024", REPO_PATH / "outputs" / "Fill" / "H4lFill_Data2024.root"},
    {"2025", REPO_PATH / "outputs" / "Fill" / "H4lFill_Data2025.root"},
};
#endif
