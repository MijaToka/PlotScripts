#ifndef PATHS_HXX
#define PATHS_HXX

#include <filesystem>
#include <map>
#include <string>

const std::map<std::string, std::filesystem::path> MC_PATHS = {
    {"2025", "/eos/user/m/mtokarev/samplesHZZ/MC_2025_d25e5b4/"},
};
const std::map<std::string, std::filesystem::path> DATA_PATHS = {
    {"2025", "/eos/user/m/mtokarev/samplesHZZ/Data_2025_d25e5b4/"},
};

const std::map<std::string, std::filesystem::path> MC_FILL_FILES = {
    {"2025",
     "/eos/user/m/mtokarev/work_area/currentCMSSW/output/H4lFill_MC2025.root"}};

const std::map<std::string, std::filesystem::path> DATA_FILL_FILES = {
    {"2025", "/eos/user/m/mtokarev/work_area/currentCMSSW/output/"
             "H4lFill_Data2025.root"}};

const std::filesystem::path
    REPO_PATH("/eos/user/m/mtokarev/work_area/currentCMSSW/src/");

#endif
