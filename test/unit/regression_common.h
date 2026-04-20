#ifndef LIBJUPITERMAG_TEST_UNIT_REGRESSION_COMMON_H_
#define LIBJUPITERMAG_TEST_UNIT_REGRESSION_COMMON_H_

#include <array>
#include <filesystem>
#include <string>
#include <vector>

#include "jupitermag.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

inline constexpr double kDeg2Rad = M_PI / 180.0;
inline constexpr double kFieldAbsTol = 1e-6;
inline constexpr double kFieldRelTol = 1e-6;
#ifdef _MSC_VER
inline constexpr double kTraceAbsTol = 1e-4;
inline constexpr double kTraceRelTol = 1e-4;
#else
inline constexpr double kTraceAbsTol = 1e-6;
inline constexpr double kTraceRelTol = 1e-6;
#endif
inline constexpr double kGeomAbsTol = 1e-5;
inline constexpr double kGeomRelTol = 1e-5;

using Footprint49 = std::array<double, 49>;

bool NearlyEqual(double expected, double actual, double absTol, double relTol);
std::vector<std::vector<double>> ReadNumericCSV(const std::filesystem::path &file,
                                                size_t expectedColumns);
std::vector<std::string> ReadCSVHeader(const std::filesystem::path &file);
void ConfigureModelsForBaseline();
std::filesystem::path TestDataDir();
std::vector<std::array<double, 3>> RegressionPoints();
bool RunTraceWithDir(const std::vector<std::array<double, 3>> &starts, int traceDir,
                     std::vector<int> *nstepOut, std::vector<double *> *sOut,
                     std::vector<double *> *rnormOut,
                     std::vector<double *> *fpOut);
std::vector<Footprint49> RunTraceAndGetFootprints(
    const std::vector<std::array<double, 3>> &starts);

#endif
