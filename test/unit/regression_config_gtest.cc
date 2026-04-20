#include <gtest/gtest.h>

#include <string>
#include <vector>

#include "regression_common.h"

TEST(Regressions, InternalConfigRoundTrip) {
    ConfigureModelsForBaseline();

    char model[64] = {0};
    bool cartIn = false;
    bool cartOut = false;
    int maxDeg = 0;

    JupitermagGetInternalCFG(model, &cartIn, &cartOut, &maxDeg);
    const int requestedMaxDeg = (maxDeg > 3) ? (maxDeg - 2) : maxDeg;

    JupitermagSetInternalCFG("jrm33", false, true, requestedMaxDeg);

    char gotModel[64] = {0};
    bool gotCartIn = true;
    bool gotCartOut = false;
    int gotMaxDeg = 0;
    JupitermagGetInternalCFG(gotModel, &gotCartIn, &gotCartOut, &gotMaxDeg);

    EXPECT_STREQ("jrm33", gotModel);
    EXPECT_FALSE(gotCartIn);
    EXPECT_TRUE(gotCartOut);
    EXPECT_EQ(requestedMaxDeg, gotMaxDeg);

    JupitermagSetInternalCFG("jrm33", true, true, maxDeg);
}

TEST(Regressions, BaselineCSVSchema) {
    const auto fieldHeader = ReadCSVHeader(TestDataDir() / "field_baseline.csv");
    const std::vector<std::string> expectedFieldHeader = {
        "x",          "y",          "z",          "internal_bx",
        "internal_by", "internal_bz", "external_bx", "external_by",
        "external_bz", "combined_bx", "combined_by", "combined_bz"};
    EXPECT_EQ(expectedFieldHeader, fieldHeader);

    const auto traceHeader = ReadCSVHeader(TestDataDir() / "trace_summary.csv");
    const std::vector<std::string> expectedTraceHeader = {
        "trace_index",     "start_x",        "start_y",        "start_z",
        "ion_north_x",     "ion_north_y",    "ion_north_z",    "ion_south_x",
        "ion_south_y",     "ion_south_z",    "surf_north_x",   "surf_north_y",
        "surf_north_z",    "surf_south_x",   "surf_south_y",   "surf_south_z",
        "eq_mag_x",        "eq_mag_y",       "eq_mag_z",       "ion_north_lat",
        "ion_north_lon",   "ion_south_lat",  "ion_south_lon",  "surf_north_lat",
        "surf_north_lon",  "surf_south_lat", "surf_south_lon", "eq_mlon",
        "eq_lshell"};
    EXPECT_EQ(expectedTraceHeader, traceHeader);

    const auto fieldRows = ReadNumericCSV(TestDataDir() / "field_baseline.csv", 12);
    const auto traceRows = ReadNumericCSV(TestDataDir() / "trace_summary.csv", 29);
    EXPECT_FALSE(fieldRows.empty());
    EXPECT_FALSE(traceRows.empty());
}

TEST(Regressions, Con2020ConfigRoundTrip) {
    ConfigureModelsForBaseline();

    const double mui = 140.5;
    const double irho = 17.2;
    const double r0 = 8.1;
    const double r1 = 49.7;
    const double d = 3.5;
    const double xt = 8.9 * kDeg2Rad;
    const double xp = 158.2 * kDeg2Rad;
    const bool edwards = false;
    const bool errChk = true;
    const bool cartIn = false;
    const bool cartOut = true;
    const bool smooth = true;
    const double deltaRho = 0.9;
    const double deltaZ = 0.12;
    const double g = 417000.0;
    const double wOOpen = 0.2;
    const double wOOm = 0.4;
    const double thetamm = 15.8 * kDeg2Rad;
    const double dthetamm = 0.6 * kDeg2Rad;
    const double thetaoc = 10.9 * kDeg2Rad;
    const double dthetaoc = 0.15 * kDeg2Rad;

    JupitermagSetCon2020Params(mui, irho, r0, r1, d, xt, xp, "integral",
                               edwards, errChk, cartIn, cartOut, smooth,
                               deltaRho, deltaZ, g, "lmic", wOOpen, wOOm,
                               thetamm, dthetamm, thetaoc, dthetaoc);

    double gotMui = 0.0, gotIrho = 0.0, gotR0 = 0.0, gotR1 = 0.0, gotD = 0.0;
    double gotXt = 0.0, gotXp = 0.0;
    char gotEqtype[32] = {0};
    bool gotEdwards = true, gotErrChk = false, gotCartIn = true,
         gotCartOut = false, gotSmooth = false;
    double gotDeltaRho = 0.0, gotDeltaZ = 0.0, gotG = 0.0;
    char gotAzfunc[32] = {0};
    double gotWOOpen = 0.0, gotWOOm = 0.0, gotThetaMM = 0.0, gotDThetaMM = 0.0,
           gotThetaOC = 0.0, gotDThetaOC = 0.0;

    JupitermagGetCon2020Params(&gotMui, &gotIrho, &gotR0, &gotR1, &gotD,
                               &gotXt, &gotXp, gotEqtype, &gotEdwards,
                               &gotErrChk, &gotCartIn, &gotCartOut,
                               &gotSmooth, &gotDeltaRho, &gotDeltaZ, &gotG,
                               gotAzfunc, &gotWOOpen, &gotWOOm, &gotThetaMM,
                               &gotDThetaMM, &gotThetaOC, &gotDThetaOC);

    EXPECT_TRUE(NearlyEqual(mui, gotMui, kFieldAbsTol, kFieldRelTol));
    EXPECT_TRUE(NearlyEqual(irho, gotIrho, kFieldAbsTol, kFieldRelTol));
    EXPECT_TRUE(NearlyEqual(r0, gotR0, kFieldAbsTol, kFieldRelTol));
    EXPECT_TRUE(NearlyEqual(r1, gotR1, kFieldAbsTol, kFieldRelTol));
    EXPECT_TRUE(NearlyEqual(d, gotD, kFieldAbsTol, kFieldRelTol));
    EXPECT_TRUE(NearlyEqual(xt, gotXt, kFieldAbsTol, kFieldRelTol));
    EXPECT_TRUE(NearlyEqual(xp, gotXp, kFieldAbsTol, kFieldRelTol));
    EXPECT_STREQ("integral", gotEqtype);
    EXPECT_EQ(edwards, gotEdwards);
    EXPECT_EQ(errChk, gotErrChk);
    EXPECT_EQ(cartIn, gotCartIn);
    EXPECT_EQ(cartOut, gotCartOut);
    EXPECT_EQ(smooth, gotSmooth);
    EXPECT_TRUE(NearlyEqual(deltaRho, gotDeltaRho, kFieldAbsTol, kFieldRelTol));
    EXPECT_TRUE(NearlyEqual(deltaZ, gotDeltaZ, kFieldAbsTol, kFieldRelTol));
    EXPECT_TRUE(NearlyEqual(g, gotG, kFieldAbsTol, kFieldRelTol));
    EXPECT_STREQ("lmic", gotAzfunc);
    EXPECT_TRUE(NearlyEqual(wOOpen, gotWOOpen, kFieldAbsTol, kFieldRelTol));
    EXPECT_TRUE(NearlyEqual(wOOm, gotWOOm, kFieldAbsTol, kFieldRelTol));
    EXPECT_TRUE(NearlyEqual(thetamm, gotThetaMM, kFieldAbsTol, kFieldRelTol));
    EXPECT_TRUE(
        NearlyEqual(dthetamm, gotDThetaMM, kFieldAbsTol, kFieldRelTol));
    EXPECT_TRUE(NearlyEqual(thetaoc, gotThetaOC, kFieldAbsTol, kFieldRelTol));
    EXPECT_TRUE(
        NearlyEqual(dthetaoc, gotDThetaOC, kFieldAbsTol, kFieldRelTol));

    ConfigureModelsForBaseline();
}
