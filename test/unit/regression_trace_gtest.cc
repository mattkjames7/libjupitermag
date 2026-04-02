#include <gtest/gtest.h>

#include <cmath>
#include <vector>

#include "regression_common.h"

TEST(Regressions, TraceFootprintBaselineCSV) {
    ConfigureModelsForBaseline();

    const auto file = TestDataDir() / "trace_summary.csv";
    const auto rows = ReadNumericCSV(file, 29);
    ASSERT_FALSE(rows.empty());

    std::vector<std::array<double, 3>> starts(rows.size());
    for (size_t i = 0; i < rows.size(); i++) {
        starts[i] = {rows[i][1], rows[i][2], rows[i][3]};
    }
    const auto fps = RunTraceAndGetFootprints(starts);

    for (size_t i = 0; i < rows.size(); i++) {
        const auto &row = rows[i];
        const auto &traceFp = fps[i];

        EXPECT_TRUE(NearlyEqual(row[4], traceFp[0], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[5], traceFp[1], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[6], traceFp[2], kTraceAbsTol, kTraceRelTol))
            << i;

        EXPECT_TRUE(NearlyEqual(row[7], traceFp[3], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[8], traceFp[4], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[9], traceFp[5], kTraceAbsTol, kTraceRelTol))
            << i;

        EXPECT_TRUE(NearlyEqual(row[10], traceFp[12], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[11], traceFp[13], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[12], traceFp[14], kTraceAbsTol, kTraceRelTol))
            << i;

        EXPECT_TRUE(NearlyEqual(row[13], traceFp[15], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[14], traceFp[16], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[15], traceFp[17], kTraceAbsTol, kTraceRelTol))
            << i;

        EXPECT_TRUE(NearlyEqual(row[16], traceFp[27], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[17], traceFp[28], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[18], traceFp[29], kTraceAbsTol, kTraceRelTol))
            << i;

        EXPECT_TRUE(NearlyEqual(row[19], traceFp[31], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[20], traceFp[30], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[21], traceFp[35], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[22], traceFp[34], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[23], traceFp[39], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[24], traceFp[38], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[25], traceFp[43], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[26], traceFp[42], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[27], traceFp[47], kTraceAbsTol, kTraceRelTol))
            << i;
        EXPECT_TRUE(NearlyEqual(row[28], traceFp[46], kTraceAbsTol, kTraceRelTol))
            << i;
    }
}

TEST(Regressions, TraceWrapperTraceDirContracts) {
    ConfigureModelsForBaseline();

    const std::vector<std::array<double, 3>> starts = {{5.0, 0.0, 0.0}};

    std::vector<int> nstep;
    std::vector<double *> s;
    std::vector<double *> rnorm;
    std::vector<double *> fp;

    EXPECT_TRUE(RunTraceWithDir(starts, -1, &nstep, &s, &rnorm, &fp));
    ASSERT_EQ(1, static_cast<int>(nstep.size()));
    EXPECT_GT(nstep[0], 0);
    EXPECT_TRUE(std::isnan(s[0][0]));
    EXPECT_TRUE(std::isnan(rnorm[0][0]));
    EXPECT_TRUE(std::isnan(fp[0][0]));
    delete[] s[0];
    delete[] rnorm[0];
    delete[] fp[0];

    EXPECT_TRUE(RunTraceWithDir(starts, 1, &nstep, &s, &rnorm, &fp));
    ASSERT_EQ(1, static_cast<int>(nstep.size()));
    EXPECT_GT(nstep[0], 0);
    EXPECT_TRUE(std::isnan(s[0][0]));
    EXPECT_TRUE(std::isnan(rnorm[0][0]));
    EXPECT_TRUE(std::isnan(fp[0][0]));
    delete[] s[0];
    delete[] rnorm[0];
    delete[] fp[0];

    EXPECT_TRUE(RunTraceWithDir(starts, 0, &nstep, &s, &rnorm, &fp));
    ASSERT_EQ(1, static_cast<int>(nstep.size()));
    EXPECT_GT(nstep[0], 0);
    EXPECT_FALSE(std::isnan(s[0][0]));
    EXPECT_FALSE(std::isnan(fp[0][0]));
    delete[] s[0];
    delete[] rnorm[0];
    delete[] fp[0];
}
