#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <vector>

#include "regression_common.h"

TEST(Regressions, CoordinateRoundTripSIIIAndMag) {
    const std::vector<std::array<double, 3>> vectors = {
        {1.0, 0.0, 0.0},
        {0.0, 1.0, 0.0},
        {0.0, 0.0, 1.0},
        {3.0, -2.0, 5.0},
        {-4.5, 1.3, 2.2},
    };
    const std::vector<std::array<double, 2>> angles = {
        {0.0, 0.0},
        {10.25 * kDeg2Rad, 163.62 * kDeg2Rad},
        {5.0 * kDeg2Rad, 45.0 * kDeg2Rad},
        {15.0 * kDeg2Rad, 270.0 * kDeg2Rad},
    };

    for (const auto &v : vectors) {
        for (const auto &a : angles) {
            double xm = 0.0, ym = 0.0, zm = 0.0;
            double x3 = 0.0, y3 = 0.0, z3 = 0.0;
            SIIItoMag(v[0], v[1], v[2], a[0], a[1], &xm, &ym, &zm);
            MagtoSIII(xm, ym, zm, a[0], a[1], &x3, &y3, &z3);

            EXPECT_TRUE(NearlyEqual(v[0], x3, kGeomAbsTol, kGeomRelTol));
            EXPECT_TRUE(NearlyEqual(v[1], y3, kGeomAbsTol, kGeomRelTol));
            EXPECT_TRUE(NearlyEqual(v[2], z3, kGeomAbsTol, kGeomRelTol));
        }
    }
}

TEST(Regressions, FootprintGeometryAndRanges) {
    ConfigureModelsForBaseline();

    const auto rows = ReadNumericCSV(TestDataDir() / "trace_summary.csv", 29);
    ASSERT_FALSE(rows.empty());

    const double as = 1.0;
    const double bs = 0.93513;
    const double ai = 0.94212;
    const double bi = 0.94212;

    for (size_t i = 0; i < rows.size(); i++) {
        const auto &r = rows[i];

        const double rin = std::sqrt(r[4] * r[4] + r[5] * r[5] + r[6] * r[6]);
        const double ris = std::sqrt(r[7] * r[7] + r[8] * r[8] + r[9] * r[9]);
        const double rsn = std::sqrt(r[10] * r[10] + r[11] * r[11] + r[12] * r[12]);
        const double rss = std::sqrt(r[13] * r[13] + r[14] * r[14] + r[15] * r[15]);

        EXPECT_TRUE(NearlyEqual(ai, rin, 2e-3, 2e-3)) << i;
        EXPECT_TRUE(NearlyEqual(bi, ris, 2e-3, 2e-3)) << i;

        const double sNorth =
            ((r[10] * r[10] + r[11] * r[11]) / (as * as)) +
            ((r[12] * r[12]) / (bs * bs));
        const double sSouth =
            ((r[13] * r[13] + r[14] * r[14]) / (as * as)) +
            ((r[15] * r[15]) / (bs * bs));
        EXPECT_TRUE(NearlyEqual(1.0, sNorth, 5e-3, 5e-3)) << i;
        EXPECT_TRUE(NearlyEqual(1.0, sSouth, 5e-3, 5e-3)) << i;

        for (int j = 19; j <= 26; j++) {
            EXPECT_TRUE(std::isfinite(r[j])) << "row=" << i << " col=" << j;
        }
        EXPECT_GE(r[19], -90.0);
        EXPECT_LE(r[19], 90.0);
        EXPECT_GE(r[21], -90.0);
        EXPECT_LE(r[21], 90.0);
        EXPECT_GE(r[23], -90.0);
        EXPECT_LE(r[23], 90.0);
        EXPECT_GE(r[25], -90.0);
        EXPECT_LE(r[25], 90.0);

        EXPECT_GE(r[20], -180.0);
        EXPECT_LE(r[20], 180.0);
        EXPECT_GE(r[22], -180.0);
        EXPECT_LE(r[22], 180.0);
        EXPECT_GE(r[24], -180.0);
        EXPECT_LE(r[24], 180.0);
        EXPECT_GE(r[26], -180.0);
        EXPECT_LE(r[26], 180.0);
    }
}
