#include <gtest/gtest.h>

#include <cmath>
#include <string>
#include <vector>

#include "regression_common.h"

TEST(Regressions, DeterministicFieldAndTraceFootprints) {
    ConfigureModelsForBaseline();

    const auto points = RegressionPoints();
    for (const auto &p : points) {
        double b0a = 0.0, b1a = 0.0, b2a = 0.0;
        double b0b = 0.0, b1b = 0.0, b2b = 0.0;
        ModelField(p[0], p[1], p[2], "jrm33", "Con2020", true, true, &b0a, &b1a,
                   &b2a);
        ModelField(p[0], p[1], p[2], "jrm33", "Con2020", true, true, &b0b, &b1b,
                   &b2b);
        EXPECT_TRUE(NearlyEqual(b0a, b0b, kFieldAbsTol, kFieldRelTol));
        EXPECT_TRUE(NearlyEqual(b1a, b1b, kFieldAbsTol, kFieldRelTol));
        EXPECT_TRUE(NearlyEqual(b2a, b2b, kFieldAbsTol, kFieldRelTol));
    }

    const std::vector<std::array<double, 3>> starts = {{5.0, 0.0, 0.0}};
    const auto fpA = RunTraceAndGetFootprints(starts);
    const auto fpB = RunTraceAndGetFootprints(starts);
    ASSERT_EQ(fpA.size(), fpB.size());
    for (size_t i = 0; i < fpA.size(); i++) {
        for (int j = 0; j < 49; j++) {
            EXPECT_TRUE(NearlyEqual(fpA[i][j], fpB[i][j], kTraceAbsTol,
                                    kTraceRelTol))
                << "row=" << i << " col=" << j;
        }
    }
}

TEST(Regressions, EdgeBehaviorAndInvalidModelPointer) {
    ConfigureModelsForBaseline();

    const modelFieldPtr invalid = getModelFieldPtr("definitely_not_a_model");
    EXPECT_EQ(nullptr, invalid);

    const std::vector<std::array<double, 3>> edgePoints = {
        {0.0, 0.0, 0.0},
        {1e-6, -1e-6, 1e-6},
        {100.0, -100.0, 30.0},
    };
    for (size_t i = 0; i < edgePoints.size(); i++) {
        const auto &p = edgePoints[i];
        double b0 = 0.0, b1 = 0.0, b2 = 0.0;
        ModelField(p[0], p[1], p[2], "none", "none", true, true, &b0, &b1, &b2);
        EXPECT_TRUE(NearlyEqual(0.0, b0, kFieldAbsTol, kFieldRelTol));
        EXPECT_TRUE(NearlyEqual(0.0, b1, kFieldAbsTol, kFieldRelTol));
        EXPECT_TRUE(NearlyEqual(0.0, b2, kFieldAbsTol, kFieldRelTol));

        ModelField(p[0], p[1], p[2], "jrm33", "Con2020", true, true, &b0, &b1,
                   &b2);
        EXPECT_FALSE(std::isinf(b0));
        EXPECT_FALSE(std::isinf(b1));
        EXPECT_FALSE(std::isinf(b2));

        if (i == edgePoints.size() - 1) {
            EXPECT_TRUE(std::isfinite(b0));
            EXPECT_TRUE(std::isfinite(b1));
            EXPECT_TRUE(std::isfinite(b2));
        }
    }
}

TEST(Regressions, AllModelsSmokeAtSafePoint) {
    ConfigureModelsForBaseline();

    const std::vector<std::string> modelNames = {
        "jrm33", "jrm09", "vip4", "vipal", "vit4", "o4", "o6",
        "gsfc15ev", "gsfc15evs", "gsfc13ev", "u17ev", "sha", "v117ev",
        "jpl15ev", "jpl15evs", "p11a", "isaac"};

    const double x = 10.0;
    const double y = 0.0;
    const double z = 0.0;

    for (const auto &name : modelNames) {
        modelFieldPtr ptr = getModelFieldPtr(name.c_str());
        ASSERT_NE(nullptr, ptr) << name;

        double bx = 0.0, by = 0.0, bz = 0.0;
        ptr(x, y, z, &bx, &by, &bz);
        EXPECT_FALSE(std::isinf(bx)) << name;
        EXPECT_FALSE(std::isinf(by)) << name;
        EXPECT_FALSE(std::isinf(bz)) << name;
    }
}
