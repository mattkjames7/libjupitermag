#include <gtest/gtest.h>

#include <cstddef>
#include <vector>

#include "regression_common.h"

TEST(Regressions, ScalarArrayConsistency) {
    ConfigureModelsForBaseline();

    const auto points = RegressionPoints();
    const int n = static_cast<int>(points.size());
    std::vector<double> p0(n), p1(n), p2(n);
    std::vector<double> iArr0(n), iArr1(n), iArr2(n);
    std::vector<double> eArr0(n), eArr1(n), eArr2(n);
    std::vector<double> mArr0(n), mArr1(n), mArr2(n);

    for (int i = 0; i < n; i++) {
        p0[i] = points[i][0];
        p1[i] = points[i][1];
        p2[i] = points[i][2];
    }

    InternalField(n, p0.data(), p1.data(), p2.data(), iArr0.data(), iArr1.data(),
                  iArr2.data());
    Con2020FieldArray(n, p0.data(), p1.data(), p2.data(), eArr0.data(),
                      eArr1.data(), eArr2.data());
    ModelFieldArray(n, p0.data(), p1.data(), p2.data(), "jrm33", "Con2020",
                    true, true, mArr0.data(), mArr1.data(), mArr2.data());

    modelFieldPtr internalScalar = getModelFieldPtr("jrm33");
    ASSERT_NE(nullptr, internalScalar);

    for (int i = 0; i < n; i++) {
        double is0 = 0.0, is1 = 0.0, is2 = 0.0;
        double es0 = 0.0, es1 = 0.0, es2 = 0.0;
        double ms0 = 0.0, ms1 = 0.0, ms2 = 0.0;

        internalScalar(p0[i], p1[i], p2[i], &is0, &is1, &is2);
        Con2020Field(p0[i], p1[i], p2[i], &es0, &es1, &es2);
        ModelField(p0[i], p1[i], p2[i], "jrm33", "Con2020", true, true, &ms0,
                   &ms1, &ms2);

        EXPECT_TRUE(NearlyEqual(is0, iArr0[i], kFieldAbsTol, kFieldRelTol)) << i;
        EXPECT_TRUE(NearlyEqual(is1, iArr1[i], kFieldAbsTol, kFieldRelTol)) << i;
        EXPECT_TRUE(NearlyEqual(is2, iArr2[i], kFieldAbsTol, kFieldRelTol)) << i;

        EXPECT_TRUE(NearlyEqual(es0, eArr0[i], kFieldAbsTol, kFieldRelTol)) << i;
        EXPECT_TRUE(NearlyEqual(es1, eArr1[i], kFieldAbsTol, kFieldRelTol)) << i;
        EXPECT_TRUE(NearlyEqual(es2, eArr2[i], kFieldAbsTol, kFieldRelTol)) << i;

        EXPECT_TRUE(NearlyEqual(ms0, mArr0[i], kFieldAbsTol, kFieldRelTol)) << i;
        EXPECT_TRUE(NearlyEqual(ms1, mArr1[i], kFieldAbsTol, kFieldRelTol)) << i;
        EXPECT_TRUE(NearlyEqual(ms2, mArr2[i], kFieldAbsTol, kFieldRelTol)) << i;
    }
}

TEST(Regressions, ModelCompositionIdentity) {
    ConfigureModelsForBaseline();

    const auto points = RegressionPoints();
    for (size_t i = 0; i < points.size(); i++) {
        const double x = points[i][0];
        const double y = points[i][1];
        const double z = points[i][2];

        double bi0 = 0.0, bi1 = 0.0, bi2 = 0.0;
        double be0 = 0.0, be1 = 0.0, be2 = 0.0;
        double bc0 = 0.0, bc1 = 0.0, bc2 = 0.0;
        double boi0 = 0.0, boi1 = 0.0, boi2 = 0.0;
        double boe0 = 0.0, boe1 = 0.0, boe2 = 0.0;

        double xx = x, yy = y, zz = z;
        InternalField(1, &xx, &yy, &zz, &bi0, &bi1, &bi2);
        Con2020Field(x, y, z, &be0, &be1, &be2);
        ModelField(x, y, z, "jrm33", "Con2020", true, true, &bc0, &bc1, &bc2);

        ModelField(x, y, z, "jrm33", "none", true, true, &boi0, &boi1, &boi2);
        ModelField(x, y, z, "none", "Con2020", true, true, &boe0, &boe1,
                   &boe2);

        EXPECT_TRUE(
            NearlyEqual(bi0 + be0, bc0, 10.0 * kFieldAbsTol, 10.0 * kFieldRelTol))
            << i;
        EXPECT_TRUE(
            NearlyEqual(bi1 + be1, bc1, 10.0 * kFieldAbsTol, 10.0 * kFieldRelTol))
            << i;
        EXPECT_TRUE(
            NearlyEqual(bi2 + be2, bc2, 10.0 * kFieldAbsTol, 10.0 * kFieldRelTol))
            << i;

        EXPECT_TRUE(NearlyEqual(bi0, boi0, kFieldAbsTol, kFieldRelTol)) << i;
        EXPECT_TRUE(NearlyEqual(bi1, boi1, kFieldAbsTol, kFieldRelTol)) << i;
        EXPECT_TRUE(NearlyEqual(bi2, boi2, kFieldAbsTol, kFieldRelTol)) << i;

        EXPECT_TRUE(NearlyEqual(be0, boe0, kFieldAbsTol, kFieldRelTol)) << i;
        EXPECT_TRUE(NearlyEqual(be1, boe1, kFieldAbsTol, kFieldRelTol)) << i;
        EXPECT_TRUE(NearlyEqual(be2, boe2, kFieldAbsTol, kFieldRelTol)) << i;
    }
}

TEST(Regressions, FieldBaselineCSV) {
    ConfigureModelsForBaseline();

    const auto file = TestDataDir() / "field_baseline.csv";
    const auto rows = ReadNumericCSV(file, 12);
    ASSERT_FALSE(rows.empty());

    for (size_t i = 0; i < rows.size(); i++) {
        const auto &r = rows[i];

        double x = r[0];
        double y = r[1];
        double z = r[2];

        double bix = 0.0, biy = 0.0, biz = 0.0;
        double bex = 0.0, bey = 0.0, bez = 0.0;
        double bcx = 0.0, bcy = 0.0, bcz = 0.0;

        InternalField(1, &x, &y, &z, &bix, &biy, &biz);
        Con2020Field(x, y, z, &bex, &bey, &bez);
        ModelField(x, y, z, "jrm33", "Con2020", true, true, &bcx, &bcy, &bcz);

        EXPECT_TRUE(NearlyEqual(r[3], bix, kFieldAbsTol, kFieldRelTol)) << i;
        EXPECT_TRUE(NearlyEqual(r[4], biy, kFieldAbsTol, kFieldRelTol)) << i;
        EXPECT_TRUE(NearlyEqual(r[5], biz, kFieldAbsTol, kFieldRelTol)) << i;

        EXPECT_TRUE(NearlyEqual(r[6], bex, kFieldAbsTol, kFieldRelTol)) << i;
        EXPECT_TRUE(NearlyEqual(r[7], bey, kFieldAbsTol, kFieldRelTol)) << i;
        EXPECT_TRUE(NearlyEqual(r[8], bez, kFieldAbsTol, kFieldRelTol)) << i;

        EXPECT_TRUE(NearlyEqual(r[9], bcx, kFieldAbsTol, kFieldRelTol)) << i;
        EXPECT_TRUE(NearlyEqual(r[10], bcy, kFieldAbsTol, kFieldRelTol)) << i;
        EXPECT_TRUE(NearlyEqual(r[11], bcz, kFieldAbsTol, kFieldRelTol)) << i;
    }
}

TEST(Regressions, InputArraysNotModified) {
    ConfigureModelsForBaseline();

    const auto points = RegressionPoints();
    const int n = static_cast<int>(points.size());

    std::vector<double> p0(n), p1(n), p2(n);
    for (int i = 0; i < n; i++) {
        p0[i] = points[i][0];
        p1[i] = points[i][1];
        p2[i] = points[i][2];
    }

    const auto p0Orig = p0;
    const auto p1Orig = p1;
    const auto p2Orig = p2;

    std::vector<double> o0(n), o1(n), o2(n);
    InternalField(n, p0.data(), p1.data(), p2.data(), o0.data(), o1.data(),
                  o2.data());
    EXPECT_EQ(p0Orig, p0);
    EXPECT_EQ(p1Orig, p1);
    EXPECT_EQ(p2Orig, p2);

    Con2020FieldArray(n, p0.data(), p1.data(), p2.data(), o0.data(), o1.data(),
                      o2.data());
    EXPECT_EQ(p0Orig, p0);
    EXPECT_EQ(p1Orig, p1);
    EXPECT_EQ(p2Orig, p2);

    ModelFieldArray(n, p0.data(), p1.data(), p2.data(), "jrm33", "Con2020",
                    true, true, o0.data(), o1.data(), o2.data());
    EXPECT_EQ(p0Orig, p0);
    EXPECT_EQ(p1Orig, p1);
    EXPECT_EQ(p2Orig, p2);
}
