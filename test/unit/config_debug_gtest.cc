#include <gtest/gtest.h>

#include <array>
#include <cmath>
#include <cstring>
#include <string>
#include <vector>

#include "regression_common.h"

extern "C" void GetCon2020Params(
    double *mui, double *irho, double *r0, double *r1, double *d, double *xt,
    double *xp, char *eqtype, bool *Edwards, bool *ErrChk, bool *CartIn,
    bool *CartOut, bool *smooth, double *DeltaRho, double *DeltaZ, double *g,
    char *azfunc, double *wO_open, double *wO_om, double *thetamm,
    double *dthetamm, double *thetaoc, double *dthetaoc);

namespace {

struct Con2020Snapshot {
    double mui;
    double irho;
    double r0;
    double r1;
    double d;
    double xt;
    double xp;
    std::string eqtype;
    bool edwards;
    bool errChk;
    bool cartIn;
    bool cartOut;
    bool smooth;
    double deltaRho;
    double deltaZ;
    double g;
    std::string azfunc;
    double wO_open;
    double wO_om;
    double thetamm;
    double dthetamm;
    double thetaoc;
    double dthetaoc;
};

Con2020Snapshot ReadCon2020Snapshot() {
    Con2020Snapshot s{};
    char eqtype[64] = {0};
    char azfunc[64] = {0};

    GetCon2020Params(&s.mui, &s.irho, &s.r0, &s.r1, &s.d, &s.xt, &s.xp,
                     eqtype, &s.edwards, &s.errChk, &s.cartIn, &s.cartOut,
                     &s.smooth, &s.deltaRho, &s.deltaZ, &s.g, azfunc,
                     &s.wO_open, &s.wO_om, &s.thetamm, &s.dthetamm,
                     &s.thetaoc, &s.dthetaoc);

    s.eqtype = eqtype;
    s.azfunc = azfunc;
    return s;
}

Con2020Snapshot ReadCon2020SnapshotViaJupitermag() {
    Con2020Snapshot s{};
    char eqtype[64] = {0};
    char azfunc[64] = {0};

    JupitermagGetCon2020Params(&s.mui, &s.irho, &s.r0, &s.r1, &s.d, &s.xt,
                               &s.xp, eqtype, &s.edwards, &s.errChk,
                               &s.cartIn, &s.cartOut, &s.smooth, &s.deltaRho,
                               &s.deltaZ, &s.g, azfunc, &s.wO_open, &s.wO_om,
                               &s.thetamm, &s.dthetamm, &s.thetaoc,
                               &s.dthetaoc);

    s.eqtype = eqtype;
    s.azfunc = azfunc;
    return s;
}

void ExpectNear(double a, double b, const char *name) {
    EXPECT_NEAR(a, b, 1e-12) << name;
}

void ExpectSameExceptCartFlags(const Con2020Snapshot &expectedBase,
                               const Con2020Snapshot &actual,
                               bool expectedCartIn,
                               bool expectedCartOut) {
    ExpectNear(expectedBase.mui, actual.mui, "mui");
    ExpectNear(expectedBase.irho, actual.irho, "irho");
    ExpectNear(expectedBase.r0, actual.r0, "r0");
    ExpectNear(expectedBase.r1, actual.r1, "r1");
    ExpectNear(expectedBase.d, actual.d, "d");
    ExpectNear(expectedBase.xt, actual.xt, "xt");
    ExpectNear(expectedBase.xp, actual.xp, "xp");
    EXPECT_EQ(expectedBase.eqtype, actual.eqtype) << "eqtype";
    EXPECT_EQ(expectedBase.edwards, actual.edwards) << "edwards";
    EXPECT_EQ(expectedBase.errChk, actual.errChk) << "errChk";
    EXPECT_EQ(expectedCartIn, actual.cartIn) << "cartIn";
    EXPECT_EQ(expectedCartOut, actual.cartOut) << "cartOut";
    EXPECT_EQ(expectedBase.smooth, actual.smooth) << "smooth";
    ExpectNear(expectedBase.deltaRho, actual.deltaRho, "deltaRho");
    ExpectNear(expectedBase.deltaZ, actual.deltaZ, "deltaZ");
    ExpectNear(expectedBase.g, actual.g, "g");
    EXPECT_EQ(expectedBase.azfunc, actual.azfunc) << "azfunc";
    ExpectNear(expectedBase.wO_open, actual.wO_open, "wO_open");
    ExpectNear(expectedBase.wO_om, actual.wO_om, "wO_om");
    ExpectNear(expectedBase.thetamm, actual.thetamm, "thetamm");
    ExpectNear(expectedBase.dthetamm, actual.dthetamm, "dthetamm");
    ExpectNear(expectedBase.thetaoc, actual.thetaoc, "thetaoc");
    ExpectNear(expectedBase.dthetaoc, actual.dthetaoc, "dthetaoc");
}

}  // namespace

TEST(ConfigDebug, ModelFieldConfigRoundTrip) {
    ConfigureModelsForBaseline();

    const Con2020Snapshot base = ReadCon2020SnapshotViaJupitermag();
    const std::array<std::array<bool, 2>, 5> cartModes = {
        {{true, true}, {false, false}, {true, false}, {false, true}, {true, true}}};

    for (size_t i = 0; i < cartModes.size(); i++) {
        const bool cartIn = cartModes[i][0];
        const bool cartOut = cartModes[i][1];

        double b0 = 0.0, b1 = 0.0, b2 = 0.0;
        ModelField(5.0, 0.0, 0.0, "jrm33", "Con2020", cartIn, cartOut, &b0, &b1,
                   &b2);

        const Con2020Snapshot after = ReadCon2020SnapshotViaJupitermag();
        SCOPED_TRACE(i);
        ExpectSameExceptCartFlags(base, after, cartIn, cartOut);
    }
}

TEST(ConfigDebug, ModelFieldArrayConfigRoundTrip) {
    ConfigureModelsForBaseline();

    const Con2020Snapshot base = ReadCon2020SnapshotViaJupitermag();
    const std::array<std::array<bool, 2>, 4> cartModes = {
        {{true, true}, {false, true}, {true, false}, {false, false}}};

    const int n = 2;
    double p0[n] = {5.0, 10.0};
    double p1[n] = {0.0, 2.0};
    double p2[n] = {0.0, -1.0};
    double b0[n] = {0.0, 0.0};
    double b1[n] = {0.0, 0.0};
    double b2[n] = {0.0, 0.0};

    for (size_t i = 0; i < cartModes.size(); i++) {
        const bool cartIn = cartModes[i][0];
        const bool cartOut = cartModes[i][1];

        ModelFieldArray(n, p0, p1, p2, "jrm33", "Con2020", cartIn, cartOut, b0,
                        b1, b2);

        const Con2020Snapshot after = ReadCon2020SnapshotViaJupitermag();
        SCOPED_TRACE(i);
        ExpectSameExceptCartFlags(base, after, cartIn, cartOut);
    }
}

TEST(ConfigDebug, DirectApiMatchesJupitermagApi) {
    ConfigureModelsForBaseline();

    double b0 = 0.0, b1 = 0.0, b2 = 0.0;
    ModelField(5.0, 0.0, 0.0, "jrm33", "Con2020", false, false, &b0, &b1, &b2);

    const Con2020Snapshot direct = ReadCon2020Snapshot();
    const Con2020Snapshot viaJupitermag = ReadCon2020SnapshotViaJupitermag();

    EXPECT_EQ(viaJupitermag.cartIn, direct.cartIn) << "cartIn mismatch between APIs";
    EXPECT_EQ(viaJupitermag.cartOut, direct.cartOut) << "cartOut mismatch between APIs";
}
