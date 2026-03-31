#include <gtest/gtest.h>

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

extern "C" {

void InternalField(int n, double *p0, double *p1, double *p2, double *B0,
                   double *B1, double *B2);
void SetInternalCFG(const char *Model, bool CartIn, bool CartOut, int MaxDeg);
void GetInternalCFG(char *Model, bool *CartIn, bool *CartOut, int *MaxDeg);

void Con2020Field(double p0, double p1, double p2, double *B0, double *B1,
                  double *B2);
void SetCon2020Params(double mui, double irho, double r0, double r1, double d,
                      double xt, double xp, const char *eqtype, bool Edwards,
                      bool ErrChk, bool CartIn, bool CartOut, bool smooth,
                      double DeltaRho, double DeltaZ, double g,
                      const char *azfunc, double wO_open, double wO_om,
                      double thetamm, double dthetamm, double thetaoc,
                      double dthetaoc);

void ModelField(double p0, double p1, double p2, const char *internal,
                const char *external, bool CartIn, bool CartOut, double *B0,
                double *B1, double *B2);

bool TraceField(int n, double *x0, double *y0, double *z0, const char *IntFunc,
                int nExt, char **ExtFunc, int MaxLen, double MaxStep,
                double InitStep, double MinStep, double ErrMax, double Delta,
                bool Verbose, int TraceDir, double as, double bs, double ai,
                double bi, int *nstep, double **x, double **y, double **z,
                double **Bx, double **By, double **Bz, double **R, double **S,
                double **Rnorm, int **traceRegion, double **FP, int nalpha,
                double *alpha, double *halpha);
}

namespace {

constexpr double kDeg2Rad = M_PI / 180.0;
constexpr double kFieldAbsTol = 1e-6;
constexpr double kFieldRelTol = 1e-6;
constexpr double kTraceAbsTol = 1e-6;
constexpr double kTraceRelTol = 1e-6;

bool NearlyEqual(double expected, double actual, double absTol, double relTol) {
    const double diff = std::fabs(expected - actual);
    if (diff <= absTol) {
        return true;
    }
    const double scale = std::max(std::fabs(expected), std::fabs(actual));
    return diff <= (relTol * scale);
}

std::vector<std::string> SplitCSV(const std::string &line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, ',')) {
        out.push_back(item);
    }
    return out;
}

std::vector<std::vector<double>> ReadNumericCSV(const std::filesystem::path &file,
                                                size_t expectedColumns) {
    std::ifstream in(file);
    if (!in) {
        throw std::runtime_error("Could not open: " + file.string());
    }

    std::string line;
    std::getline(in, line);  // header

    std::vector<std::vector<double>> rows;
    while (std::getline(in, line)) {
        if (line.empty()) {
            continue;
        }
        std::vector<std::string> fields = SplitCSV(line);
        if (fields.size() != expectedColumns) {
            throw std::runtime_error("Unexpected column count in: " + file.string());
        }
        std::vector<double> row;
        row.reserve(fields.size());
        for (const auto &f : fields) {
            row.push_back(std::stod(f));
        }
        rows.push_back(std::move(row));
    }
    return rows;
}

void ConfigureModelsForBaseline() {
    char model[64] = {0};
    bool cartIn = true;
    bool cartOut = true;
    int maxDeg = 0;

    GetInternalCFG(model, &cartIn, &cartOut, &maxDeg);
    SetInternalCFG("jrm33", true, true, maxDeg);

    SetCon2020Params(139.6, 16.7, 7.8, 51.4, 3.6, 9.3 * kDeg2Rad,
                     155.8 * kDeg2Rad, "hybrid", true, true, true, true, true,
                     1.0, 0.1, 417659.38364764, "lmic", 0.1, 0.35,
                     16.1 * kDeg2Rad, 0.5 * kDeg2Rad, 10.716 * kDeg2Rad,
                     0.125 * kDeg2Rad);
}

std::filesystem::path TestDataDir() {
    const auto cwd = std::filesystem::current_path();
    return cwd / "data";
}

}  // namespace

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

TEST(Regressions, TraceFootprintBaselineCSV) {
    ConfigureModelsForBaseline();

    const auto file = TestDataDir() / "trace_summary.csv";
    const auto rows = ReadNumericCSV(file, 29);
    ASSERT_FALSE(rows.empty());

    const int n = static_cast<int>(rows.size());
    const int maxLen = 2000;

    std::vector<double> x0(n), y0(n), z0(n);
    for (int i = 0; i < n; i++) {
        x0[i] = rows[i][1];
        y0[i] = rows[i][2];
        z0[i] = rows[i][3];
    }

    const char *internal = "jrm33";
    int nExt = 1;
    char extName[] = "Con2020";
    char *externalNames[1] = {extName};

    std::vector<int> nstep(n);
    std::vector<double *> x(n), y(n), z(n), bx(n), by(n), bz(n), r(n), s(n),
        rnorm(n), fp(n);
    std::vector<int *> traceRegion(n);

    for (int i = 0; i < n; i++) {
        x[i] = new double[maxLen];
        y[i] = new double[maxLen];
        z[i] = new double[maxLen];
        bx[i] = new double[maxLen];
        by[i] = new double[maxLen];
        bz[i] = new double[maxLen];
        r[i] = new double[maxLen];
        s[i] = new double[maxLen];
        rnorm[i] = new double[maxLen];
        traceRegion[i] = new int[maxLen];
        fp[i] = new double[49];
    }

    const bool ok = TraceField(n, x0.data(), y0.data(), z0.data(), internal, nExt,
                               externalNames, maxLen, 1.0, 0.5, 0.001, 1e-4,
                               0.05, false, 0, 1.0, 0.93513, 0.94212, 0.94212,
                               nstep.data(), x.data(), y.data(), z.data(),
                               bx.data(), by.data(), bz.data(), r.data(),
                               s.data(), rnorm.data(), traceRegion.data(),
                               fp.data(), 0, nullptr, nullptr);

    ASSERT_TRUE(ok);

    for (int i = 0; i < n; i++) {
        const auto &row = rows[i];
        const double *traceFp = fp[i];

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

    for (int i = 0; i < n; i++) {
        delete[] x[i];
        delete[] y[i];
        delete[] z[i];
        delete[] bx[i];
        delete[] by[i];
        delete[] bz[i];
        delete[] r[i];
        delete[] s[i];
        delete[] rnorm[i];
        delete[] traceRegion[i];
        delete[] fp[i];
    }
}
