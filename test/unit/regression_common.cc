#include "regression_common.h"

#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace {

std::vector<std::string> SplitCSV(const std::string &line) {
    std::vector<std::string> out;
    std::stringstream ss(line);
    std::string item;
    while (std::getline(ss, item, ',')) {
        out.push_back(item);
    }
    return out;
}

}  // namespace

bool NearlyEqual(double expected, double actual, double absTol, double relTol) {
    const double diff = std::fabs(expected - actual);
    if (diff <= absTol) {
        return true;
    }
    const double scale = std::max(std::fabs(expected), std::fabs(actual));
    return diff <= (relTol * scale);
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

std::vector<std::string> ReadCSVHeader(const std::filesystem::path &file) {
    std::ifstream in(file);
    if (!in) {
        throw std::runtime_error("Could not open: " + file.string());
    }
    std::string header;
    std::getline(in, header);
    if (header.empty()) {
        throw std::runtime_error("Empty header in: " + file.string());
    }
    return SplitCSV(header);
}

void ConfigureModelsForBaseline() {
    char model[64] = {0};
    bool cartIn = true;
    bool cartOut = true;
    int maxDeg = 0;

    JupiterMagGetInternalCFG(model, &cartIn, &cartOut, &maxDeg);
    JupiterMagSetInternalCFG("jrm33", true, true, maxDeg);

    JupiterMagSetCon2020Params(139.6, 16.7, 7.8, 51.4, 3.6, 9.3 * kDeg2Rad,
                               155.8 * kDeg2Rad, "hybrid", true, true, true,
                               true, true, 1.0, 0.1, 417659.38364764, "lmic",
                               0.1, 0.35, 16.1 * kDeg2Rad, 0.5 * kDeg2Rad,
                               10.716 * kDeg2Rad, 0.125 * kDeg2Rad);
}

std::filesystem::path TestDataDir() {
    const auto cwd = std::filesystem::current_path();
    return cwd / "data";
}

std::vector<std::array<double, 3>> RegressionPoints() {
    return {
        {1.5, 0.1, -0.2},
        {5.0, 0.0, 0.0},
        {10.0, 2.0, -1.0},
        {20.0, -10.0, 3.0},
        {30.0, 0.0, 15.0},
        {40.0, 5.0, -20.0},
    };
}

bool RunTraceWithDir(const std::vector<std::array<double, 3>> &starts, int traceDir,
                     std::vector<int> *nstepOut, std::vector<double *> *sOut,
                     std::vector<double *> *rnormOut,
                     std::vector<double *> *fpOut) {
    const int n = static_cast<int>(starts.size());
    const int maxLen = 1000;

    std::vector<double> x0(n), y0(n), z0(n);
    for (int i = 0; i < n; i++) {
        x0[i] = starts[i][0];
        y0[i] = starts[i][1];
        z0[i] = starts[i][2];
    }

    const char *internal = "jrm33";
    int nExt = 1;
    char extName[] = "Con2020";
    char *externalNames[1] = {extName};

    std::vector<int> nstep(n);
    std::vector<double *> x(n), y(n), z(n), bx(n), by(n), bz(n), r(n), s(n),
        rnorm(n), fp(n);
    std::vector<int *> traceRegion(n);

    const double sentinel = std::numeric_limits<double>::quiet_NaN();

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

        for (int j = 0; j < maxLen; j++) {
            s[i][j] = sentinel;
            rnorm[i][j] = sentinel;
        }
        for (int j = 0; j < 49; j++) {
            fp[i][j] = sentinel;
        }
    }

    const bool ok = TraceField(n, x0.data(), y0.data(), z0.data(), internal, nExt,
                               externalNames, maxLen, 1.0, 0.5, 0.001, 1e-4,
                               0.05, false, traceDir, 1.0, 0.93513, 0.94212,
                               0.94212, nstep.data(), x.data(), y.data(),
                               z.data(), bx.data(), by.data(), bz.data(),
                               r.data(), s.data(), rnorm.data(), traceRegion.data(),
                               fp.data(), 0, nullptr, nullptr);

    if (nstepOut != nullptr) {
        *nstepOut = nstep;
    }
    if (sOut != nullptr) {
        *sOut = s;
    }
    if (rnormOut != nullptr) {
        *rnormOut = rnorm;
    }
    if (fpOut != nullptr) {
        *fpOut = fp;
    }

    for (int i = 0; i < n; i++) {
        delete[] x[i];
        delete[] y[i];
        delete[] z[i];
        delete[] bx[i];
        delete[] by[i];
        delete[] bz[i];
        delete[] r[i];
        delete[] traceRegion[i];
    }

    if (sOut == nullptr) {
        for (int i = 0; i < n; i++) {
            delete[] s[i];
        }
    }
    if (rnormOut == nullptr) {
        for (int i = 0; i < n; i++) {
            delete[] rnorm[i];
        }
    }
    if (fpOut == nullptr) {
        for (int i = 0; i < n; i++) {
            delete[] fp[i];
        }
    }

    return ok;
}

std::vector<Footprint49> RunTraceAndGetFootprints(
    const std::vector<std::array<double, 3>> &starts) {
    const int n = static_cast<int>(starts.size());
    const int maxLen = 2000;

    std::vector<double> x0(n), y0(n), z0(n);
    for (int i = 0; i < n; i++) {
        x0[i] = starts[i][0];
        y0[i] = starts[i][1];
        z0[i] = starts[i][2];
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

    if (!ok) {
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
        throw std::runtime_error("TraceField failed in RunTraceAndGetFootprints");
    }

    std::vector<Footprint49> out(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < 49; j++) {
            out[i][j] = fp[i][j];
        }
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
    return out;
}
