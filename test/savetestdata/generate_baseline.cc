#include <cstring>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
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

struct Point3 {
    double x;
    double y;
    double z;
};

void WriteFieldBaseline(const std::filesystem::path &outDir) {
    std::vector<Point3> points = {
        {1.5, 0.1, -0.2},
        {5.0, 0.0, 0.0},
        {10.0, 2.0, -1.0},
        {20.0, -10.0, 3.0},
        {30.0, 0.0, 15.0},
        {40.0, 5.0, -20.0},
    };

    std::ofstream file(outDir / "field_baseline.csv");
    file << std::setprecision(17);

    file << "x,y,z,"
         << "internal_bx,internal_by,internal_bz,"
         << "external_bx,external_by,external_bz,"
         << "combined_bx,combined_by,combined_bz\n";

    for (const auto &p : points) {
        double bix = 0.0, biy = 0.0, biz = 0.0;
        double bex = 0.0, bey = 0.0, bez = 0.0;
        double bcx = 0.0, bcy = 0.0, bcz = 0.0;

        double x = p.x;
        double y = p.y;
        double z = p.z;

        InternalField(1, &x, &y, &z, &bix, &biy, &biz);
        Con2020Field(p.x, p.y, p.z, &bex, &bey, &bez);
        ModelField(p.x, p.y, p.z, "jrm33", "Con2020", true, true, &bcx, &bcy,
                   &bcz);

        file << p.x << ',' << p.y << ',' << p.z << ',' << bix << ',' << biy << ','
             << biz << ',' << bex << ',' << bey << ',' << bez << ',' << bcx << ','
             << bcy << ',' << bcz << '\n';
    }
}

void WriteTraceBaseline(const std::filesystem::path &outDir) {
    std::vector<Point3> starts = {
        {5.0, 0.0, 0.0},
        {8.0, 3.0, 1.0},
        {12.0, -4.0, 2.0},
    };

    const int n = static_cast<int>(starts.size());
    const int maxLen = 2000;
    std::vector<double> x0(n), y0(n), z0(n);
    for (int i = 0; i < n; i++) {
        x0[i] = starts[i].x;
        y0[i] = starts[i].y;
        z0[i] = starts[i].z;
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

    bool ok = TraceField(n, x0.data(), y0.data(), z0.data(), internal, nExt,
                         externalNames, maxLen, 1.0, 0.5, 0.001, 1e-4, 0.05,
                         false, 0, 1.0, 0.93513, 0.94212, 0.94212, nstep.data(),
                         x.data(), y.data(), z.data(), bx.data(), by.data(),
                         bz.data(), r.data(), s.data(), rnorm.data(),
                         traceRegion.data(), fp.data(), 0, nullptr, nullptr);

    if (!ok) {
        std::cerr << "TraceField failed while generating baseline." << std::endl;
        std::exit(1);
    }

    std::ofstream summary(outDir / "trace_summary.csv");
    summary << std::setprecision(17);
    summary << "trace_index,start_x,start_y,start_z,";
    summary << "ion_north_x,ion_north_y,ion_north_z,";
    summary << "ion_south_x,ion_south_y,ion_south_z,";
    summary << "surf_north_x,surf_north_y,surf_north_z,";
    summary << "surf_south_x,surf_south_y,surf_south_z,";
    summary << "eq_mag_x,eq_mag_y,eq_mag_z,";
    summary << "ion_north_lat,ion_north_lon,ion_south_lat,ion_south_lon,";
    summary << "surf_north_lat,surf_north_lon,surf_south_lat,surf_south_lon,";
    summary << "eq_mlon,eq_lshell\n";

    for (size_t i = 0; i < starts.size(); i++) {
        double *traceFp = fp[i];
        summary << i << ',' << starts[i].x << ',' << starts[i].y << ',' << starts[i].z
            << ',' << traceFp[0] << ',' << traceFp[1] << ',' << traceFp[2] << ','
            << traceFp[3] << ',' << traceFp[4] << ',' << traceFp[5] << ','
            << traceFp[12] << ',' << traceFp[13] << ',' << traceFp[14] << ','
            << traceFp[15] << ',' << traceFp[16] << ',' << traceFp[17] << ','
            << traceFp[27] << ',' << traceFp[28] << ',' << traceFp[29] << ','
            << traceFp[31] << ',' << traceFp[30] << ',' << traceFp[35] << ','
            << traceFp[34] << ',' << traceFp[39] << ',' << traceFp[38] << ','
            << traceFp[43] << ',' << traceFp[42] << ',' << traceFp[47] << ','
            << traceFp[46] << '\n';
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

}  // namespace

int main() {
    constexpr double deg2rad = M_PI / 180.0;

    char model[64] = {0};
    bool cartIn = true;
    bool cartOut = true;
    int maxDeg = 0;

    GetInternalCFG(model, &cartIn, &cartOut, &maxDeg);

    (void)model;
    (void)cartIn;
    (void)cartOut;
    SetInternalCFG("jrm33", true, true, maxDeg);
    SetCon2020Params(139.6, 16.7, 7.8, 51.4, 3.6, 9.3 * deg2rad, 155.8 * deg2rad,
                     "hybrid", true, true, true, true, true, 1.0, 0.1, 417659.38364764,
                     "lmic", 0.1, 0.35, 16.1 * deg2rad, 0.5 * deg2rad,
                     10.716 * deg2rad, 0.125 * deg2rad);

    const auto outDir = std::filesystem::path("../data");
    std::filesystem::create_directories(outDir);

    WriteFieldBaseline(outDir);
    WriteTraceBaseline(outDir);

    std::cout << "Saved baseline files in test/data" << std::endl;
    return 0;
}
