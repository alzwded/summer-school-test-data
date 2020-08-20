// compile with:
//    cl.exe /EHa /O2 /openmp /std:c++17 genplate.cpp
//    gcc -fopenmp -o genplate -O3 genplate.cpp -lm
#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <map>
#include <vector>
#include <array>
#include <omp.h>

struct PT
{
    double x, y, z;
};

typedef std::array<int64_t, 8> HE;

struct Mesh
{
    std::vector<PT> nodes;
    std::vector<int64_t> nodeIds;
    std::vector<HE> connect;
};

struct Mode
{
    int64_t number;
    std::vector<PT> displacement;
};

Mesh genMesh(int64_t N, double L)
{
    FILE*g = fopen("debug2.txt", "w");
    Mesh rval;
    rval.nodeIds.resize(N * N * N, -1);
    for(int64_t k = 0; k < N; ++k) {
      printf("...generating 1/2  %f%%\n", 100.0f*(float)k/N);
      for(int64_t i = 0; i < N; ++i) {
        for(int64_t j = 0; j < N; ++j) {
            // poke a hole in the mesh if it's big enough:-)
            if((N >= 16
                &&! (i >= 4.0 / 8 * N && i < 5.0 / 8 * N &&
                     j >= 4.0 / 8 * N && j < 7.0 / 8 * N))
            || N < 16)
            {
                PT p = { -L/2.0 + L / (N - 1) * i,
                         -L/2.0 + L / (N - 1) * j,
                         -L/2.0 + L / (N - 1) * k };
                rval.nodeIds[k * N * N + i * N + j] = rval.nodes.size();
                rval.nodes.push_back(p);
            }
            if(i < N - 1 && j < N - 1 && k < N - 1) {
                rval.connect.push_back(HE{
                    k * N * N + i * N + j,
                    k * N * N + i * N + j + 1,
                    k * N * N + (i + 1) * N + j + 1,
                    k * N * N + (i + 1) * N + j,
                    (k + 1) * N * N + i * N + j,
                    (k + 1) * N * N + i * N + j + 1,
                    (k + 1) * N * N + (i + 1) * N + j + 1,
                    (k + 1) * N * N + (i + 1) * N + j,
                });
                fprintf(g, "%lld %lld %lld %lld %lld %lld %lld %lld\n",
                    rval.connect.back()[0],
                    rval.connect.back()[1],
                    rval.connect.back()[2],
                    rval.connect.back()[3],
                    rval.connect.back()[4],
                    rval.connect.back()[5],
                    rval.connect.back()[6],
                    rval.connect.back()[7]);
            }
        }
      }
    }
    FILE* f = fopen("debug.txt", "w");
    for(int64_t i = 0; i < N*N*N; ++i) {
        fprintf(f, "%lld*%lld*%lld=%lld %lld\n", i / N / N, i / N, i % N, i, rval.nodeIds[i]);
    }
    decltype(rval.connect) newConnect;
    for(int64_t i = 0; i < rval.connect.size(); ++i) {
        if(i % (N*N) == 0) printf("...generating 2/2  %f%%\n", 100.0f*(float)i/rval.connect.size());
        bool erase = false;
        for(int64_t j = 0; j < rval.connect[i].size(); ++j) {
            fprintf(g, "%lld => %lld\n", rval.connect[i][j], rval.nodeIds[rval.connect[i][j]]);
            rval.connect[i][j] = rval.nodeIds[rval.connect[i][j]];
            if(rval.connect[i][j] == -1) erase = true;
        }
        if(!erase) newConnect.push_back(rval.connect[i]);
    }
    std::swap(rval.connect, newConnect);
    fclose(f); fclose(g);
    return rval;
}

void output_csv(Mesh const& mesh) {}
void output_tec(Mesh const& mesh) {}
void output_txt(Mesh const& mesh)
{
    FILE* f;
    f = fopen("nodes.txt", "w");
    int cc, cn;
    cc = 0;
    cn = mesh.nodes.size();
    for(auto&& n : mesh.nodes) {
        if((cc++ % (cn / 10)) == 0) printf("...writing nodes %f%%\n", 100.0f * cc / cn);
        fprintf(f, "%lf %lf %lf\n", n.x, n.y, n.z);
    }
    fclose(f);
    f = fopen("connect.txt", "w");
    cc = 0;
    cn = mesh.connect.size();
    for(auto&& h : mesh.connect)
    {
        if((cc++ % (cn / 10)) == 0) printf("...writing connect %f%%\n", 100.0f * cc / cn);
        fprintf(f, "HE %lld %lld %lld %lld %lld %lld %lld %lld\n", h[0], h[1], h[2], h[3], h[4], h[5], h[6], h[7]);
    }
    fclose(f);
}
void output(Mesh const& mesh)
{
    //output_csv(mesh);
    output_txt(mesh);
    output_tec(mesh);
}

void output_csv(Mesh const& mesh, Mode const& mode)
{
    char buf[33];
    snprintf(buf, 32, "mode-%lld.csv", mode.number);
    FILE* f;
    f = fopen(buf, "w");
    for(int64_t i = 0; i < mode.displacement.size(); ++i)
        fprintf(f, "%lf,%lf,%lf,%lf,%lf,%lf,%lf\n",
            mesh.nodes[i].x,
            mesh.nodes[i].y,
            mesh.nodes[i].z,
            sqrt(mode.displacement[i].x * mode.displacement[i].x + mode.displacement[i].y * mode.displacement[i].y + mode.displacement[i].z * mode.displacement[i].z),
            mode.displacement[i].x,
            mode.displacement[i].y,
            mode.displacement[i].z);
    fclose(f);
}
void output_txt(Mesh const& mesh, Mode const& mode)
{
    char buf[33];
    snprintf(buf, 32, "mode-%lld.txt", mode.number);
    FILE* f = fopen(buf, "w");
    for(int64_t i = 0; i < mode.displacement.size(); ++i) {
        auto&& p = mode.displacement[i];
        fprintf(f, "%lf %lf %lf\n", p.x, p.y, p.z);
    }
    fclose(f);
}
void output_tec(Mesh const& mesh, Mode const& mode)
{
    char buf[33];
    snprintf(buf, 32, "mode-%lld.tec", mode.number);
    FILE* f = fopen(buf, "w");
    fprintf(f, R"(TITLE = "Example: FE-Volume Brick Data"
VARIABLES = "X", "Y", "Z", "dx", "dy", "dz", "mag"
ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEBRICK
)", (int)mesh.nodes.size(), (int)mesh.connect.size());
    for(int64_t i = 0; i < mode.displacement.size(); ++i) {
        auto&& p = mode.displacement[i];
        auto&& n = mesh.nodes[i];
        fprintf(f, "%lf %lf %lf %lf %lf %lf %lf\n", n.x, n.y, n.z, p.x, p.y, p.z,
            sqrt(mode.displacement[i].x * mode.displacement[i].x + mode.displacement[i].y * mode.displacement[i].y + mode.displacement[i].z * mode.displacement[i].z));
    }
    fprintf(f, "\n");
    for(int i = 0; i < mesh.connect.size(); ++i) {
        fprintf(f, "%lld %lld %lld %lld %lld %lld %lld %lld\n",
            mesh.connect[i][0]+ 1,
            mesh.connect[i][1]+ 1,
            mesh.connect[i][2]+ 1,
            mesh.connect[i][3]+ 1,
            mesh.connect[i][4]+ 1,
            mesh.connect[i][5]+ 1,
            mesh.connect[i][6]+ 1,
            mesh.connect[i][7]+ 1);
    }
    fclose(f);
}
void output(Mesh const& mesh, Mode const& mode)
{
    //output_csv(mesh, mode);
    output_txt(mesh, mode);
    output_tec(mesh, mode);
}

Mode genMode(Mesh const& mesh, double L, double A, int64_t n)
{
    Mode rval;
    rval.number = n;
    int cn = mesh.nodes.size(), cc = 0;
    for(int64_t i = 0; i < mesh.nodes.size(); ++i)
    {
        if((cc++ % (cn / 10)) == 0) printf("...mode %lld processed %f%% nodes\n", n, 100.f * cc / cn);
        PT p, C = mesh.nodes[i];
        int combos[] = {
            1,
            2,
            4,
            3,
            5,
            7

        };
        int useX = combos[n % 6] & 1;
        int useY = combos[n % 6] & 2;
        int useZ = combos[n % 6] & 4;
        auto dist = [=](PT const& p1, PT const& p2) -> double {
            auto sqr = [](double d) -> double { return d*d; };
            return sqrt(sqr(p1.x - p2.x) + sqr(p1.y - p2.y) + sqr(p1.z - p2.z)) / L;
        };
        auto d = dist(C, PT{0.0, 0.0, 0.0});
        if(d < 0.1) d = 0.1;
        double sx = (A/n/(1/d)) * sin(3.14159 * dist(C, PT{0.0, 0.0, 0.0}) * n);
        p.x = sx * useX / (useX + useY + useZ);
        p.y = sx * useY / (useX + useY + useZ);
        p.z = sx * useZ / (useX + useY + useZ);
        rval.displacement.push_back(p);
    }
    return rval;
}

int main(int argc, char* argv[])
{
    if(argc < 5) {
        printf("Usage: %s N L A M\n", argv[0]);
        exit(1);
    }

    int64_t N = atoi(argv[1]);
    double L = atof(argv[2]);
    double A = atof(argv[3]);
    int64_t M = atoi(argv[4]);
    printf("Probably using %d\n", omp_get_max_threads());

    printf("Generating mesh\n");
    auto mesh = genMesh(N, L);
    printf("Writing mesh\n");
    output(mesh);

    int i;
    #pragma omp parallel for
    for(i = 0; i < M; ++i)
    {
        printf("Generating mode %ld\n", i + 1);
        auto mode = genMode(mesh, L, A, i + 1);
        printf("Writing mode %ld\n", i + 1);
        output(mesh, mode);
    }
}
