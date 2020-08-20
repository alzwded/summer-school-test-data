// compile with:
//    cl.exe /EHa /O2 /openmp /std:c++17 genplate.cpp
//    gcc -fopenmp -o genplate -O3 genplate.cpp -lm
#include <cstdio>
#include <map>
#include <vector>
#include <cstdlib>
#include <cmath>
#include <array>

struct PT
{
    double x, y, z;
};

typedef std::array<int, 4> QD;

struct Mesh
{
    std::vector<PT> nodes;
    std::vector<int> nodeIds;
    std::vector<QD> connect;
};

struct Mode
{
    int number;
    std::vector<PT> displacement;
};

Mesh genMesh(int N, double L)
{
    FILE*g = fopen("debug2.txt", "w");
    Mesh rval;
    rval.nodeIds.resize(N * N, -1);
    for(int i = 0; i < N; ++i) {
        for(int j = 0; j < N; ++j) {
            // poke a hole in the mesh if it's big enough:-)
            if((N >= 16
                &&! (i >= 4.0 / 8 * N && i < 5.0 / 8 * N &&
                     j >= 4.0 / 8 * N && j < 7.0 / 8 * N))
            || N < 16)
            {
                PT p = { -L/2.0 + L / (N - 1) * i,
                         -L/2.0 + L / (N - 1) * j,
                         0.0 };
                rval.nodeIds[i * N + j] = rval.nodes.size();
                rval.nodes.push_back(p);
            }
            if(i < N - 1 && j < N - 1) {
                rval.connect.push_back({
                    i * N + j,
                    i * N + j + 1,
                    (i + 1) * N + j,
                    (i + 1) * N + j + 1
                });
                fprintf(g, "%d %d %d %d\n",
                    rval.connect.back()[0],
                    rval.connect.back()[1],
                    rval.connect.back()[2],
                    rval.connect.back()[3]);
            }
        }
    }
    FILE* f = fopen("debug.txt", "w");
    for(int i = 0; i < N*N; ++i) {
        fprintf(f, "%d*%d=%d %d\n", i / N, i % N, i, rval.nodeIds[i]);
    }
    // fixup nodeIds
    decltype(rval.connect) newConnect;
    for(int64_t i = 0; i < rval.connect.size(); ++i) {
        bool erase = false;
        for(int64_t j = 0; j < rval.connect[i].size(); ++j) {
            fprintf(g, "%d => %d\n", rval.connect[i][j], rval.nodeIds[rval.connect[i][j]]);
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
    for(auto&& n : mesh.nodes)
        fprintf(f, "%lf %lf %lf\n", n.x, n.y, n.z);
    fclose(f);
    f = fopen("connect.txt", "w");
    for(auto&& q : mesh.connect)
        fprintf(f, "QU %d %d %d %d\n", q[0], q[1], q[3], q[2]);
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
    snprintf(buf, 32, "mode-%d.csv", mode.number);
    FILE* f;
    f = fopen(buf, "w");
    for(int i = 0; i < mode.displacement.size(); ++i)
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
    snprintf(buf, 32, "mode-%d.txt", mode.number);
    FILE* f = fopen(buf, "w");
    for(int i = 0; i < mode.displacement.size(); ++i) {
        auto&& p = mode.displacement[i];
        fprintf(f, "%lf %lf %lf\n", p.x, p.y, p.z);
    }
    fclose(f);
}
void output_tec(Mesh const& mesh, Mode const& mode)
{
    char buf[33];
    snprintf(buf, 32, "mode-%d.tec", mode.number);
    FILE* f = fopen(buf, "w");
    fprintf(f, R"(TITLE = "Example: FE-Volume Brick Data"
VARIABLES = "X", "Y", "Z", "dX", "dY", "dZ", "mag"
ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL
)", (int)mesh.nodes.size(), (int)mesh.connect.size());
    for(int64_t i = 0; i < mode.displacement.size(); ++i) {
        auto&& p = mode.displacement[i];
        auto&& n = mesh.nodes[i];
        fprintf(f, "%lf %lf %lf %lf %lf %lf %lf\n", n.x, n.y, n.z, p.x, p.y, p.z,
            sqrt(mode.displacement[i].x * mode.displacement[i].x + mode.displacement[i].y * mode.displacement[i].y + mode.displacement[i].z * mode.displacement[i].z));
    }
    fprintf(f, "\n");
    for(int i = 0; i < mesh.connect.size(); ++i) {
        fprintf(f, "%d %d %d %d\n",
            mesh.connect[i][0]+ 1,
            mesh.connect[i][1]+ 1,
            mesh.connect[i][3]+ 1,
            mesh.connect[i][2]+ 1);
    }
    fclose(f);
}
void output(Mesh const& mesh, Mode const& mode)
{
    //output_csv(mesh, mode);
    output_txt(mesh, mode);
    output_tec(mesh, mode);
}

Mode genMode(Mesh const& mesh, double L, double A, int n)
{
    Mode rval;
    rval.number = n;
    for(int i = 0; i < mesh.nodes.size(); ++i)
    {
        PT p, C = mesh.nodes[i];
        p.x = p.y = 0.0;
        double d = sqrt(C.x * C.x + C.y * C.y + C.z * C.z);
        p.z = (A / n) * sin(3.14159 * d / L * n);
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

    int N = atoi(argv[1]);
    double L = atof(argv[2]);
    double A = atof(argv[3]);
    int M = atoi(argv[4]);

    printf("Generating mesh\n");
    auto mesh = genMesh(N, L);
    printf("Writing mesh\n");
    output(mesh);

    int i;
    #pragma omp parallel for
    for(i = 0; i < M; ++i)
    {
        printf("Generating mode %d\n", i + 1);
        auto mode = genMode(mesh, L, A, i + 1);
        printf("Writing mode %d\n", i + 1);
        output(mesh, mode);
    }
}
