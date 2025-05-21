#pragma once

#include <sstream>
#include <iostream>
#include <fstream>
#include "Polygonal.hpp"

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;

namespace PolygonalLibrary
{
bool ImportMesh(Polygonal& mesh, int q);
bool ImportCell0Ds(Polygonal& mesh, int q);
bool ImportCell1Ds(Polygonal& mesh, int q);
bool ImportCell2Ds(Polygonal& mesh, int q);
}

struct Vector3dComparator {
    bool operator()(const Eigen::Vector3d& a, const Eigen::Vector3d& b) const {
        double eps = 1e-14;
        if (std::abs(a[0] - b[0]) > eps) return a[0] < b[0];
        if (std::abs(a[1] - b[1]) > eps) return a[1] < b[1];
        if (std::abs(a[2] - b[2]) > eps) return a[2] < b[2];
        return false; // considerati uguali
        std::cout << "Comparo:\nA = " << a.transpose() << "\nB = " << b.transpose() << std::endl;

    }
};

int VerificaEInserisci(Vector3d& Coord, map<Vector3d, int, Vector3dComparator>& mappa, Polygonal& mesh);
bool TriangTotC_1(int b, int c, Polygonal& mesh, Polygonal& meshTriang);
