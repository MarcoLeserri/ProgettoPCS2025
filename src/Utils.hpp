#pragma once

#include <sstream>
#include <iostream>
#include <fstream>
#include "Polygonal.hpp"

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;

using coppia = pair<int, double>;
using grafo =  vector<vector<coppia>>;

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
    }
};

struct Vector2iComparator {
    bool operator()(const Eigen::Vector2i& a, const Eigen::Vector2i& b) const {
        if (a[0] != b[0]) return a[0] < b[0];
        return a[1] < b[1];
    }
};


int VerificaEInserisci(Vector3d& Coord, map<Vector3d, int, Vector3dComparator>& mappa, Polygonal& mesh);

int VerificaEInserisciDual(Vector3d& Coord, map<Vector3d, int, Vector3dComparator>& mappa, PolygonalDual& mesh);

array<unsigned int, 3> VerificaEInserisci2(array<unsigned int, 3> NewFace, map<Vector2i, int, Vector2iComparator>& mappa, Polygonal& mesh, array<unsigned int, 3>& Face);

vector<unsigned int> VerificaEInserisci2Dual(vector<unsigned int> NewFace, map<Vector2i, int, Vector2iComparator>& mappa, Eigen::MatrixXi& Cell1DsExtrema, PolygonalDual& mesh, vector<unsigned int>& Face);

bool TriangTotC_1(int b, int c, Polygonal& mesh, Polygonal& meshTriang);

vector<array<unsigned int, 3>> SearchFaces( int id, vector<array<unsigned int, 3>> VerticesF, int length);

bool DualTot(Polygonal& meshTriang, PolygonalDual& meshDual);

vector<array<unsigned int, 3>> SortFaces(int id, int IdOrd, vector<array<unsigned int,3>> FacesAd);

vector<Vector2i> SearchEdges(int id1, vector<Vector2i> Lati);

vector<coppia> SearchPoint(int id1, vector<Vector2i> EdgesList, map<Vector2i, double, Vector2iComparator> EdgesLength);

grafo BuildGraph(PolygonalDual& mesh, vector<Vector2i> EdgesList, map<Vector2i, double, Vector2iComparator> EdgesLength);

vector<unsigned int> ShortestPath(int id1, int id2, PolygonalDual& meshDual);

vector<double> ParaviewPoints(vector<unsigned int> percorso, PolygonalDual& mesh);

vector<double> ParaviewEdges(vector<unsigned int> percorso, PolygonalDual& mesh);

void FileTxt(const Polygonal& mesh);