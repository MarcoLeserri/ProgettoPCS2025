#pragma once

#include <iostream>
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

bool TriangTotC_1(int b, int c, Polygonal& mesh, Polygonal& meshTriang);
bool TriangFaceC_1(Polygonal& meshTriang, int IdFace, map<char, array<double, 3>> VertFace, int n);