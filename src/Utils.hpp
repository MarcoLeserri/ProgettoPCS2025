#pragma once

#include <iostream>
#include "Polygonal.hpp"

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;

namespace PolygonalLibrary
{
bool ImportMesh(Polygonal& mesh);
bool ImportCell0Ds(Polygonal& mesh);
bool ImportCell1Ds(Polygonal& mesh);
bool ImportCell2Ds(Polygonal& mesh);
}