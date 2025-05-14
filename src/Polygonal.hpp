#pragma once

#include <iostream>
#include "Eigen/Eigen"s

using namespace std;
using namespace Eigen;

namespace PolygonalLibrary {

struct Polygonal
{	
	unsigned int NumCell0Ds = 0; // num punti
	unsigned int NumCell1Ds = 0; // num segmenti
	unsigned int NumCell2Ds = 0; // num poligoni 
	
	vector<unsigned int> Cell0DsID = {}; //registra ID punti
	vector<unsigned int> Cell1DsID = {};
	vector<unsigned int> Cell2DsID = {};
	
	Eigen::MatrixXd Cell0DsCoordinates = {}; 
	Eigen::MatrixXi Cell1DsExtrema = {};
	vector<array<unsigned int,3>> Cell2DsVertices = {};
	vector<array<unsigned int,3>> Cell2DsEdges = {};
	
};