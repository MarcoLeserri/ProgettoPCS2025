#include <sstream>
#include <iostream>
#include <fstream>
#include "Utils.hpp"
#include<map>

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;

namespace PolygonalLibrary{

bool ImportCell0Ds(Polygonal& mesh, int q)
{	
	double err = 1.0e-16;
	string Filename;
	if (abs(q-3) < err) {
		Filename = "./Cell0TDs.csv"
	}
	else if (abs(q-3) < err) {
		Filename = "./Cell0ODs.csv"
	}
	else if (abs(q-5) < err) {
		Filename = "./Cell0IDs.csv"
	}	
	
	
    ifstream file(Filename);

	if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
		listLines.push_back(line);
	
    file.close();
    listLines.pop_front();
    mesh.NumCell0Ds = listLines.size();

    if (mesh.NumCell0Ds == 0) 
    {
        cerr << "There is no cell 0D" << endl;
        return false;
    }
	
	mesh.Cell0DsID.reserve(mesh.NumCell0Ds);
    mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line); // istringstream legge id,  coord x e y 
		stringstream s(line);
		string value;
		getline(s, value, ';');
		unsigned int id = stoi(value);
		mesh.Cell0DsID.push_back(id);
		getline(s, value, ';');
		mesh.Cell0DsCoordinates(0, id) = stod(value);
		getline(s, value, ';');
		mesh.Cell0DsCoordinates(1, id) = stod(value);

    }
    return true;
}
// ***************************************************************************

bool ImportCell1Ds(Polygonal& mesh)
{
    double err = 1.0e-16;
	string Filename;
	if (abs(q-3) < err) {
		Filename = "./Cell1TDs.csv"
	}
	else if (abs(q-3) < err) {
		Filename = "./Cell1ODs.csv"
	}
	else if (abs(q-5) < err) {
		Filename = "./Cell1IDs.csv"
	}	

    ifstream file(Filename);

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();
    listLines.pop_front();

    mesh.NumCell1Ds = listLines.size();

    if (mesh.NumCell1Ds == 0)
    {
        cerr << "There is no cell 1D" << endl;
        return false;
    }

    mesh.Cell1DsID.reserve(mesh.NumCell1Ds);
    mesh.Cell1DsExtrema = Eigen::MatrixXi(2, mesh.NumCell1Ds);

    for (const string& line : listLines)
    {	istringstream converter(line); // istringstream legge id,  coord x e y 
		stringstream s(line);
		string value;
		getline(s, value, ';');
		unsigned int id = stoi(value);
		mesh.Cell1DsID.push_back(id);
		getline(s, value, ';');
		mesh.Cell1DsExtrema(0, id) = stoi(value);
		getline(s, value, ';');
		mesh.Cell1DsExtrema(1, id) = stoi(value);

    }

    return true;
}

// ***************************************************************************
bool ImportCell2Ds(Polygonal& mesh)
{
	
	double err = 1.0e-16;
	string Filename;
	if (abs(q-3) < err) {
		Filename = "./Cell2TDs.csv"
	}
	else if (abs(q-3) < err) {
		Filename = "./Cell2ODs.csv"
	}
	else if (abs(q-5) < err) {
		Filename = "./Cell2IDs.csv"
	}	

    ifstream file;
    file.open(Filename);

    if(file.fail())
        return false;

    list<string> listLines;
    string line;
    while (getline(file, line))
        listLines.push_back(line);

    file.close();
    listLines.pop_front();

    mesh.NumCell2Ds = listLines.size();

    if (mesh.NumCell2Ds == 0)
    {
        cerr << "There is no cell 2D" << endl;
        return false;
    }

    mesh.Cell2DsId.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsVertices.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsEdges.reserve(mesh.NumCell2Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        array<unsigned int, 3> vertices;
        array<unsigned int, 3> edges;
		unsigned int m;

        converter >>  id;
		converter >> m;
        for(unsigned int i = 0; i < 3; i++)
            converter >> vertices[i];
		converter >> m;
        for(unsigned int i = 0; i < 3; i++)
            converter >> edges[i];

        mesh.Cell2DsId.push_back(id);
        mesh.Cell2DsVertices.push_back(vertices);
        mesh.Cell2DsEdges.push_back(edges);
    }

    return true;
}

}

bool TriangTotC_1(int b, int c, Polygonal& mesh, Polygonal& meshTriang){
	
	int n;
	double err=1.0e-16;
	if(abs(b) < err)
		n = c;
	else if(abs(c) < err)
		n = b;
		
	int numFaces = length(mesh.Cell2DsVertices); //numero di facce nel poligono
	array<double, 3> CoordFace;
	
	for( unsigned int i = 0; i < numFaces; i++){
		int IdFace = i;
		map<int, array<double, 3>> VertFace;
		for( unsigned int j = 0; j < 3; j++){
			int IdVert = mesh.Cell2DsVertices[i][j];
			double X = mesh.Cell0DsCoordinates[IdVert][0];
			double Y = mesh.Cell0DsCoordinates[IdVert][1];
			double Z = mesh.Cell0DsCoordinates[IdVert][2];
			CoordFace[0] = X;
			CoordFace[1] = Y;
			CoordFace[2] = Z;
			VertFace[to_string(i + 1)] = CoordFace;
		}
		TriangFaceC_1(meshTriang, IdFace, CoordFace, n);
	}
	//stampate i txt
	
	
	return true;
}

bool TrianfFaceC_1(Polygonal& meshTriang, int IdFace, map<int, array<double, 3>> VertFace, int n){
	int numPunti;
	for( unsigned int i = 0; i < b + 1; i++){
		numPunti += (i + 1);
	}
	
	int numEdges;
	for( unsigned int i = 1; i < b + 1; i++){
		numEdges += (3 * i);
	}
	
	int numFaces = n * n;
	
	meshTriang.numCell0Ds = numPunti;
	meshTriang.numCell1Ds = numEdges;
	meshTriang.numCell2Ds = numFaces;
	
	for( unsigned int i = 0; i < numPunti; i++){
		meshTriang.Cell0DsID.push_back(i);
	}
	
	for( unsigned int i = 0; i < numEdges; i++){
		meshTriang.Cell1DsID.push_back(i);
	}
	
	for( unsigned int i = 0; i < numFaces; i++){
		meshTriang.Cell2DsID.push_back(i);
	}
	
	array<double, 3> CoordA = VertFace["1"]
	array<double, 3> CoordB = VertFace["2"]
	array<double, 3> CoordC = VertFace["3"]
	
	//triangolazione lati
	for( unsigned int i = 0; i < n ;i++){
		double AddX = abs(CoordA[0] - CoordB[0]) / n;
		double AddY = abs(CoordA[1] - CoordB[1]) / n;
		double AddZ = abs(CoordA[2] - CoordB[2]) / n;
		
		
	}
}