#include <sstream>
#include <iostream>
#include <fstream>
#include "Utils.hpp"
#include "Eigen/Eigen"
#include <map>

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;

namespace PolygonalLibrary{

bool ImportMesh(PolygonalMesh& mesh)
{

    if(!ImportCell0Ds(mesh))
        return false;

    if(!ImportCell1Ds(mesh))
        return false;

    if(!ImportCell2Ds(mesh))
        return false;

    return true;

}


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
			VertFace[to_string(i)] = CoordFace;
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
	
	//triangolazione lati
	Cell0DsCoordinates[0][0] = VertFace["1"][0];
	Cell0DsCoordinates[1][0] = VertFace["1"][1];
	Cell0DsCoordinates[2][0] = VertFace["1"][2];

	double AddX01 = abs(VertFace["0"][0] - VertFace["1"][0]) / n;
	double AddY01 = abs(VertFace["0"][1] - VertFace["1"][1]) / n;
	double AddZ01 = abs(VertFace["0"][2] - VertFace["1"][2]) / n;

	double AddX12 = abs(VertFace["2"][0] - VertFace["1"][0]) / n;
	double AddY12 = abs(VertFace["2"][1] - VertFace["1"][1]) / n;
	double AddZ12 = abs(VertFace["2"][2] - VertFace["1"][2]) / n;
	
	double X0 = VertFace["1"][0];
	double Y0 = VertFace["1"][1];
	double Z0 = VertFace["1"][2];
	int Id = 1;
	
	for( unsigned int i = 1; i < n + 1; i++){
		
		if(VertFace["1"][0] < VertFace["0"][0]){
			double XEstrS = X0 + ( i * AddX01);
		}
		else{
			double XEstrS = X0 + ( (n - i) * AddX01);
		}
		
		if(VertFace["1"][1] < VertFace["0"][1]){
			double YEstrS = Y0 + ( i * AddY01);
		}
		else{
			double YEstrS = Y0 + ( (n - i) * AddY01);
		}
		
		if(VertFace["1"][2] < VertFace["0"][2]){
			double ZEstrS = Z0 + ( i * AddZ01);
		}
		else{
			double ZEstrS = Z0 + ( (n - i) * AddZ01);
		}
		
		if(VertFace["1"][0] < VertFace["2"][0]){
			double XEstrD = X0 + ( i * AddX12);
		}
		else{
			double XEstrD = X0 + ( (n - i) * AddX12);
		}
		
		if(VertFace["1"][1] < VertFace["2"][1]){
			double YEstrD = Y0 + ( i * AddY12);
		}
		else{
			double YEstrD = Y0 + ( (n - i) * AddY12);
		}
		
		if(VertFace["1"][2] < VertFace["2"][2]){
			double ZEstrD = Z0 + ( i * AddZ12);
		}
		else{
			double ZEstrD = Z0 + ( (n - i) * AddZ12);
		}
		
		Cell0DsCoordinates[0][Id] = XEstrS;
		Cell0DsCoordinates[1][Id] = YEstrS;
		Cell0DsCoordinates[2][Id] = ZEstrS;
		Id += 1;
		
		double DistEstX = abs( XEstrS - XEstrD) / i;
		double DistEstY = abs( YEstrS - YEstrD) / i;
		double DistEstZ = abs( ZEstrS - ZEstrD) / i;
			
		for( unsigned int j = 1; j < i; j++){
			
			if(XEstrS < XEstrD){
				Cell0DsCoordinates[0][Id] = XEstrS + ( j * DistEstX );
			}
			else{
				Cell0DsCoordinates[0][Id] = XEstrS + ( (i - j) * DistEstX );
			}
			
			if(YEstrS < YEstrD){
				Cell0DsCoordinates[1][Id] = YEstrS + ( j * DistEstY );
			}
			else{
				Cell0DsCoordinates[1][Id] = YEstrS + ( (i - j) * DistEstY );
			}
			
			if(ZEstrS < ZEstrD){
				Cell0DsCoordinates[2][Id] = ZEstrS + ( j * DistEstZ );
			}
			else{
				Cell0DsCoordinates[2][Id] = ZEstrS + ( (i - j) * DistEstZ );
			}
			Id += 1;
		}
		
		Cell0DsCoordinates[0][Id] = XEstrD;
		Cell0DsCoordinates[1][Id] = YEstrD;
		Cell0DsCoordinates[2][Id] = ZEstrD;
		Id += 1;
	
	}
	return true;

}