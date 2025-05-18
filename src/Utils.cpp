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

bool ImportMesh(Polygonal& mesh, int q)
{

    if(!ImportCell0Ds(mesh, q))
        return false;

    if(!ImportCell1Ds(mesh, q))
        return false;

    if(!ImportCell2Ds(mesh, q))
        return false;

    return true;

}


bool ImportCell0Ds(Polygonal& mesh, int q)
{	
	double err = 1.0e-16;
	string Filename;
	if (abs(q-3) < err) {
		Filename = "./Cell0TDs.csv";
	}
	else if (abs(q-4) < err) {
		Filename = "./Cell0ODs.csv";
	}
	else if (abs(q-5) < err) {
		Filename = "./Cell0IDs.csv";
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
        istringstream converter(line); // istringstream legge id,  coord x e y e z 
		stringstream s(line);
		string value;
		getline(s, value, ';');
		unsigned int id = stoi(value);
		mesh.Cell0DsID.push_back(id);
		getline(s, value, ';');
		mesh.Cell0DsCoordinates(0, id) = stod(value);
		getline(s, value, ';');
		mesh.Cell0DsCoordinates(1, id) = stod(value);
		getline(s, value, ';');
		mesh.Cell0DsCoordinates(2, id) = stod(value);

    }
    return true;
}
// ***************************************************************************

bool ImportCell1Ds(Polygonal& mesh, int q)
{
    double err = 1.0e-16;
	string Filename;
	if (abs(q-3) < err) {
		Filename = "./Cell1TDs.csv";
	}
	else if (abs(q-4) < err) {
		Filename = "./Cell1ODs.csv";
	}
	else if (abs(q-5) < err) {
		Filename = "./Cell1IDs.csv";
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
bool ImportCell2Ds(Polygonal& mesh, int q)
{
	
	double err = 1.0e-16;
	string Filename;
	if (abs(q-3) < err) {
		Filename = "./Cell2TDs.csv";
	}
	else if (abs(q-4) < err) {
		Filename = "./Cell2ODs.csv";
	}
	else if (abs(q-5) < err) {
		Filename = "./Cell2IDs.csv";
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

    mesh.Cell2DsID.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsVertices.reserve(mesh.NumCell2Ds);
    mesh.Cell2DsEdges.reserve(mesh.NumCell2Ds);

    for (const string& line : listLines)
    {
        istringstream converter(line);

        unsigned int id;
        array<unsigned int, 3> vertices;
        array<unsigned int, 3> edges;

        string val;
        getline(converter, val, ';'); 
        id = stoi(val);
        
		getline(converter, val, ';');
        for( int i = 0; i < 3; i++){
        	getline(converter, val, ';');
            vertices[i] = stoi(val);
        }
		getline(converter, val, ';');
        for( int i = 0; i < 3; i++){
            getline(converter, val, ';');
            edges[i] = stoi(val);
        }

        mesh.Cell2DsID.push_back(id);
        mesh.Cell2DsVertices.push_back(vertices);
        mesh.Cell2DsEdges.push_back(edges);
    }

    return true;
}

}

bool TriangFaceC_1(Polygonal& meshTriang, int IdFace, map<string, array<double, 3>> VertFace, int n, int Face){
	int numPunti = 0;
	for( int i = 0; i < n + 1; i++){
		numPunti += (i + 1);
	}
	
	int numEdges = 0;
	for( int i = 1; i < n + 1; i++){
		numEdges += (3 * i);
	}
	
	int numFaces = n * n;
	
	meshTriang.NumCell0Ds = numPunti * numFaces;
	meshTriang.NumCell1Ds = numEdges * numFaces;
	meshTriang.NumCell2Ds = numFaces;
	
	meshTriang.Cell0DsID.reserve(meshTriang.NumCell0Ds);
    meshTriang.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, meshTriang.NumCell0Ds);
    meshTriang.Cell1DsID.reserve(meshTriang.NumCell1Ds);
    meshTriang.Cell1DsExtrema = Eigen::MatrixXi(2, meshTriang.NumCell1Ds);
    meshTriang.Cell2DsID.reserve(meshTriang.NumCell2Ds);
    meshTriang.Cell2DsVertices.reserve(meshTriang.NumCell2Ds);
    meshTriang.Cell2DsEdges.reserve(meshTriang.NumCell2Ds);

	
	for( int i = 0; i < numPunti; i++){
		meshTriang.Cell0DsID.push_back(i + numPunti * Face);
	}
	
	for( int i = 0; i < numEdges; i++){
		meshTriang.Cell1DsID.push_back(i);
	}
	
	for( int i = 0; i < numFaces; i++){
		meshTriang.Cell2DsID.push_back(i);
	}
	
	//triangolazione lati
	int Id = 0;
	
	meshTriang.Cell0DsCoordinates(0,Id + numPunti * Face) = VertFace["1"][0];
	meshTriang.Cell0DsCoordinates(1,Id + numPunti * Face) = VertFace["1"][1];
	meshTriang.Cell0DsCoordinates(2,Id + numPunti * Face) = VertFace["1"][2];

	double AddX01 = abs(VertFace["0"][0] - VertFace["1"][0]) / n;
	double AddY01 = abs(VertFace["0"][1] - VertFace["1"][1]) / n;
	double AddZ01 = abs(VertFace["0"][2] - VertFace["1"][2]) / n;

	double AddX12 = abs(VertFace["2"][0] - VertFace["1"][0]) / n;
	double AddY12 = abs(VertFace["2"][1] - VertFace["1"][1]) / n;
	double AddZ12 = abs(VertFace["2"][2] - VertFace["1"][2]) / n;
	
	double X0 = VertFace["1"][0];
	double Y0 = VertFace["1"][1];
	double Z0 = VertFace["1"][2];
	Id += 1;
	
	cout << X0<<Y0<<Z0<<endl;
	
	for( int i = 1; i < n + 1; i++){
		
		double XEstrS;
		double YEstrS;
		double ZEstrS;
		double XEstrD;
		double YEstrD;
		double ZEstrD;
		
		if(VertFace["1"][0] < VertFace["0"][0]){
			XEstrS = X0 + ( i * AddX01);
		}
		else{
			XEstrS = X0 - ( i * AddX01);
		}
		
		if(VertFace["1"][1] < VertFace["0"][1]){
			YEstrS = Y0 + ( i * AddY01);
		}
		else{
			YEstrS = Y0 - ( i * AddY01);
		}
		
		if(VertFace["1"][2] < VertFace["0"][2]){
			ZEstrS = Z0 + ( i * AddZ01);
		}
		else{
			ZEstrS = Z0 - ( i * AddZ01);
		}
		
		if(VertFace["1"][0] < VertFace["2"][0]){
			XEstrD = X0 + ( i * AddX12);
		}
		else{
			XEstrD = X0 - ( i * AddX12);
		}
		
		if(VertFace["1"][1] < VertFace["2"][1]){
			YEstrD = Y0 + ( i * AddY12);
		}
		else{
			YEstrD = Y0 - ( i * AddY12);
		}
		
		if(VertFace["1"][2] < VertFace["2"][2]){
			ZEstrD = Z0 + ( i * AddZ12);
		}
		else{
			ZEstrD = Z0 - ( i * AddZ12);
		}
		
		meshTriang.Cell0DsCoordinates(0,Id + numPunti * Face) = XEstrS;
		meshTriang.Cell0DsCoordinates(1,Id + numPunti * Face) = YEstrS;
		meshTriang.Cell0DsCoordinates(2,Id + numPunti * Face) = ZEstrS;
		Id += 1;
		
		double DistEstX = abs( XEstrS - XEstrD) / i;
		double DistEstY = abs( YEstrS - YEstrD) / i;
		double DistEstZ = abs( ZEstrS - ZEstrD) / i;
			
		for( int j = 1; j < i; j++){
			
			if(XEstrS < XEstrD){
				meshTriang.Cell0DsCoordinates(0,Id + numPunti * Face) = XEstrS + ( j * DistEstX );
			}
			else{
				meshTriang.Cell0DsCoordinates(0,Id + numPunti * Face) = XEstrS - ( j * DistEstX );
			}
			
			if(YEstrS < YEstrD){
				meshTriang.Cell0DsCoordinates(1,Id + numPunti * Face) = YEstrS + ( j * DistEstY );
			}
			else{
				meshTriang.Cell0DsCoordinates(1,Id + numPunti * Face) = YEstrS - ( j* DistEstY );
			}
			
			if(ZEstrS < ZEstrD){
				meshTriang.Cell0DsCoordinates(2,Id + numPunti * Face) = ZEstrS + ( j * DistEstZ );
			}
			else{
				meshTriang.Cell0DsCoordinates(2,Id + numPunti * Face) = ZEstrS - ( j * DistEstZ );
			}
			Id += 1;
		}
		
		meshTriang.Cell0DsCoordinates(0,Id + numPunti * Face) = XEstrD;
		meshTriang.Cell0DsCoordinates(1,Id + numPunti * Face) = YEstrD;
		meshTriang.Cell0DsCoordinates(2,Id + numPunti * Face) = ZEstrD;
		Id += 1;
	
	}
	return true;

}

bool TriangTotC_1(int b, int c, Polygonal& mesh, Polygonal& meshTriang){
	
	int n;
	double err=1.0e-16;
	if(abs(b) < err)
		n = c;
	else if(abs(c) < err)
		n = b;
		
	int numFaces = mesh.Cell2DsVertices.size(); //numero di facce nel poligono
	array<double, 3> CoordFace;
	
	cout <<numFaces<<endl;
	
	for( int i = 0; i < numFaces; i++){
		int IdFace = i;
		map<string, array<double, 3>> VertFace;
		for( int j = 0; j < 3; j++){
			int IdVert = mesh.Cell2DsVertices[i][j];
			double X = mesh.Cell0DsCoordinates(0,IdVert);
			double Y = mesh.Cell0DsCoordinates(1,IdVert);
			double Z = mesh.Cell0DsCoordinates(2,IdVert);
			
			CoordFace[0] = X;
			CoordFace[1] = Y;
			CoordFace[2] = Z;
			VertFace[to_string(j)] = CoordFace;
		}
		
		int Face = i;
		TriangFaceC_1(meshTriang, IdFace, VertFace, n, Face);
	}
	//stampate i txt
	
	
	return true;
}

