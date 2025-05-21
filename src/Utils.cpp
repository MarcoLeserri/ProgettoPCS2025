#include <sstream>
#include <iostream>
#include <fstream>
#include "Utils.hpp"
#include "Eigen/Eigen"
#include <map>
#include <vector>
#include <iomanip>

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
		cout << value << endl;
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

int VerificaEInserisci(Vector3d& Coord, map<Vector3d, int, Vector3dComparator>& mappa, Polygonal& mesh) {
    Coord = Coord / Coord.norm(); // normalizzazione
    auto it = mappa.find(Coord);
    if (it != mappa.end()) {
	    std::cout << "Già presente: " << Coord.transpose() << " → id " << it->second << std::endl;
        return it->second;
    } else {
        int id = mappa.size();
        mappa[Coord] = id;
        std::cout << fixed << setprecision(16) << "Inserisco punto " << id << ": " << Coord.transpose() << std::endl;
        mesh.Cell0DsCoordinates(0, id) = Coord[0];
        mesh.Cell0DsCoordinates(1, id) = Coord[1];
        mesh.Cell0DsCoordinates(2, id) = Coord[2];
        mesh.Cell0DsID.push_back(id);
        return id;
    }
}


bool TriangTotC_1(int b, int c, Polygonal& mesh, Polygonal& meshTriang){
	
	int n;
	double err=1.0e-16;
	if(abs(b) < err)
		n = c;
	else if(abs(c) < err)
		n = b;
	else
		cerr << "Error: both b and c are not zero" << endl;
		
	cout << "n è: " << n <<endl;
		
	int numFaces = mesh.Cell2DsVertices.size(); //numero di facce nel poligono
	cout << "il numero di facce del poligono è: " << numFaces << endl;
	int NumPuntiF = ( ( n + 1 ) * ( n + 2 ) ) / 2; //numero totale di punti in una faccia
	cout << "il numero di punti in una faccia è: " << NumPuntiF << endl;
	meshTriang.NumCell0Ds = NumPuntiF * numFaces; //da cambiare
	meshTriang.NumCell2Ds = n * n * numFaces;
	cout << "il numero di facce triangolate è: " << meshTriang.NumCell2Ds << endl;
	
	meshTriang.Cell0DsID.reserve(meshTriang.NumCell0Ds);
    meshTriang.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, meshTriang.NumCell0Ds);
    meshTriang.Cell1DsID.reserve(meshTriang.NumCell1Ds);
    meshTriang.Cell1DsExtrema = Eigen::MatrixXi(2, meshTriang.NumCell1Ds);
    meshTriang.Cell2DsID.reserve(meshTriang.NumCell2Ds);
    meshTriang.Cell2DsVertices.reserve(meshTriang.NumCell2Ds);
    meshTriang.Cell2DsEdges.reserve(meshTriang.NumCell2Ds);
	
	
	map<Vector3d, int, Vector3dComparator> ControlloPunti;
	Vector3d CoordA; //coordinate A
	Vector3d CoordB; //coordinate B
	Vector3d CoordC; //coordinate C
	Vector3d CoordP; //coordinate nuovo punto
	array<unsigned int, 3> NewFace;
	//map< Vector3d, int> ControlloPunti;
		
	for(int k = 0; k < numFaces; k++){
		
		int numfaccia = k;
		
		int id;
		id = mesh.Cell2DsVertices[numfaccia][0];
		cout << id << endl;
		CoordA[0] = mesh.Cell0DsCoordinates(0,id);
		CoordA[1] = mesh.Cell0DsCoordinates(1,id);
		CoordA[2] = mesh.Cell0DsCoordinates(2,id);
		cout << CoordA[0] << " " << CoordA[1] << " " << CoordA[2] << endl;
		
		id = mesh.Cell2DsVertices[numfaccia][1];
		cout << id << endl;
		CoordB[0] = mesh.Cell0DsCoordinates(0,id);
		CoordB[1] = mesh.Cell0DsCoordinates(1,id);
		CoordB[2] = mesh.Cell0DsCoordinates(2,id);
		cout << CoordB[0] << " " << CoordB[1] << " " << CoordB[2] << endl;
		
		id = mesh.Cell2DsVertices[numfaccia][2];
		cout << id << endl;
		CoordC[0] = mesh.Cell0DsCoordinates(0,id);
		CoordC[1] = mesh.Cell0DsCoordinates(1,id);
		CoordC[2] = mesh.Cell0DsCoordinates(2,id);
		cout << CoordC[0] << " " << CoordC[1] << " " << CoordC[2] << endl;
		
		for(int i = 0; i < n; i++){ //conteggio piani
			cout << "Triangoli dritti" << endl;
			for(int j = 0; j < (n - i); j++){
				CoordP = ((double(i)/n) * CoordA) + ((double(n - (i + j))/n) * CoordB) + ((double(j)/n) * CoordC);
				cout << fixed << setprecision(16) << CoordP[0] << " " << CoordP[1] << " " << CoordP[2] << endl;
				id = VerificaEInserisci(CoordP, ControlloPunti, meshTriang);
				NewFace[0] = id;
				
				CoordP = ((double(i)/n) * CoordA) + ((double(n - (i + j + 1))/n) * CoordB) + ((double(j + 1)/n) * CoordC);
				cout << fixed << setprecision(16) << CoordP[0] << " " << CoordP[1] << " " << CoordP[2] << endl;
				id = VerificaEInserisci(CoordP, ControlloPunti, meshTriang);
				NewFace[1] = id;
				
				CoordP = ((double(i + 1)/n) * CoordA) + ((double(n - (i + j + 1))/n) * CoordB) + ((double(j)/n) * CoordC);
				cout << fixed << setprecision(16) << CoordP[0] << " " << CoordP[1] << " " << CoordP[2] << endl;
				id = VerificaEInserisci(CoordP, ControlloPunti, meshTriang);
				NewFace[2] = id;

				meshTriang.Cell2DsVertices.push_back(NewFace);
			}
			
			cout << "Triangoli inversi" << endl;
			for(int j = 0; j < (n - (i + 1)); j++){
				CoordP = ((double(i)/n) * CoordA) + ((double(n - (i + j + 1))/n) * CoordB) + ((double(j + 1)/n) * CoordC);
				cout << fixed << setprecision(16) << CoordP[0] << " " << CoordP[1] << " " << CoordP[2] << endl;
				id = VerificaEInserisci(CoordP, ControlloPunti, meshTriang);
				NewFace[0] = id;
				
				CoordP = ((double(i + 1)/n) * CoordA) + ((double(n - (i + j + 2))/n) * CoordB) + ((double(j + 1)/n) * CoordC);
				cout << fixed << setprecision(16) << CoordP[0] << " " << CoordP[1] << " " << CoordP[2] << endl;
				id = VerificaEInserisci(CoordP, ControlloPunti, meshTriang);
				NewFace[1] = id;
				
				CoordP = ((double(i + 1)/n) * CoordA) + ((double(n - (i + j + 1))/n) * CoordB) + ((double(j)/n) * CoordC);
				cout << fixed << setprecision(16) << CoordP[0] << " " << CoordP[1] << " " << CoordP[2] << endl;
				id = VerificaEInserisci(CoordP, ControlloPunti, meshTriang);
				NewFace[2] = id;
				
				meshTriang.Cell2DsVertices.push_back(NewFace);
			}
			
		}
	}
		
	return true;
}
