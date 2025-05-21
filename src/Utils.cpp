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
        return it->second;
    } else {
        int id = mappa.size();
        mappa[Coord] = id;
        mesh.Cell0DsCoordinates(0, id) = Coord[0];
        mesh.Cell0DsCoordinates(1, id) = Coord[1];
        mesh.Cell0DsCoordinates(2, id) = Coord[2];
        mesh.Cell0DsID.push_back(id);
        return id;
    }
}

int VerificaEInserisciDual(Vector3d& Coord, map<Vector3d, int, Vector3dComparator>& mappa, PolygonalDual& mesh) {
    Coord = Coord / Coord.norm(); // normalizzazione
    auto it = mappa.find(Coord);
    if (it != mappa.end()) {
        return it->second;
    } else {
        int id = mappa.size();
        mappa[Coord] = id;
        mesh.Cell0DsCoordinates(0, id) = Coord[0];
        mesh.Cell0DsCoordinates(1, id) = Coord[1];
        mesh.Cell0DsCoordinates(2, id) = Coord[2];
        mesh.Cell0DsID.push_back(id);
        return id;
    }
}

array<unsigned int, 3> VerificaEInserisci2(array<unsigned int, 3> NewFace, map<Vector2i, int, Vector2iComparator>& mappa, Polygonal& mesh, array<unsigned int, 3>& Face) {
    Vector2i Lato;
    for( int i = 0; i < 3; i++){
	    int a = NewFace[i % 3];
		int b = NewFace[(i + 1) % 3];
		Lato << std::min(a, b), std::max(a, b);

	    
	    auto it = mappa.find(Lato);
		if (it != mappa.end()) {
			Face[i] = mappa[Lato];
		} else {
			int id = mappa.size();
			mappa[Lato] = id;
			mesh.Cell1DsExtrema(0, id) = Lato[0];
			mesh.Cell1DsExtrema(1, id) = Lato[1];
			mesh.Cell1DsID.push_back(id);
			Face[i] = mappa[Lato];
    		}
	    
    }
	return Face;
}

vector<unsigned int> VerificaEInserisci2Dual(vector<unsigned int> NewFace, map<Vector2i, int, Vector2iComparator>& mappa, Eigen::MatrixXi& Cell1DsExtrema, PolygonalDual& mesh, vector<unsigned int>& Face) {
	int n = NewFace.size();
    Vector2i Lato;
    for( int i = 0; i < n; i++){
	    int a = NewFace[i % n];
		int b = NewFace[(i + 1) % n];
		Lato << std::min(a, b), std::max(a, b);

	    
	    auto it = mappa.find(Lato);
		if (it != mappa.end()) {
			Face.insert(Face.begin() + i, mappa[Lato]);
		} else {
			int id = mappa.size();
			mappa[Lato] = id;
			Cell1DsExtrema(0, id) = Lato[0];
			Cell1DsExtrema(1, id) = Lato[1];
			mesh.Cell1DsID.push_back(id);
			Face.insert(Face.begin() + i, mappa[Lato]);
    		}
	    
    }
	return Face;
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
		
	int numFaces = mesh.Cell2DsVertices.size();
	int NumPunti = 0;
	int NumLati = 0;
	NumLati += mesh.NumCell1Ds * n;
	for( int i = 1; i < n; i++){
		NumLati += i * 3 * numFaces;
	}
	
	for( int i = 1; i < n - 1; i ++){
		NumPunti += i * numFaces;
	}
	NumPunti += mesh.Cell0DsID.size();
	NumPunti += mesh.Cell1DsID.size() * (n - 1);
	meshTriang.NumCell0Ds = NumPunti;
	meshTriang.NumCell1Ds = NumLati;
	meshTriang.NumCell2Ds = n * n * numFaces;
	
	meshTriang.Cell0DsID.reserve(meshTriang.NumCell0Ds);
    meshTriang.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, meshTriang.NumCell0Ds);
    meshTriang.Cell1DsID.reserve(meshTriang.NumCell1Ds);
    meshTriang.Cell1DsExtrema = Eigen::MatrixXi(2, meshTriang.NumCell1Ds);
    meshTriang.Cell2DsID.reserve(meshTriang.NumCell2Ds);
    meshTriang.Cell2DsVertices.reserve(meshTriang.NumCell2Ds);
    meshTriang.Cell2DsEdges.reserve(meshTriang.NumCell2Ds);
	
	for ( unsigned int i = 0; i < meshTriang.NumCell2Ds; i++){
		meshTriang.Cell2DsID.push_back(i);
	}
	
	map<Vector3d, int, Vector3dComparator> ControlloPunti;
	map<Vector2i, int, Vector2iComparator> ControlloEdges;
	Vector3d CoordA; //coordinate A
	Vector3d CoordB; //coordinate B
	Vector3d CoordC; //coordinate C
	Vector3d CoordP; //coordinate nuovo punto
	array<unsigned int, 3> NewFace;
	array<unsigned int, 3> Face;
		
	for(int k = 0; k < numFaces; k++){
		
		int numfaccia = k;
		
		int id;
		id = mesh.Cell2DsVertices[numfaccia][0];
		CoordA[0] = mesh.Cell0DsCoordinates(0,id);
		CoordA[1] = mesh.Cell0DsCoordinates(1,id);
		CoordA[2] = mesh.Cell0DsCoordinates(2,id);
		
		id = mesh.Cell2DsVertices[numfaccia][1];
		CoordB[0] = mesh.Cell0DsCoordinates(0,id);
		CoordB[1] = mesh.Cell0DsCoordinates(1,id);
		CoordB[2] = mesh.Cell0DsCoordinates(2,id);
		
		id = mesh.Cell2DsVertices[numfaccia][2];
		CoordC[0] = mesh.Cell0DsCoordinates(0,id);
		CoordC[1] = mesh.Cell0DsCoordinates(1,id);
		CoordC[2] = mesh.Cell0DsCoordinates(2,id);
		
		for(int i = 0; i < n; i++){ //conteggio piani
			for(int j = 0; j < (n - i); j++){
				CoordP = ((double(i)/n) * CoordA) + ((double(n - (i + j))/n) * CoordB) + ((double(j)/n) * CoordC);
				id = VerificaEInserisci(CoordP, ControlloPunti, meshTriang);
				NewFace[0] = id;
				
				CoordP = ((double(i)/n) * CoordA) + ((double(n - (i + j + 1))/n) * CoordB) + ((double(j + 1)/n) * CoordC);
				id = VerificaEInserisci(CoordP, ControlloPunti, meshTriang);
				NewFace[1] = id;
				
				CoordP = ((double(i + 1)/n) * CoordA) + ((double(n - (i + j + 1))/n) * CoordB) + ((double(j)/n) * CoordC);
				id = VerificaEInserisci(CoordP, ControlloPunti, meshTriang);
				NewFace[2] = id;

				meshTriang.Cell2DsVertices.push_back(NewFace);
				Face = VerificaEInserisci2(NewFace, ControlloEdges, meshTriang, Face);
				meshTriang.Cell2DsEdges.push_back(Face);
			}
			
			for(int j = 0; j < (n - (i + 1)); j++){
				CoordP = ((double(i)/n) * CoordA) + ((double(n - (i + j + 1))/n) * CoordB) + ((double(j + 1)/n) * CoordC);
				id = VerificaEInserisci(CoordP, ControlloPunti, meshTriang);
				NewFace[0] = id;
				
				CoordP = ((double(i + 1)/n) * CoordA) + ((double(n - (i + j + 2))/n) * CoordB) + ((double(j + 1)/n) * CoordC);
				id = VerificaEInserisci(CoordP, ControlloPunti, meshTriang);
				NewFace[1] = id;
				
				CoordP = ((double(i + 1)/n) * CoordA) + ((double(n - (i + j + 1))/n) * CoordB) + ((double(j)/n) * CoordC);
				id = VerificaEInserisci(CoordP, ControlloPunti, meshTriang);
				NewFace[2] = id;
				
				meshTriang.Cell2DsVertices.push_back(NewFace);
				Face = VerificaEInserisci2(NewFace, ControlloEdges, meshTriang, Face);
				meshTriang.Cell2DsEdges.push_back(Face);
			}
		}
	}
}

vector<array<unsigned int, 3>> SearchFaces(int id, vector<array<unsigned int, 3>> VerticesF, int length) {
	
	vector<array<unsigned int, 3>> FacceAdiacenti; 
	double err = 1e-16;
	for( int i = 0; i < length; i++) {
		
		array<unsigned int, 3> Face = VerticesF[i];
		int V0 = Face[0];
		int V1 = Face[1];
		int V2 = Face[2];
		if (abs(id-V0) < err || abs(id-V1) < err || abs(id-V2) < err) {
			FacceAdiacenti.push_back(Face);
		}
	}
	return FacceAdiacenti;
}


vector<array<unsigned int, 3>> SortFaces(int id, int IdOrd, vector<array<unsigned int,3>> FacesAd){
	
	double err =1e-16;
	vector<array<unsigned int,3>> SortedFaces;
	int n = FacesAd.size();
	
	for( int i = 0; i < n; i++){
		
		vector<array<unsigned int,3>> VFacciaAd = SearchFaces( IdOrd, FacesAd ,n);
		array<unsigned int, 3> FacciaAd = VFacciaAd[0];
		SortedFaces.push_back(FacciaAd);
		FacesAd.erase(
			remove(FacesAd.begin(), FacesAd.end(), FacciaAd),
			FacesAd.end()
		);	
		
		
		int V0 = FacciaAd[0];
		int V1 = FacciaAd[1];
		int V2 = FacciaAd[2];
		
		if(id != V0 && IdOrd != V0){
			IdOrd = V0;
		}
		else if(id != V1 && IdOrd != V1){
			IdOrd = V1;
		}else{
			IdOrd = V2;
		}
		
		
	}
	return SortedFaces;
	
	
}


	
bool DualTot(Polygonal& meshTriang, PolygonalDual& meshDual) {
	
	map<Vector3d, int, Vector3dComparator> ControlloPunti;
	map<Vector2i, int, Vector2iComparator> ControlloEdges;
	
	Vector3d CoordA; //coordinate A
	Vector3d CoordB; //coordinate B
	Vector3d CoordC; //coordinate C
	Vector3d CoordP; //coordinate nuovo punto
	vector<unsigned int> NewFace;
	vector<unsigned int> Face;
	Eigen::MatrixXi Cell1DsExtremaBig(2, 10000);
		
	vector<array<unsigned int, 3>> VerticesF = meshTriang.Cell2DsVertices;
	
	int numFaces = meshTriang.NumCell0Ds;
	int NumPunti = meshTriang.NumCell2Ds;
	meshDual.NumCell0Ds = NumPunti;
	meshDual.NumCell2Ds = numFaces;
	
	meshDual.Cell0DsID.reserve(meshDual.NumCell0Ds);
    meshDual.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, meshDual.NumCell0Ds);
    meshDual.Cell2DsID.reserve(meshDual.NumCell2Ds);
    meshDual.Cell2DsVertices.reserve(meshDual.NumCell2Ds);
    meshDual.Cell2DsEdges.reserve(meshDual.NumCell2Ds);
    
    for ( unsigned int i = 0; i < meshDual.NumCell2Ds; i++){
		meshDual.Cell2DsID.push_back(i);
	}
	
	for( int i = 0; i < meshTriang.NumCell0Ds; i++){
		NewFace.clear();
		Face.clear();
		int id = i;
		int length = VerticesF.size();
		vector<array<unsigned int, 3>> FacesAd = SearchFaces(id, VerticesF, length);
		
		vector<array<unsigned int, 3>> SortedFaces;
		int IdOrd = FacesAd[1][1];
		if (id != IdOrd) {
			SortedFaces = SortFaces(id, IdOrd, FacesAd);
			
		}else{
			IdOrd = FacesAd[1][2];
			SortedFaces = SortFaces(id, IdOrd, FacesAd);
		}
		
		for( int i = 0; i < SortedFaces.size(); i++){
			
			int idPunto;
			idPunto = SortedFaces[i][0];
			CoordA[0] = meshTriang.Cell0DsCoordinates(0,idPunto);
			CoordA[1] = meshTriang.Cell0DsCoordinates(1,idPunto);
			CoordA[2] = meshTriang.Cell0DsCoordinates(2,idPunto);
			
			idPunto = SortedFaces[i][1];
			CoordB[0] = meshTriang.Cell0DsCoordinates(0,idPunto);
			CoordB[1] = meshTriang.Cell0DsCoordinates(1,idPunto);
			CoordB[2] = meshTriang.Cell0DsCoordinates(2,idPunto);
			
			idPunto = SortedFaces[i][2];
			CoordC[0] = meshTriang.Cell0DsCoordinates(0,idPunto);
			CoordC[1] = meshTriang.Cell0DsCoordinates(1,idPunto);
			CoordC[2] = meshTriang.Cell0DsCoordinates(2,idPunto);
			
			CoordP = (CoordA + CoordB + CoordC) / 3;
			idPunto = VerificaEInserisciDual(CoordP, ControlloPunti, meshDual);
			NewFace.push_back(idPunto);
		}
		meshDual.Cell2DsVertices.push_back(NewFace);
		Face = VerificaEInserisci2Dual(NewFace, ControlloEdges, Cell1DsExtremaBig, meshDual, Face);
		meshDual.Cell2DsEdges.push_back(Face);
		
	}
	int m = meshDual.Cell1DsID.size();
	meshDual.Cell1DsExtrema = Eigen::MatrixXi(2, m);
	meshDual.Cell1DsExtrema = Cell1DsExtremaBig.leftCols(m);
	meshDual.NumCell1Ds = m;
	
	return true;
}
	
	
	
