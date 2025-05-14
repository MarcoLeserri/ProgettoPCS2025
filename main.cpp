#include <iostream>
#include <fstream>
#include <iomanip>
#include "Eigen/Eigen"
#include "PolygonalMesh.hpp"
#include "Utils.hpp"
#include "UCDUtilities.hpp"
#include <iostream>
#include <sstream>
#include <vector>

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;


int main(int argc, char *argv[])
{	
	double err = 1.0e-16
	if (argc < 5)
	{
		cerr << "Error: insert at least four values" << endl; 
		return 1; 
	}
	
	int p = stoi(argv[1]);
	int q = stoi(argv[2]);
	int b = stoi(argv[3]);
	int c = stoi(argv[4]);
	
	if (abs(p-3) < err and  q >= 3)
	{
		if ((abs(b) < err and  c >= 1) or  (abs(c) < err and  b >= 1)) // Class I
		{
			Polygonal Polygon;
			Polygonal PolygTriang;
			
			//importa i dati
			if(!ImportMesh(Polygon)){
				cerr << "Error: file not found" << endl;
				return 1;
			}
			
			//triangolazione poligono
			TriangTotC_1(b, c, Polygon, PolygTriang);
		
		
		
		
		
		
		
		
		}
		else if(abs(b-c) < err) {
		// polyhedrons of II Class 
		}
		else {
			cerr << "Error: polyhedron belongs to III class " << endl;			
		}	
	}
	else {
		cerr << "Error: values out of range" << endl;
	}
		
	



	Polygonal mesh