#include <iostream>
#include <fstream>
#include <iomanip>
#include "Polygonal.hpp"
#include "Eigen/Eigen"
#include "Utils.hpp"
#include <iostream>
#include <sstream>
#include <vector>
#include "UCDUtilities.hpp"


using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;


int main(int argc, char *argv[])
{	
	Polygonal Polygon;
	Polygonal PolygTriang;
	PolygonalDual PolygDual;
	
	double err = 1.0e-16;
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
		if ((abs(b) < err and  c >= 1) or (abs(c) < err and  b >= 1)) // Class I
		{	
			//importa i dati
			if(!ImportMesh(Polygon, q)){
				cerr << "Error: file not found" << endl;
				return 1;
			}
			
			//triangolazione poligono
			TriangTotC_1(b, c, Polygon, PolygTriang);
			
			if(abs(q-3) < err) {
				DualTot(PolygTriang, PolygDual);
			}
			
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
	
	Gedim::UCDUtilities utilities;
    {
        utilities.ExportPoints("./Cell0Ds.inp",
                               PolygTriang.Cell0DsCoordinates
                               );
    }

    {
        utilities.ExportSegments("./Cell1Ds.inp",
                                 PolygTriang.Cell0DsCoordinates,
                                 PolygTriang.Cell1DsExtrema
                                 );
    }
	
	
	return 0;
}
		
	