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
			FileTxt(PolygTriang);
			
			if(abs(q-3) < err) {
				DualTot(PolygTriang, PolygDual);
			}
			
			if(abs(argc - 7) < err){
				int id1 = stoi(argv[5]);
				int id2 = stoi(argv[6]);
				vector<unsigned int> percorso = ShortestPath(id1, id2, PolygDual);
				vector<double> ParaViewPunti = ParaviewPoints(percorso, PolygDual);
				vector<double> ParaViewEdges = ParaviewEdges(percorso, PolygDual);
				
				Gedim::UCDUtilities utilities;
	
				vector<Gedim::UCDProperty<double>> propsVertices(1);
				propsVertices[0].Label = "ShortPath";
				propsVertices[0].UnitLabel = "-";
				propsVertices[0].NumComponents = 1;
				propsVertices[0].Data = ParaViewPunti.data();
			
			
				vector<Gedim::UCDProperty<double>> propsEdges(1);
				propsEdges[0].Label = "ShortPath";
				propsEdges[0].UnitLabel = "-";
				propsEdges[0].NumComponents = 1;
				propsEdges[0].Data = ParaViewEdges.data();
				
				{
					utilities.ExportPoints("./Cell0Ds.inp",
										   PolygDual.Cell0DsCoordinates,
										   propsVertices
										   );
				}
			
				{
					utilities.ExportSegments("./Cell1Ds.inp",
											 PolygDual.Cell0DsCoordinates,
											 PolygDual.Cell1DsExtrema,
											 propsVertices,  // richiesto anche per segmenti
											 propsEdges
											 );
				}

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
	
	/*Gedim::UCDUtilities utilities;
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
    }*/
	
	return 0;
}
		
	