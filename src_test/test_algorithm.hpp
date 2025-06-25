#pragma once

#include <iostream>
#include <vector>

#include <gtest/gtest.h>
#include "Utils.hpp"

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;

using coppia = pair<int, double>;
using grafo =  vector<vector<coppia>>;

namespace SearchLibrary {

TEST(TestSearching, TestSearchPoint)
{
    int id = 5;
    vector<Vector2i> EdgesList = {{0,1}, {0,2}, {0,3}, {0,4}, {1,2}, {2,3}, {3,4}, {4,1}, {1,5}, {2,5}, {3,5}, {4,5}};
    map<Vector2i, double, Vector2iComparator> EdgesLength = {
	    {{0,1}, 1.0},
	    {{0,2}, 1.1},
	    {{0,3}, 1.2},
	    {{0,4}, 1.3},
	    {{1,2}, 1.4},
	    {{2,3}, 1.5},
	    {{3,4}, 1.6},
	    {{4,1}, 1.7},
	    {{1,5}, 1.8},
	    {{2,5}, 1.9},
	    {{3,5}, 2.0},
	    {{4,5}, 2.1}};
    vector<pair<int, double>> SearchedPoints = {
	    {1, 1.8},
	    {2, 1.9},
	    {3, 2.0},
	    {4, 2.1}};
    EXPECT_EQ(SearchPoint(id, EdgesList, EdgesLength), SearchedPoints);
}

TEST(TestSearching, TestSearchEdges)
{
    int id = 5;
    vector<Vector2i> EdgesList = {{0,1}, {0,2}, {0,3}, {0,4}, {1,2}, {2,3}, {3,4}, {4,1}, {1,5}, {2,5}, {3,5}, {4,5}};
    vector<Vector2i> SearchedEdges = {{1,5}, {2,5}, {3,5}, {4,5}};
    EXPECT_EQ(SearchEdges(id, EdgesList), SearchedEdges);
}

TEST(TestSearching, TestSearchFaces)
{
    int id = 5;
    vector<array<unsigned int, 3>> VerticesF = {{0,1,2}, {0,2,3}, {0,3,4}, {0,4,1}, {1,2,5}, {2,3,5}, {3,4,5}, {4,1,5}};
    int length = 8;
    vector<array<unsigned int, 3>> SearchedFaces = {{1,2,5}, {2,3,5}, {3,4,5}, {4,1,5}};
    EXPECT_EQ(SearchFaces(id, VerticesF, length), SearchedFaces);
}

}


namespace SortLibrary {

TEST(TestSorting, TestSortFaces)
{
    int id = 5;
    int idOrd = 1;
	vector<array<unsigned int, 3>> FacceAd = {{1,2,5}, {3,4,5}, {2,3,5}, {4,1,5}};
	vector<array<unsigned int, 3>> SortedFaces = {{1,2,5}, {2,3,5}, {3,4,5}, {4,1,5}};
    EXPECT_EQ(SortFaces(id, idOrd, FacceAd), SortedFaces);
}

}

namespace PathLibrary {

TEST(TestPath, TestShortestPath)
{
   int id1 = 1;
   int id2 = 2;
   PolygonalDual mesh;
   mesh.NumCell0Ds = 10;
   mesh.NumCell1Ds = 24;
   mesh.NumCell2Ds = 16;

   mesh.Cell0DsID = {0,1,2,3,4,5,6,7,8,9};
   mesh.Cell1DsID = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23};
   mesh.Cell2DsID = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15};
   
 
   mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds);
   mesh.Cell1DsExtrema = Eigen::MatrixXi::Zero(2, mesh.NumCell1Ds);
   
   mesh.Cell0DsCoordinates << 1.0, 1.0, -1.0, -1.0, 1.0, 0.0, 0.0, 0.0, -1.0, 0.0, 
							1.0, -1.0, 1.0, -1.0, 0.0, 1.0, 0.0, 0.0, 0.0, -1.0,
							1.0, -1.0, -1.0, 1.0, 0.0, 0.0, -1.0, 1.0, 0.0, 0.0;

   mesh.Cell1DsExtrema << 1, 4, 0, 5, 2, 6, 4, 4, 5, 5, 8, 7, 0, 7, 3, 8, 1, 9, 6, 9, 8, 4, 9, 7,
							4, 0, 5, 2, 6, 1, 6, 5, 6, 8, 7, 5, 7, 3, 8, 2, 9, 3, 9, 8, 6, 9, 7,4;
   
   mesh.Cell2DsVertices = {{1,6,4}, {6,2,5}, {4,6,5}, {0,4,5}, {2,8,5}, {7,8,3}, {5,8,7}, {0,5,7}, {1,6,9}, {6,2,8}, {9,6,8}, {9,8,3},
	   {4,1,9}, {7,9,3}, {4,9,7}, {0,4,7}   
   };
   
   mesh.Cell2DsEdges = {{0,5,6}, {8,4,3}, {6,8,7}, {1,7,2}, {3,15,9}, {10,14,13}, {9,10,11}, {2,11,12}, {16,5,18}, {20,4,15}, {18,20,19}, {17,19,14},
	   {0,16,21}, {22,17,13}, {21,22,23}, {1,23,12}
   };

   vector<unsigned int> Path = {1,6,2}; 
   EXPECT_EQ(ShortestPath(id1, id2, mesh), Path);
}

}
namespace VerifyLibrary {

TEST(TestVerify, TestVerify_Insert)
{
   PolygonalDual mesh;
   Vector3d a = {1.0, 1.0, 1.0};
   a = a/ a.norm();
   Vector3d b = {1.0, -1.0, -1.0};
   b = b/ b.norm();
   Vector3d c = {1.0, 0.0, 0.0};
   c = c/ c.norm();
   
   map<Vector3d, int, Vector3dComparator> mappa = {
	    {a, 0},
	    {b, 1},
        {c, 2}};
	
   Vector3d Coord = {-1.0, 0.0, 1.0};
   Vector3d Coord2 = {1.0, -1.0, -1.0};
   mesh.NumCell0Ds = 10;
   mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds);
   int id = 3;
   int id2 = 1;
   EXPECT_EQ(VerificaEInserisciDual(Coord, mappa, mesh), id);
   EXPECT_EQ(VerificaEInserisciDual(Coord2, mappa, mesh), id2);
}


TEST(TestVerify2, TestVerify_Insert2)
{
   PolygonalDual mesh;
   vector<unsigned int> NewFace = {4, 6, 5};
   
   map<Vector2i, int, Vector2iComparator> mappa = {
	    {{4, 6}, 0},
	    {{5, 6}, 1},
        {{1, 2}, 2}};
		
   mesh.NumCell1Ds = 24;
   mesh.Cell1DsExtrema = Eigen::MatrixXi::Zero(2, mesh.NumCell1Ds);
   vector<unsigned int> Face;
   vector<unsigned int> ExpectedFace = {0, 1, 3};
   EXPECT_EQ(VerificaEInserisci2Dual(NewFace, mappa, mesh.Cell1DsExtrema, mesh, Face), ExpectedFace);
}

}

namespace TriangulationLibrary{
	
TEST(TestTriangulation, TestTriangulationC_1)
{
   int b = 5;
   int c = 0;
   Polygonal mesh;
   Polygonal meshTriang;
   mesh.NumCell0Ds = 4;
   mesh.NumCell1Ds = 6;
   mesh.NumCell2Ds = 4;

   mesh.Cell0DsID = {0,1,2,3};
   mesh.Cell1DsID = {0,1,2,3,4,5};
   mesh.Cell2DsID = {0,1,2,3};
 
   mesh.Cell0DsCoordinates = Eigen::MatrixXd::Zero(3, mesh.NumCell0Ds);
   mesh.Cell1DsExtrema = Eigen::MatrixXi::Zero(2, mesh.NumCell1Ds);
   
   
   
   mesh.Cell0DsCoordinates << 1.0, 1.0, -1.0, -1.0, 
							1.0, -1.0, 1.0, -1.0,
							1.0, -1.0, -1.0, 1.0;
   mesh.Cell1DsExtrema << 0, 0, 0, 1, 1, 2,
                          1, 2, 3, 2, 3, 3;
						  
   mesh.Cell2DsVertices = {{0,1,2}, {1,2,3}, {0,2,3}, {0,1,3}};
   
   mesh.Cell2DsEdges = {{0,3,1}, {3,5,4}, {1,5,2}, {0,4,2}};
   
   int ExpectNumCell0Ds = 52;
   int ExpectNumCell1Ds = 150;
   int ExpectNumCell2Ds = 100;
   
   TriangTotC_1(b, c, mesh, meshTriang);
   
   EXPECT_EQ(meshTriang.NumCell0Ds, ExpectNumCell0Ds);
   EXPECT_EQ(meshTriang.NumCell1Ds, ExpectNumCell1Ds);
   EXPECT_EQ(meshTriang.NumCell2Ds, ExpectNumCell2Ds);
	
}
	
} 