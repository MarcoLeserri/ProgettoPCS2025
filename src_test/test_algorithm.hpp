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
    int id = 5;
    int idOrd = 1;
	vector<array<unsigned int, 3>> FacceAd = {{1,2,5}, {3,4,5}, {2,3,5}, {4,1,5}};
	vector<array<unsigned int, 3>> SortedFaces = {{1,2,5}, {2,3,5}, {3,4,5}, {4,1,5}};
    EXPECT_EQ(SortFaces(id, idOrd, FacceAd), SortedFaces);
}

}