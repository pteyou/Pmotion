#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <iterator>
#include "FMM/grid.h"
#define INF 100000 
#define epsi 1.0e-6

std::vector<double> getVectorFromFile(std::string &fname);

void putVectorInAFile(std::vector<double> &outv, std::string &fname);

bool onSegment(gridPoint p, gridPoint q, gridPoint r);
int orientation(gridPoint p, gridPoint q, gridPoint r);
bool doIntersect(gridPoint p1, gridPoint q1, gridPoint p2, gridPoint q2);
bool isInside(gridPoint* polygon, int n, gridPoint p);
