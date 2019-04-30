#pragma once

#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>
#include <string>
#include "math.h"
#include <random>
#include <chrono>
#include "linearInterp.h"
#include "utilities.h"
// Define Infinite (Using INT_MAX caused overflow problems) 
#define dr 0.5

void placePeople(std::vector<double> &initX, std::vector<double> &initY);
void putVectorInAFile(std::vector<double> &outv, std::string fname);
std::vector<double> getVectorFromFile(std::string fname);
