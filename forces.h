#pragma once
#include <vector>
#include "myTypes.h"
#include "polygonalWall.h"

#define epss 1.0e-6

std::vector<double> desiredSpeedMotivationSFM (Arguments &args, Wall &wall);
std::vector<double> desiredSpeedMotivationHeuristics(Arguments &args, Wall &wall);
std::vector<double> wallRepulsionSFM (Arguments &args, Wall &wall);
std::vector<double> wallRepulsionHeuristics (Arguments &args, Wall &wall);
std::vector<double> individualRepulsionCrowdTurbulence (Arguments &args, Wall &wall);
std::vector<double> individualRepulsionSFM (Arguments &args, Wall &wall);
std::vector<double> individualRepulsionHeuristics(Arguments &args, Wall &wall);