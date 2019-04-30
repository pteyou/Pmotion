#pragma once
#include <vector>
#include <functional>
#include "forces.h"
#include "polygonalWall.h"

class totalForces {
public:
    totalForces(std::vector<std::function<std::vector<double>(Arguments &, Wall&)>> &forces) : 
        forces_(forces) 
    {
        nb_Forces_ = forces_.size();
    }

    int size() { return nb_Forces_; }

    std::vector<double> computeTotalForces(Arguments &args, Wall &wall) {
        std::vector<double> result(2);
        for (int i=0; i<nb_Forces_; i++) {
            std::vector<double> tmp {forces_[i](args, wall)};
            result[0] += tmp[0];
            result[1] += tmp[1];
        }
        return result;
    }
private:
    std::vector<std::function<std::vector<double>(Arguments &, Wall&)>> forces_;
    int nb_Forces_;
};