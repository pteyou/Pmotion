#pragma once
#include <iostream>
#include <algorithm>
#include <string>
#include <vector>
#include <map>
#include <assert.h>
#include <omp.h>

class gridPoint{
public:
    double x, y;
    double T_;
    int i , j;
    int status;
    double gradX, gradY;
    gridPoint(double x, double y, double T) : x(x), y(y), T_(T), gradX(0), gradY(0) {}
    gridPoint() : x(0), y(0), T_(0), gradX(0), gradY(0) {}
};

class door {
public:
    gridPoint a_, b_;
    door(gridPoint &a, gridPoint &b) : a_(a), b_(b) {}
    ~door(){}
};

class grid {
public:
    int iMax_, jMax_; 
    std::vector<gridPoint> points_;
    std::map<std::pair<int, int>, int> pointsMap_;
    std::map<int, std::pair<int, int> > rPointsMap_;
    std::vector<int> boundary_;
    std::vector<int> exits_;
    std::vector<door> doors_;
    double dx_, dy_;
    double xMin_, yMin_;
    grid() {}
    grid(std::vector<gridPoint> &points, std::vector<int> &boundary, std::vector<int> &exits, 
        double dx, double dy, std::vector<door> doors, double xMin=0.0, double yMin=0) : 
        points_(points), iMax_(0), jMax_(0), boundary_(boundary), exits_(exits), dx_(dx), 
        dy_(dy), doors_(doors), xMin_(xMin), yMin_(yMin) {
        for (int i=0; i<points_.size(); i++) {
            std::pair<int, int> tmp = std::make_pair(points_[i].i, points_[i].j);
            pointsMap_[tmp] = i;
            rPointsMap_[i] = tmp;
            iMax_ = std::max(iMax_, tmp.first);
            jMax_ = std::max(jMax_, tmp.second);
        }
    }
    gridPoint& get(int i, int j) {
        assert(i<=iMax_ && j<=jMax_);
        std::pair<int, int> tmp = std::make_pair(i,j);
        int idx;
        #pragma omp critical
        idx = pointsMap_[tmp];
        return points_[idx];
    }
    gridPoint& get(int i) {
        return points_[i];
    }
};