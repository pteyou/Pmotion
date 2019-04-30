#pragma once
#include <iostream>
#include <assert.h>

class linearInterp {
public:
    linearInterp(std::vector<double> &x, std::vector<double>& y) : 
        size_(x.size()), x_min(x[0]), x_max(x.back()), x_(x), y_(y) {
            for(int i=0; i<size_-1; i++) assert(x_[i+1] > x_[i]);
        }
    double operator()(double x) {
        int idx = locateLeft(x);
        if (idx == size_) return y_[size_-1];
        if(x>= x_max) return y_.back();
        if (x < x_min) return y_[0]; 
        return (y_[idx+1] - y_[idx]) / (x_[idx+1] - x_[idx]) * (x - x_[idx]) + y_[idx];
    }
private:
    int locateLeft(double x) {
        int res = 0;
        while (x_[res] < x) {
            res++;
        }
        return --res;
    }
    int size_;
    double x_min, x_max;
    std::vector<double> x_, y_;
};