#pragma once
#include "grid.h"
#include "math.h"
#include <fstream>
#include "rootFinding.h"
#include "../cases.h"

#define INFI 1.0e8
#define defaultT 1000000
#define guess 10000.0
#define nbExternalBoundaryExtraLayers 2
#define nbInternalBoundaryExtraLayers 1
#define outFluxLimiter 0


auto pointCompare = [](gridPoint &a, gridPoint &b) -> bool {return a.T_ > b.T_;};

class FMM {
private:
    grid &g_;
    std::vector<int> narrowBand_;
    double F_;
public:
    FMM(grid &g) : g_(g) {
        //set T=infinity to all points on boundaries
        for(auto i : g_.boundary_) {
            g_.get(i).T_ = INFI;
            g_.get(i).status = 0;
            //narrowBand_.push_back(i);
        }
        // set T=0 on all exits points
        for(auto i : g_.exits_) { 
            g_.get(i).T_ = 0.0;
            g_.get(i).status = 0;
            narrowBand_.push_back(i);
        }
        // set other values to default T
        for(int i=0; i < g_.points_.size(); i++) {
            if (std::find(g_.exits_.begin(), g_.exits_.end(), i) == g.exits_.end() &&
                std::find(g_.boundary_.begin(), g_.boundary_.end(), i) == g_.boundary_.end()) {
                g_.get(i).T_ = defaultT;
                g_.get(i).status = -1;
            }
        }
        F_ = 1.0;
    }

    void step() {
        // update the band
        std::vector<gridPoint> bandMinHeap;
        for(auto i : narrowBand_) {
            bandMinHeap.push_back(g_.get(i));
        }
        std::make_heap(bandMinHeap.begin(), bandMinHeap.end(), pointCompare);
        gridPoint pivot = bandMinHeap.front();
        g_.get(pivot.i, pivot.j).status = 0;
        std::pair<int, int> tmp {pivot.i, pivot.j};
        narrowBand_.erase(std::remove(narrowBand_.begin(), narrowBand_.end(), 
            g_.pointsMap_[tmp]), narrowBand_.end());
        //left
        if (g_.get((pivot.i - 1) % (g_.iMax_ + 1), pivot.j).status != 0) {
            if (g_.get((pivot.i - 1) % (g_.iMax_ + 1), pivot.j).status == -1)
                narrowBand_.push_back(g_.pointsMap_[std::make_pair((pivot.i - 1) % (g_.iMax_ + 1), pivot.j)]);
            g_.get((pivot.i - 1) % (g_.iMax_ + 1), pivot.j).status = 1;
            /* double a = std::min(g_.get((pivot.i - 2) % (g_.iMax_ + 1), pivot.j).T_, 
                g_.get((pivot.i), pivot.j).T_);
            double b = std::min(g_.get((pivot.i - 1) % (g_.iMax_ + 1), (pivot.j - 1) % (g_.jMax_ + 1)).T_, 
                g_.get((pivot.i - 1) % (g_.iMax_ + 1) , (pivot.j + 1) % (g_.jMax_+1)).T_);
            if (F_ > std::abs(a - b)) 
                g_.get((pivot.i - 1) % (g_.iMax_ + 1), pivot.j).T_ = 0.5 * (a+b+sqrt(2.0*F_*F_ - (a-b) * (a-b)));
            else 
                g_.get((pivot.i - 1) % (g_.iMax_ + 1), pivot.j).T_ = F_ + std::min(a,b);*/
            g_.get((pivot.i - 1) % (g_.iMax_ + 1), pivot.j).T_ = solver(guess, g_, (pivot.i - 1) % (g_.iMax_ + 1), pivot.j);
        }
        //right
        if (g_.get((pivot.i + 1) % (g_.iMax_ + 1), pivot.j).status != 0) {
            if (g_.get((pivot.i + 1) % (g_.iMax_ + 1), pivot.j).status == -1)
                narrowBand_.push_back(g_.pointsMap_[std::make_pair((pivot.i + 1) % (g_.iMax_ + 1), pivot.j)]);
            g_.get((pivot.i + 1) % (g_.iMax_ + 1), pivot.j).status = 1;
            /* double a = std::min(g_.get(pivot.i, pivot.j).T_, 
                g_.get((pivot.i+2) % (g_.iMax_ + 1), pivot.j).T_);
            double b = std::min(g_.get((pivot.i+1) % (g_.iMax_ + 1), (pivot.j - 1) % (g_.jMax_ + 1)).T_, 
                g_.get((pivot.i+1) % (g_.iMax_ + 1), (pivot.j + 1) % (g_.jMax_ + 1)).T_);
            if (F_ > std::abs(a - b)) 
                g_.get((pivot.i + 1) % (g_.iMax_ + 1), pivot.j).T_ = 0.5 * (a+b+sqrt(2.0*F_*F_ - (a-b) * (a-b)));
            else 
                g_.get((pivot.i + 1) % (g_.iMax_ + 1), pivot.j).T_ = F_ + std::min(a,b); */
                g_.get((pivot.i + 1) % (g_.iMax_ + 1), pivot.j).T_ = solver(guess, g_, (pivot.i + 1) % (g_.iMax_ + 1), pivot.j);
        }
        //up
        if (g_.get(pivot.i, (pivot.j + 1) % (g_.jMax_ + 1)).status != 0) {
            if (g_.get(pivot.i, (pivot.j + 1) % (g_.jMax_ + 1)).status == -1)
                narrowBand_.push_back(g_.pointsMap_[std::make_pair(pivot.i, (pivot.j + 1) % (g_.jMax_ + 1))]);
            g_.get(pivot.i, (pivot.j + 1) % (g_.jMax_ + 1)).status = 1;
            /* double a = std::min(g_.get((pivot.i - 1) % (g_.iMa=x_ + 1), (pivot.j + 1) % (g_.jMax_ + 1)).T_, 
                g_.get((pivot.i + 1) % (g_.iMax_ + 1), (pivot.j + 1) % (g_.jMax_ + 1)).T_);
            double b = std::min(g_.get(pivot.i, pivot.j).T_, 
                g_.get(pivot.i, (pivot.j + 2) % (g_.jMax_ + 1)).T_);
            if (F_ > std::abs(a - b)) 
                g_.get(pivot.i, (pivot.j + 1) % (g_.jMax_ + 1)).T_ = 0.5 * (a+b+sqrt(2.0*F_*F_ - (a-b) * (a-b)));
            else 
                g_.get(pivot.i, (pivot.j + 1) % (g_.jMax_ + 1)).T_ = F_ + std::min(a,b); */
                g_.get(pivot.i, (pivot.j + 1) % (g_.jMax_ + 1)).T_ = solver(guess, g_, pivot.i, (pivot.j + 1) % (g_.jMax_ + 1)); 
        }
        //down
        if (g_.get(pivot.i, (pivot.j - 1) % (g_.jMax_ + 1)).status != 0) {
            if (g_.get(pivot.i, (pivot.j - 1) % (g_.jMax_ + 1)).status == -1)
                narrowBand_.push_back(g_.pointsMap_[std::make_pair(pivot.i, (pivot.j - 1) % (g_.jMax_ + 1))]);
            g_.get(pivot.i, (pivot.j - 1) % (g_.jMax_ + 1)).status = 1;
            /* double a = std::min(g_.get((pivot.i - 1) % (g_.iMax_ + 1), (pivot.j - 1) % (g_.jMax_ + 1)).T_, 
                g_.get((pivot.i + 1) % (g_.iMax_ + 1), (pivot.j - 1) % (g_.jMax_ + 1)).T_);
            double b = std::min(g_.get(pivot.i, (pivot.j - 2) % (g_.jMax_ + 1)).T_, 
                g_.get(pivot.i, pivot.j).T_);
            if (F_ > std::abs(a - b)) 
                g_.get(pivot.i, (pivot.j - 1) % (g_.jMax_ + 1)).T_ = 0.5 * (a+b+sqrt(2.0*F_*F_ - (a-b) * (a-b)));
            else 
                g_.get(pivot.i, (pivot.j - 1) % (g_.jMax_ + 1)).T_ = F_ + std::min(a,b); */
            g_.get(pivot.i, (pivot.j - 1) % (g_.jMax_ + 1)).T_ = solver(guess, g_, pivot.i, (pivot.j - 1) % (g_.jMax_ + 1));
        }

    }

    void run() {
        while(!narrowBand_.empty()){
            step();
        } 
        // store the disired directions
        for(auto t : g_.points_) {
            if (std::find(g_.boundary_.begin(), g_.boundary_.end(), g_.pointsMap_[std::make_pair(t.i, t.j)]) == g_.boundary_.end()
             && std::find(g_.exits_.begin(), g_.exits_.end(), g_.pointsMap_[std::make_pair(t.i, t.j)]) == g_.exits_.end()) {
                 double right = g_.get(t.i + 1, t.j).T_;
                 double left = g_.get(t.i - 1, t.j).T_;
                 double up = g_.get(t.i, t.j + 1).T_;
                 double down = g_.get(t.i, t.j - 1).T_;
                 double gradX = (right - left) / (2.0 * g_.dx_);
                 double gradY = (up - down) / (2.0 * g_.dy_);
                 double norm = sqrt(pow(gradX, 2.0) + pow(gradY, 2.0));
                 g_.get(t.i, t.j).gradX = -gradX / norm;
                 g_.get(t.i, t.j).gradY = -gradY / norm;
             }
        }
    }

    void output(std::string &fname, std::string &vfname, std::string &bname, std::string &ename) {
        std::ofstream ofs(fname.c_str());
        for(auto t : g_.points_) {
            if (std::find(g_.boundary_.begin(), g_.boundary_.end(), g_.pointsMap_[std::make_pair(t.i, t.j)]) == g_.boundary_.end())
                ofs << t.x << "\t" << t.y << "\t" << t.T_ << std::endl;
        }
        ofs.close();
        ofs.open(vfname.c_str());
        for(auto t : g_.points_) {
            if (std::find(g_.boundary_.begin(), g_.boundary_.end(), g_.pointsMap_[std::make_pair(t.i, t.j)]) == g_.boundary_.end()
             && std::find(g_.exits_.begin(), g_.exits_.end(), g_.pointsMap_[std::make_pair(t.i, t.j)]) == g_.exits_.end()) {
                ofs << t.x << "\t" << t.y << "\t" << t.gradX << "\t" << t.gradY << std::endl;
             }
        }
        ofs.close();
        ofs.open(bname.c_str());
        for (auto t : g_.boundary_) {
            ofs << g_.get(t).x << "\t" << g_.get(t).y << std::endl;
        }
        ofs.close();

        ofs.open(ename.c_str());
        for (auto t : g_.exits_) {
            ofs << g_.get(t).x << "\t" << g_.get(t).y << std::endl;
        }
        ofs.close();
    }

};
