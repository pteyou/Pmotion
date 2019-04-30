#pragma once
#include <vector>
#include "assert.h"
#include <numeric>
#include "myTypes.h"
#include <string>
#include <fstream>
#include <iterator>

class polygonalWall {
public:
    polygonalWall(std::vector<double> &verticesX, std::vector<double> &verticesY, 
        std::vector<int> &nbSubPoints) : 
        verticesX_(verticesX), verticesY_(verticesY), nbSubPoints_(nbSubPoints) {
        assert((verticesX_.size() == nbSubPoints_.size() + 1) && (verticesX_.size() == verticesY_.size()));
        size_ = std::accumulate(nbSubPoints_.begin(), nbSubPoints_.end(), 0) - verticesX_.size() + 2;
        pointsX.resize(size_); pointsY.resize(size_);
        int cursor = 0;
        for(int i=0; i < verticesX_.size()-1; i++){
            int nbPointsOnEdge = nbSubPoints_[i];
            for(int k=0; k < nbPointsOnEdge-1; k++){
                pointsX[cursor] = verticesX_[i] + (verticesX_[i+1] - verticesX_[i]) * k / (nbPointsOnEdge-1);
                pointsY[cursor] = verticesY_[i] + (verticesY_[i+1] - verticesY_[i]) * k / (nbPointsOnEdge-1);
                cursor++;
            }
        }
        pointsX[cursor] = verticesX_.back();
        pointsY[cursor] = verticesY_.back();
        nbLines_ = nbSubPoints_.size();
    }

    int size(){ return size_; }
    int size() const { return size_; }
    int nbLines() const { return nbLines_; }
    /* Start numbering at 0 */
    std::vector<double> point(int i){
        assert (i>=0 && i<size_);
        std::vector<double> res { pointsX[i], pointsY[i] };
        return res;
    }

    void output(std::string &fnameX, std::string &fnameY) {
        std::ofstream ofsX(fnameX.c_str());
        std::ostream_iterator<double> out_itX (ofsX, "\n");
        std::copy ( pointsX.begin(), pointsX.end(), out_itX );
        ofsX.close();
        std::ofstream ofsY(fnameY.c_str());
        std::ostream_iterator<double> out_itY (ofsY, "\n");
        std::copy ( pointsY.begin(), pointsY.end(), out_itY );
        ofsY.close();
    }

    std::vector<double> getVerticesX() { return verticesX_;}
    std::vector<double> getVerticesY() { return verticesY_;}
    
    friend void defineInitialPositionsMeetingRoom(DVector &, DVector &, polygonalWall &);

    ~polygonalWall() {}

private:
    std::vector<double> verticesX_, verticesY_;
    std::vector<int> nbSubPoints_;
    int size_;
    std::vector<double> pointsX, pointsY;
    int nbLines_;
};

class Wall {
    friend class EulerBase;
    friend class Euler;
    template<class T>
    friend class EulerFMM;
public:
    Wall(std::vector<polygonalWall> &walls ) : walls_(walls) {
        size_ = 0;
        nbLines_ = 0;
        std::vector<int> sizes(walls_.size());
        partialSum_.resize(walls_.size());
        int tmp = 0;
        for (const polygonalWall& it : walls_){
            size_ += it.size();
            sizes[tmp++] = it.size();
            nbLines_ += it.nbLines();
        }
        std::partial_sum(sizes.begin(), sizes.end(), partialSum_.begin());
    }

    int size() { return size_; }
    int nbLines() const { return nbLines_; }

    void addPolygonalWall(const polygonalWall &w){
        walls_.push_back(w);
        partialSum_.push_back(partialSum_.back() + w.size());
        size_ += w.size();
        nbLines_ += w.nbLines();
        storeCoordiantes();
    }

    std::vector<double> point(int i){
        assert (i>=0 && i<size_);
        int iWall, iPoint, idx;
        for(idx=0; idx<partialSum_.size(); idx++){
            if(partialSum_[idx] >= i+1) break;
        }
        iWall = idx;
        iPoint = partialSum_[idx] - i - 1;
        return walls_[iWall].point(iPoint);
    }

    void storeCoordiantes() {
        wallX_.resize(size_);
        wallY_.resize(size_);
        for(int i=0; i<size_; i++) {
            std::vector<double> coords {point(i)};
            wallX_[i] = coords[0];
            wallY_[i] = coords[1];
        }
    }

    void output(std::string wallFname){
        std::ofstream wallf(wallFname.c_str());  
        for(int i=0; i<size_; i++){
            std::vector<double> point{this->point(i)};
            wallf << point[0]  << "\t" << point[1] << std::endl;
        }
        wallf.close();
    }

    friend double collisionDistanceWall(double, double, double, double, double, double, Wall &);
    friend std::vector<int> crossedWallPointsRegularization (const DVector &, const DVector &, double *, 
        double *, const DVector&, const DVector &, double *, double *, DVector *, Wall &wall);
    friend std::vector<int> crossedWallPointsRegul (const DVector &, const DVector &, double *, 
        double *, const DVector&, const DVector &, double *, double *, DVector *, Wall &wall);
    ~Wall(){}

private:
    std::vector<polygonalWall> walls_;
    int size_;
    std::vector<int> partialSum_;
    DVector wallX_, wallY_;
    int nbLines_;
};