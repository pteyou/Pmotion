#include <iostream>
#include <algorithm>
#include <vector>
#include <string>
#include <fstream>
#include "../polygonalWall.h"
#include "FMM.h"
#include "../utilities.h"
#include <set>



void addADoorAsExit(door &d, grid &g, std::vector<int> &exits){
    std::vector<double> X{d.a_.x, d.b_.x}, Y{d.a_.y, d.b_.y};
    std::vector<int> nbPoints{g.iMax_};
    polygonalWall porte(X, Y, nbPoints);
    std::set<int> newIdx;
    for(int i=0; i<porte.size(); i++){
        std::vector<double> xy {porte.point(i)};
        // locate the rectangle containing this point and mark its vertices as boundaries
        double x{xy[0]}, y{xy[1]};
        int kx = (int) std::floor((x-g.xMin_) / g.dx_);
        int ky = (int) std::floor((y-g.yMin_) / g.dy_);
        if (abs(x - (double) kx * g.dx_) < epsi) {
            // mark only one point
            newIdx.insert(g.pointsMap_[std::make_pair(kx, ky)]);
        }
        else {
            // mark the two surounding points
            newIdx.insert(g.pointsMap_[std::make_pair(kx, ky)]);
            newIdx.insert(g.pointsMap_[std::make_pair(kx+1, ky)]);
        }
        if (abs(y - (double) ky * g.dy_) < epsi) {
            // mark only one point
            newIdx.insert(g.pointsMap_[std::make_pair(kx, ky)]);
        }
        else {
            // mark the two surounding points
            newIdx.insert(g.pointsMap_[std::make_pair(kx, ky)]);
            newIdx.insert(g.pointsMap_[std::make_pair(kx, ky+1)]);
        }
    }
    exits.insert(exits.end(), newIdx.begin(), newIdx.end());
}

void placeAPolygonOnBoundary(polygonalWall &w, grid &g, 
    std::vector<int> &boundary, bool internal = false) {
    std::set<int> newIdx;
    std::vector<double> X{w.getVerticesX()};
    std::vector<double> Y{w.getVerticesY()};
    std::vector<int> nbPoints{g.iMax_};
    for (int idx=0; idx < X.size()-1; idx++){
        std::vector<double> wX{X[idx], X[idx+1]}, wY{Y[idx], Y[idx+1]};
        polygonalWall wall(wX, wY, nbPoints);
        for(int i=0; i<wall.size(); i++){
            std::vector<double> xy {wall.point(i)};
            // locate the rectangle containing this point and mark its vertices as boundaries
            double x{xy[0]}, y{xy[1]};
            int kx = (int) std::floor((x-g.xMin_) / g.dx_);
            int ky = (int) std::floor((y-g.yMin_) / g.dy_);
            if (abs(x - (double) kx * g.dx_) < epsi) {
                // mark only one point
                newIdx.insert(g.pointsMap_[std::make_pair(kx, ky)]);
            }
            else {
                // mark the two surounding points
                newIdx.insert(g.pointsMap_[std::make_pair(kx, ky)]);
                newIdx.insert(g.pointsMap_[std::make_pair(kx+1, ky)]);
            }
            if (abs(y - (double) ky * g.dy_) < epsi) {
                // mark only one point
                newIdx.insert(g.pointsMap_[std::make_pair(kx, ky)]);
            }
            else {
                // mark the two surounding points
                newIdx.insert(g.pointsMap_[std::make_pair(kx, ky)]);
                newIdx.insert(g.pointsMap_[std::make_pair(kx, ky+1)]);
            }
        }
    }
    boundary.insert(boundary.end(), newIdx.begin(), newIdx.end());
}

void initialiseGrid(std::vector<gridPoint> &gridPoints, double xmin, double xmax, 
    double ymin, double ymax, double dx, double dy) {
    assert(gridPoints.empty());
    int i=0, j=0; 
    for(double x=xmin; x<=xmax; x+=dx) {
        j=0;
        for(double y=ymin; y<=ymax; y+=dy) {
            gridPoint a (x,y,0);
            a.status = -1;
            a.i = i;
            a.j = j;
            gridPoints.push_back(a);
            j++;
        }
        i++;
    }
}

void reduceDoorsSize(std::vector<door> &doors, Wall &wall, std::vector<int> &exits, 
    grid &g, std::vector<int>& boundary, double outFluxLimiter_){
    if (outFluxLimiter > 0) {
        exits.clear();
        std::vector<int> nbSubPoints {10};
        for (int i=0; i<doors.size(); i++){
            double dx {doors[i].b_.x - doors[i].a_.x};
            std::vector<double> X1{doors[i].a_.x, doors[i].a_.x + outFluxLimiter_ / 2.0 * dx}, 
                                X2{doors[i].b_.x - outFluxLimiter_ / 2.0 * dx, doors[i].b_.x};
            double dy {doors[i].b_.y - doors[i].a_.y};
            std::vector<double> Y1{doors[i].a_.y, doors[i].a_.y + outFluxLimiter_ / 2.0 * dy}, 
                                Y2{doors[i].b_.y - outFluxLimiter_ / 2.0 * dy, doors[i].b_.y};
            polygonalWall w1(X1, Y1, nbSubPoints), w2(X2, Y2, nbSubPoints);
            wall.addPolygonalWall(w1);
            placeAPolygonOnBoundary(w1, g, boundary);
            wall.addPolygonalWall(w2);
            placeAPolygonOnBoundary(w2, g, boundary);
            gridPoint a(X1[1], Y1[1], 0.0), b(X2[0], Y2[0], 0.0);
            door d(a, b);
            addADoorAsExit(d, g, exits);
        }
    }
}

std::ostream & operator << (std::ostream &os, gridPoint & p) {
    os << p.T_ << std::endl;
    return os;
}

std::ofstream & operator << (std::ofstream &ofs, std::vector<gridPoint> &gp) {
    for (auto a : gp) {
        ofs << a.x << '\t' << a.y << std::endl;
    }
    return ofs;
}


grid Reunion::initializeGrid(Wall &wall, double dx, double dy) {
    std::string xName {"../inputs/tableX.dat"}, yName {"../inputs/tableY.dat"};
    std::vector<double> tableX {getVectorFromFile(xName)};
    std::vector<double> tableY {getVectorFromFile(yName)};
    tableX.push_back(tableX[0]);
    tableY.push_back(tableY[0]);
    int nbTableLines = tableX.size() - 1;
    std::vector<int> nbTableSubPoints(nbTableLines, 4);
    polygonalWall table(tableX, tableY, nbTableSubPoints);

    // to do :: add the door
    std::vector<double> roomX {0.0, 0.0, 8.24, 8.24, 0.0, 0.0};
    std::vector<double> roomY {0.56, 0.0, 0.0, 4.32, 4.32, 1.36};
    std::vector<int> nbRoomSubPoints{10, 50, 20, 50, 30};

    polygonalWall room(roomX, roomY, nbRoomSubPoints);

    std::vector<polygonalWall> wallVector {table, room};
    Wall w(wallVector);

    // computing the room's grid
    double xmin = *std::min_element(roomX.begin(), roomX.end());
    double xmax = *std::max_element(roomX.begin(), roomX.end());
    double ymin = *std::min_element(roomY.begin(), roomY.end());
    double ymax = *std::max_element(roomY.begin(), roomY.end());

    std::vector<gridPoint> gridPoints;
    initialiseGrid(gridPoints, xmin, xmax, ymin, ymax, dx, dy);
    std::vector<int> boundary, exits;
    std::vector<door> doors;

    grid g(gridPoints, boundary, exits, dx, dy, doors);
    for(int i=0; i<=g.iMax_; i++) {
        boundary.push_back(g.pointsMap_[std::make_pair(i, 0)]);
    }
    for(int i=0; i<=g.iMax_; i++) {
        boundary.push_back(g.pointsMap_[std::make_pair(i, g.jMax_)]);
    }
    for(int j=0; j<=g.jMax_; j++) {
        boundary.push_back(g.pointsMap_[std::make_pair(0, j)]);
    }
    for(int j=0; j<=g.jMax_; j++) {
        boundary.push_back(g.pointsMap_[std::make_pair(g.iMax_, j)]);
    }
    // defining table as boundary
    for(int i=0; i<table.size(); i++){
        std::vector<double> xy {table.point(i)};
        // locate the rectangle containing this point and mark its vertices as boundaries
        double x{xy[0]}, y{xy[1]};
        int kx = (int) std::floor(x / g.dx_);
        int ky = (int) std::floor(y / g.dy_);
        if (abs(x - (double) kx * dx) < epsi) {
            // mark only one point
            boundary.push_back(g.pointsMap_[std::make_pair(kx, ky)]);
        }
        else {
            // mark the two surounding points
            boundary.push_back(g.pointsMap_[std::make_pair(kx, ky)]);
            boundary.push_back(g.pointsMap_[std::make_pair(kx+1, ky)]);
        }
        if (abs(y - (double) ky * dy) < epsi) {
            // mark only one point
            boundary.push_back(g.pointsMap_[std::make_pair(kx, ky)]);
        }
        else {
            // mark the two surounding points
            boundary.push_back(g.pointsMap_[std::make_pair(kx, ky)]);
            boundary.push_back(g.pointsMap_[std::make_pair(kx, ky+1)]);
        }
    }

    //defining exits
    for(int j=0; j<g.jMax_; j++){
        if(g.get(0, j).y >= 0.56 && g.get(0, j).y <= 1.36)
            exits.push_back(g.pointsMap_[std::make_pair(0, j)]);
    }

    //defining doors
    gridPoint a {g.get(*std::min_element(exits.begin(), exits.end()))};
    gridPoint b {g.get(*std::max_element(exits.begin(), exits.end()))};
    doors.push_back(door(a,b));

    reduceDoorsSize(doors, wall, exits, g, boundary, outFluxLimiter);

    grid gg(gridPoints, boundary, exits, dx, dy, doors);

    FMM fmm(gg);

    fmm.run();

    std::string test("reunion.dat");
    std::string vtest("vreunion.dat");
    std::string btest("breunion.dat");
    std::string etest("ereunion.dat");
    fmm.output(test, vtest, btest, etest);
    std::cout << " ******* FMM completed *****" << std::endl;
    
    return gg;
}

grid Cinema::initializeGrid(Wall &wall, double dx, double dy) {
    double xmin{0}, ymin{0}, xmax{27.9}, ymax{16.0};
    std::vector<gridPoint> gridPoints;
    initialiseGrid(gridPoints, xmin, xmax, ymin, ymax, dx, dy);
    std::vector<int> boundary, exits;
    std::vector<door> doors;

    grid g(gridPoints, boundary, exits, dx, dy, doors);

    // defining exits and doors
    int iLeft = 0, iRight = 0, jUp = 0, jDown = 0;
    bool foundLeft {false}, foundDown{false};
    for(int i=0; i<g.iMax_; i++){
        if(g.get(i, 0).x >= 2.7 && g.get(i, 0).x <= 4.7) {
            exits.push_back(g.pointsMap_[std::make_pair(i, 0)]);
            g.get(i, 0).gradX = 0.0;
            g.get(i, 0).gradY = -1.0;
            if (!foundLeft) {iLeft = i ; foundLeft = true; }
            iRight = std::max(i, iRight);
        }
    }
    gridPoint a {g.get(iLeft, 0)};
    gridPoint b {g.get(iRight, 0)};
    doors.push_back(door(a,b));

    iLeft=0; iRight=0; foundLeft=false;
    for(int i=0; i<g.iMax_; i++){
        if(g.get(i, 0).x >= 21 && g.get(i, 0).x <= 23.4) {
            exits.push_back(g.pointsMap_[std::make_pair(i, 0)]);
            g.get(i, 0).gradX = 0.0;
            g.get(i, 0).gradY = -1.0;
            if (!foundLeft) {iLeft = i ; foundLeft = true; }
            iRight = std::max(i, iRight);
        }
    }
    gridPoint aa {g.get(iLeft, 0)};
    gridPoint ba {g.get(iRight, 0)};
    doors.push_back(door(aa,ba));

    iLeft=0; iRight=0; foundLeft=false;
    for(int i=0; i<g.iMax_; i++){
        if(g.get(i, g.jMax_).x >= 21 && g.get(i, g.jMax_).x <= 23.4) {
            exits.push_back(g.pointsMap_[std::make_pair(i, g.jMax_)]);
            g.get(i, g.jMax_).gradX = 0.0;
            g.get(i, g.jMax_).gradY = 1.0;
            if (!foundLeft) {iLeft = i ; foundLeft = true; }
            iRight = std::max(i, iRight);
        }
    }
    gridPoint aaa {g.get(iLeft, g.jMax_)};
    gridPoint baa {g.get(iRight, g.jMax_)};
    doors.push_back(door(aaa,baa));

    iLeft=0; iRight=0; foundLeft=false;
    for(int j=0; j<g.jMax_; j++){
        if(g.get(g.iMax_, j).y >= 4.0 && g.get(g.iMax_, j).y <= 5.5) {
            exits.push_back(g.pointsMap_[std::make_pair(g.iMax_, j)]);
            g.get(g.iMax_, j).gradX = 1.0;
            g.get(g.iMax_, j).gradY = 0.0;
            if(!foundDown) {jDown = j; foundDown=true; }
            jUp = std::max(jUp, j);
        }
    }
    gridPoint aaaa {g.get(g.iMax_, jDown)};
    gridPoint baaa {g.get(g.iMax_, jUp)};
    doors.push_back(door(aaaa,baaa));

    jDown=0; jUp=0; foundDown=false;
    for(int j=0; j<g.jMax_; j++){
        if(g.get(g.iMax_, j).y >= 10.5 && g.get(g.iMax_, j).y <= 12) {
            exits.push_back(g.pointsMap_[std::make_pair(g.iMax_, j)]);
            g.get(g.iMax_, j).gradX = 1.0;
            g.get(g.iMax_, j).gradY = 0.0;
            if(!foundDown) {jDown = j; foundDown=true; }
            jUp = std::max(jUp, j);
        }
    }
    gridPoint aaaaa {g.get(g.iMax_, jDown)};
    gridPoint baaaa {g.get(g.iMax_, jUp)};
    doors.push_back(door(aaaaa,baaaa));

    // defining external boundaries
    for(int i=0; i<=g.iMax_; i++) {
        if(std::find(exits.begin(), exits.end(), g.pointsMap_[std::make_pair(i, 0)]) == exits.end()) {
            boundary.push_back(g.pointsMap_[std::make_pair(i, 0)]);
            g.get(i, 0).gradX = (i==0  ? 1.0 : (i == g.iMax_ ? -1.0 : 0.0));
            g.get(i, 0).gradY = 1.0;
            for(int idx=0; idx < nbExternalBoundaryExtraLayers; idx++) {
                boundary.push_back(g.pointsMap_[std::make_pair(i, idx)]);
                g.get(i, idx).gradX = (i==0  ? 1.0 : (i == g.iMax_ ? -1.0 : 0.0));
                g.get(i, idx).gradY = 1.0;
            }
        }
    }
    for(int i=0; i<=g.iMax_; i++) {
        if(std::find(exits.begin(), exits.end(), g.pointsMap_[std::make_pair(i, g.jMax_)]) == exits.end()) {
            boundary.push_back(g.pointsMap_[std::make_pair(i, g.jMax_)]);
            g.get(i, g.jMax_).gradX = (i==0  ? 1.0 : (i == g.iMax_ ? -1.0 : 0.0));
            g.get(i, g.jMax_).gradY = -1.0;
            for(int idx=0; idx < nbExternalBoundaryExtraLayers; idx++) {
                boundary.push_back(g.pointsMap_[std::make_pair(i, g.jMax_ - idx)]);
                g.get(i, g.jMax_ - idx).gradX = (i==0  ? 1.0 : (i == g.iMax_ ? -1.0 : 0.0));
                g.get(i, g.jMax_ - idx).gradY = -1.0;
            }
        }
    }
    for(int j=0; j<=g.jMax_; j++) {
        boundary.push_back(g.pointsMap_[std::make_pair(0, j)]);
        g.get(0, j).gradX = 1.0;
        g.get(0, j).gradY = (j==0  ? 1.0 : (j == g.jMax_ ? -1.0 : 0.0));
        for(int idx=0; idx < nbExternalBoundaryExtraLayers; idx++) {
            boundary.push_back(g.pointsMap_[std::make_pair(idx, j)]);
            g.get(idx, j).gradX = 1.0;
            g.get(idx, j).gradY = (j==0  ? 1.0 : (j == g.jMax_ ? -1.0 : 0.0));
        }
    }
    for(int j=0; j<=g.jMax_; j++) {
        if(std::find(exits.begin(), exits.end(), g.pointsMap_[std::make_pair(g.iMax_, j)]) == exits.end()) {
            boundary.push_back(g.pointsMap_[std::make_pair(g.iMax_, j)]);
            g.get(g.iMax_, j).gradX = -1.0;
            g.get(g.iMax_, j).gradY = (j==0  ? 1.0 : (j == g.jMax_ ? -1.0 : 0.0));
            for(int idx=0; idx < nbExternalBoundaryExtraLayers; idx++) {
                boundary.push_back(g.pointsMap_[std::make_pair(g.iMax_ - idx, j)]);
                g.get(g.iMax_ - idx, j).gradX = -1.0;
                g.get(g.iMax_ - idx, j).gradY = (j==0  ? 1.0 : (j == g.jMax_ ? -1.0 : 0.0));
            }
        }
    }
    // adding internal boundaries
    std::set<int> bound;
    std::vector<double> rightUpX{17.4, 23.4, 23.4};
    std::vector<double> rightUpY{1.2, 1.2, 0.0};
    std::vector<int> nbPointsRightUp{g.iMax_, g.jMax_};
    polygonalWall rightUp(rightUpX, rightUpY, nbPointsRightUp);
    for(int i=0; i<rightUp.size(); i++){
        std::vector<double> xy {rightUp.point(i)};
        // locate the rectangle containing this point and mark its vertices as boundaries
        double x{xy[0]}, y{xy[1]};
        int kx = (int) std::floor(x / g.dx_);
        int ky = (int) std::floor(y / g.dy_);
        if (abs(x - (double) kx * dx) < epsi) {
            // mark only one point
            bound.insert(g.pointsMap_[std::make_pair(kx, ky)]);
        }
        else {
            // mark the two surounding points
            bound.insert(g.pointsMap_[std::make_pair(kx, ky)]);
            bound.insert(g.pointsMap_[std::make_pair(kx+1, ky)]);
        }
        if (abs(y - (double) ky * dy) < epsi) {
            // mark only one point
            bound.insert(g.pointsMap_[std::make_pair(kx, ky)]);
        }
        else {
            // mark the two surounding points
            bound.insert(g.pointsMap_[std::make_pair(kx, ky)]);
            bound.insert(g.pointsMap_[std::make_pair(kx, ky+1)]);
        }
    }

    std::vector<double> leftUpX{17.4, 23.4, 23.4};
    std::vector<double> leftUpY{14.8, 14.8, 16};
    std::vector<int> nbPointsLeftUp{g.iMax_, g.jMax_};
    polygonalWall leftUp(leftUpX, leftUpY, nbPointsLeftUp);
    for(int i=0; i<leftUp.size(); i++){
        std::vector<double> xy {leftUp.point(i)};
        // locate the rectangle containing this point and mark its vertices as boundaries
        double x{xy[0]}, y{xy[1]};
        int kx = (int) std::floor(x / g.dx_);
        int ky = (int) std::floor(y / g.dy_);
        if (abs(x - (double) kx * dx) < epsi) {
            // mark only one point
            bound.insert(g.pointsMap_[std::make_pair(kx, ky)]);
        }
        else {
            // mark the two surounding points
            bound.insert(g.pointsMap_[std::make_pair(kx, ky)]);
            bound.insert(g.pointsMap_[std::make_pair(kx+1, ky)]);
        }
        if (abs(y - (double) ky * dy) < epsi) {
            // mark only one point
            bound.insert(g.pointsMap_[std::make_pair(kx, ky)]);
        }
        else {
            // mark the two surounding points
            bound.insert(g.pointsMap_[std::make_pair(kx, ky)]);
            bound.insert(g.pointsMap_[std::make_pair(kx, ky+1)]);
        }
    }

    std::string rWallFname{"../inputs/cinemaRowWalls.dat"};
    std::vector<double> rwallpositions {getVectorFromFile(rWallFname)};
    for(int i=0; i<rwallpositions.size(); i+=60){
        // locate all the rectangles hiting the line
        std::vector<double> rowX{rwallpositions[i], rwallpositions[i+58]};
        std::vector<double> rowY{rwallpositions[i+1], rwallpositions[i+59]};
        std::vector<int> rowNbPoints{g.jMax_};
        polygonalWall row(rowX, rowY, rowNbPoints);
        for(int i=0; i<row.size(); i++){
            std::vector<double> xy {row.point(i)};
            // locate the rectangle containing this point and mark its vertices as boundaries
            double x{xy[0]}, y{xy[1]};
            int kx = (int) std::floor(x / g.dx_);
            int ky = (int) std::floor(y / g.dy_);
            if (abs(x - (double) kx * dx) < epsi) {
                // mark only one point
                bound.insert(g.pointsMap_[std::make_pair(kx, ky)]);
                g.get(kx, ky).gradX = (i==0 || i==row.size()-1) ? 0.0 : -1.0;
                g.get(kx, ky).gradY = (i==0 ? -1.0 : (i==row.size()-1 ? 1.0 : 0.0));

                for(int idx=0; idx<nbInternalBoundaryExtraLayers; idx++) {
                    bound.insert(g.pointsMap_[std::make_pair(kx-idx, ky)]);
                    g.get(kx-idx, ky).gradX = (i==0 || i==row.size()-1) ? 0.0 : -1.0;
                    g.get(kx-idx, ky).gradY = (i==0 ? -1.0 : (i==row.size()-1 ? 1.0 : 0.0));
                    bound.insert(g.pointsMap_[std::make_pair(kx+idx, ky)]);
                    g.get(kx+idx, ky).gradX = (i==0 || i==row.size()-1) ? 0.0 : 1.0;
                    g.get(kx+idx, ky).gradY = (i==0 ? -1.0 : (i==row.size()-1 ? 1.0 : 0.0));
                }
            }
            else {
                // mark the two surounding points
                bound.insert(g.pointsMap_[std::make_pair(kx, ky)]);
                g.get(kx, ky).gradX = (i==0 || i==row.size()-1) ? 0.0 : -1.0;
                g.get(kx, ky).gradY = (i==0 ? -1.0 : (i==row.size()-1 ? 1.0 : 0.0));
                bound.insert(g.pointsMap_[std::make_pair(kx+1, ky)]);
                g.get(kx, ky).gradX = (i==0 || i==row.size()-1) ? 0.0 : 1.0;
                g.get(kx, ky).gradY = (i==0 ? -1.0 : (i==row.size()-1 ? 1.0 : 0.0));
                for(int idx=0; idx<nbInternalBoundaryExtraLayers; idx++) {
                    bound.insert(g.pointsMap_[std::make_pair(kx-idx, ky)]);
                    g.get(kx-idx, ky).gradX = (i==0 || i==row.size()-1) ? 0.0 : -1.0;
                    g.get(kx-idx, ky).gradY = (i==0 ? -1.0 : (i==row.size()-1 ? 1.0 : 0.0));
                    bound.insert(g.pointsMap_[std::make_pair(kx+idx, ky)]);
                    g.get(kx+idx, ky).gradX = (i==0 || i==row.size()-1) ? 0.0 : 1.0;
                    g.get(kx+idx, ky).gradY = (i==0 ? -1.0 : (i==row.size()-1 ? 1.0 : 0.0));
                }
            }
            if (abs(y - (double) ky * dy) < epsi) {
                // mark only one point
                bound.insert(g.pointsMap_[std::make_pair(kx, ky)]);
                g.get(kx, ky).gradX = (i==0 || i==row.size()-1) ? 0.0 : -1.0;
                g.get(kx, ky).gradY = (i==0 ? -1.0 : (i==row.size()-1 ? 1.0 : 0.0));

                for(int idx=0; idx<nbInternalBoundaryExtraLayers; idx++) {
                    bound.insert(g.pointsMap_[std::make_pair(kx-idx, ky)]);
                    g.get(kx-idx, ky).gradX = (i==0 || i==row.size()-1) ? 0.0 : -1.0;
                    g.get(kx-idx, ky).gradY = (i==0 ? -1.0 : (i==row.size()-1 ? 1.0 : 0.0));
                    bound.insert(g.pointsMap_[std::make_pair(kx+idx, ky)]);
                    g.get(kx+idx, ky).gradX = (i==0 || i==row.size()-1) ? 0.0 : 1.0;
                    g.get(kx+idx, ky).gradY = (i==0 ? -1.0 : (i==row.size()-1 ? 1.0 : 0.0));
                }
            }
            else {
                // mark the two surounding points
                bound.insert(g.pointsMap_[std::make_pair(kx, ky)]);
                g.get(kx, ky).gradX = (i==0 || i==row.size()-1) ? 0.0 : -1.0;
                g.get(kx, ky).gradY = (i==0 ? -1.0 : (i==row.size()-1 ? 1.0 : 0.0));
                bound.insert(g.pointsMap_[std::make_pair(kx, ky+1)]);
                g.get(kx, ky+1).gradX = (i==0 || i==row.size()-1) ? 0.0 : -1.0;
                g.get(kx, ky+1).gradY = (i==0 ? -1.0 : (i==row.size()-1 ? 1.0 : 0.0));
                for(int idx=0; idx<nbInternalBoundaryExtraLayers; idx++) {
                    bound.insert(g.pointsMap_[std::make_pair(kx-idx, ky)]);
                    g.get(kx-idx, ky).gradX = (i==0 || i==row.size()-1) ? 0.0 : -1.0;
                    g.get(kx-idx, ky).gradY = (i==0 ? -1.0 : (i==row.size()-1 ? 1.0 : 0.0));
                    bound.insert(g.pointsMap_[std::make_pair(kx+idx, ky)]);
                    g.get(kx+idx, ky).gradX = (i==0 || i==row.size()-1) ? 0.0 : 1.0;
                    g.get(kx+idx, ky).gradY = (i==0 ? -1.0 : (i==row.size()-1 ? 1.0 : 0.0));
                }
            }
        }
    }

    boundary.insert(boundary.end(), bound.begin(), bound.end());

    reduceDoorsSize(doors, wall, exits, g, boundary, outFluxLimiter);
    // just to fix a bug of unknown origin
    for(int j=0; j<=g.jMax_; j++) {
        if(std::find(exits.begin(), exits.end(), g.pointsMap_[std::make_pair(g.iMax_, j)]) == exits.end()
       &&  std::find(boundary.begin(), boundary.end(), g.pointsMap_[std::make_pair(g.iMax_, j)]) == boundary.end()) {
            boundary.push_back(g.pointsMap_[std::make_pair(g.iMax_, j)]);
            g.get(g.iMax_, j).gradX = -1.0;
            g.get(g.iMax_, j).gradY = (j==0  ? 1.0 : (j == g.jMax_ ? -1.0 : 0.0));
        }
    }

    grid gg(g.points_, boundary, exits, dx, dy, doors);

    FMM fmm(gg);

    fmm.run();

    std::string test("cinema.dat");
    std::string vtest("vcinema.dat");
    std::string btest("bcinema.dat");
    std::string etest("ecinema.dat");
    fmm.output(test, vtest, btest, etest);
    std::cout << " ******* FMM completed *****" << std::endl;

    return gg;
}

grid Zenith::initializeGrid(Wall &wall, double dx, double dy){
    double xmin{-50}, ymin{-20}, xmax{50}, ymax{60};
    std::vector<gridPoint> gridPoints;
    initialiseGrid(gridPoints, xmin, xmax, ymin, ymax, dx, dy);
    std::vector<int> boundary, exits;
    std::vector<door> doors;

    grid g(gridPoints, boundary, exits, dx, dy, doors, xmin, ymin);

    //defining boundaries, exits and doors
    for(int i=0; i<=g.iMax_; i++){
        boundary.push_back(g.pointsMap_[std::make_pair(i, 0)]);
        boundary.push_back(g.pointsMap_[std::make_pair(i, g.jMax_)]);
    }
    for(int j=0; j<=g.jMax_; j++){
        boundary.push_back(g.pointsMap_[std::make_pair(0, j)]);
        boundary.push_back(g.pointsMap_[std::make_pair(g.iMax_, j)]);
    }

    std::string doorsFname {"../inputs/zenith_up_doors.txt"};
    std::vector<double> tmp {getVectorFromFile(doorsFname)};
    std::vector<int> nbPointsPerDoorRows (3, g.iMax_);
    for(int i=0; i<tmp.size(); i+=8) {
        std::vector<double> gateX(4), gateY(4);
        for(int j=0; j<4; j++){
            gateX[j] = tmp[i + j*2];
            gateY[j] = tmp[i + j*2 + 1];
        }
        polygonalWall porte(gateX, gateY, nbPointsPerDoorRows);
        placeAPolygonOnBoundary(porte, g, boundary, true);
        if(i < 6 * 8) // Portes du haut 
        {   
            gridPoint a(gateX[0], gateY[0], 0.0);
            gridPoint b(gateX[3], gateY[3], 0.0);
            door d (a,b);
            doors.push_back(d);
            addADoorAsExit(d, g, exits);
        }
        else { // portes milieu
            gridPoint a(gateX[0], gateY[0], 0.0);
            gridPoint b(gateX[3], gateY[3], 0.0);
            door d (a,b);
            doors.push_back(d);
            addADoorAsExit(d, g, exits);
        }
    }

    tmp.clear();
    std::string wallSupFname {"../inputs/zenith_ext.txt"};
    tmp = getVectorFromFile(wallSupFname);
    std::vector<double> verticesWallUpX(15), verticesWallUpY(15);
    int k=0;
    for(int i=0; i < 15*2; i+=2) {
        verticesWallUpX[k] = tmp[i];
        verticesWallUpY[k++] = tmp[i+1];
    }
    std::vector<int> wallUpNbPoints(14, g.iMax_);
    polygonalWall wallSup(verticesWallUpX, verticesWallUpY, wallUpNbPoints);
    placeAPolygonOnBoundary(wallSup, g, boundary, false);

    k = 0;
    std::vector<double> verticesWallDownX(8), verticesWallDownY(8);
    for(int i=15*2; i<23*2; i+=2) {
        verticesWallDownX[k] = tmp[i];
        verticesWallDownY[k++] = tmp[i+1];
    }
    std::vector<int> wallDownNbPoints(7, g.iMax_);
    polygonalWall wallDown(verticesWallDownX, verticesWallDownY, wallDownNbPoints);
    placeAPolygonOnBoundary(wallDown, g, boundary, false);

    for(int i=23*2; i<tmp.size(); i+=4) {
        std::vector<double> X, Y;
        X.push_back(tmp[i]);
        Y.push_back(tmp[i+1]);
        X.push_back(tmp[i+2]);
        Y.push_back(tmp[i+3]);
        std::vector<int> nbPoints{g.iMax_};
        polygonalWall internalWall(X, Y, nbPoints);
        placeAPolygonOnBoundary(internalWall, g, boundary, true);
    }
    gridPoint a(tmp[23*2], tmp[23*2 + 1], 0.0);
    gridPoint b(tmp[22*2], tmp[22*2 + 1], 0.0);
    door d (a,b);
    doors.push_back(d);
    addADoorAsExit(d, g, exits);
    a.x *= -1;
    b.x *= -1;
    door dd (a,b);
    doors.push_back(dd);
    addADoorAsExit(dd, g, exits);

    tmp.clear();
    std::string frontWallFname {"../inputs/zenith_front_walls.txt"};
    tmp = getVectorFromFile(frontWallFname);
    for(int wallIdx = 0; wallIdx<2; wallIdx++) {
        std::vector<double> X(4), Y(4);
        std::vector<int> nbPoints {g.iMax_, g.iMax_, g.iMax_};
        int k=0;
        for(int i=wallIdx*8; i < (wallIdx+1) * 8; i+= 2) {
            X[k] = tmp[i];
            Y[k++] = tmp[i+1];
        }
        polygonalWall fw(X,Y, nbPoints);
        placeAPolygonOnBoundary(fw, g, boundary, true);
    }
    {
        std::vector<double> X(5), Y(5);
        std::vector<int> nbPoints(4, g.iMax_);
        int k=0;
        for(int i=8*2; i<tmp.size(); i+=2) {
            X[k] = tmp[i];
            Y[k++] = tmp[i+1];
        }
        X.back() = X[0];
        Y.back() = Y[0];
        polygonalWall fw(X,Y, nbPoints);
        placeAPolygonOnBoundary(fw, g, boundary, true);
    }

    reduceDoorsSize(doors, wall, exits, g, boundary, outFluxLimiter);
    grid gg(g.points_, boundary, exits, dx, dy, doors, xmin, ymin);

    FMM fmm(gg);

    fmm.run();

    std::string test("zenith.dat");
    std::string vtest("vzenith.dat");
    std::string btest("bzenith.dat");
    std::string etest("ezenith.dat");
    fmm.output(test, vtest, btest, etest);
    std::cout << " ******* FMM completed *****" << std::endl;
    return gg;
}
