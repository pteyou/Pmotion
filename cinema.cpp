#include <vector>
#include <string>
#include <fstream>
#include <iterator>
#include "polygonalWall.h"
#include <iostream>
#include <chrono>
#include <algorithm>
#include <random>
#include "myTypes.h"
#include "eulerScheme.h"
#include "totalForces.h"
#include "FMM/FMM.h"

#define distanceToBorder 1
#define centralRowWidth 2
#define CAPKMODE 0

int closerGateIdx (double x0, double y0, std::vector<double> &gatesX, std::vector<double> &gatesY) {
    double minDist = L2norm(x0 - gatesX[0], y0 - gatesY[0]);
    int res = 0;
    for(int i=1; i < gatesX.size(); i++) {
        double dist = L2norm(x0 - gatesX[i], y0 - gatesY[i]);
        if (minDist > dist) {
            minDist = dist;
            res = i;
        }
    }
    return res;
}

inline bool movesDown(double y, double gateY, double upY, double downY) {
    if (abs(downY - gateY) < abs (upY - gateY)) return true;
    else return false;
}

void setArgumentsSFM(Arguments &args, double dt, double lambda, double Ai, double Bi, 
    double ki, double kappai, double maximalRepulsionForce, double k, double D0, double D1,
    double relaxationTime){
    int nbIndividual = args.nbIndividual;
    args.dt = dt;
    args.lambda = lambda;
    args.Ai = Ai;
    args.Bi = Bi;
    args.ki = ki;
    args.kappai = kappai;
    args.maximalRepulsionForce = maximalRepulsionForce;
    args.D0 = D0;
    args.D1 = D1;
    args.k = k;
    for(int i=0; i<nbIndividual; i++) args.relaxationTime->coeffRef(i) = relaxationTime;
}

EulerBase* EulerBase::me = nullptr;

int main(int argc, char *argv[]) {
    //omp_set_num_threads(1);
    auto start = std::chrono::high_resolution_clock::now();

    if (CAPKMODE) assert(argc == 12);

    std::vector<double> upWallX{27.9, 27.9};
    std::vector<double> upWallY{5.5, 10.5};
    std::vector<int> nbPointsUp{30};
    polygonalWall up(upWallX, upWallY, nbPointsUp);

    std::vector<double> leftWallX{21, 0.0, 0.0, 2.7}; // X3 is measured 2.82
    std::vector<double> leftWallY{16.0, 16.0, 0.0, 0.0};
    std::vector<int> nbPointsLeftWall{40, 30, 5};
    polygonalWall left(leftWallX, leftWallY, nbPointsLeftWall);

    std::vector<double> rightWallX{4.7, 21}; // approx 23.4 - 2
    std::vector<double> rightWallY{0.0, 0.0};
    std::vector<int> nbPointsRightWall{30};
    polygonalWall right(rightWallX, rightWallY, nbPointsRightWall);

    std::vector<double> rightUpX{17.4, 23.4, 23.4, 27.9, 27.9};
    std::vector<double> rightUpY{1.2, 1.2, 0.0, 0.0, 4.0 };
    std::vector<int> nbPointsRightUp{70, 10, 10, 10};
    polygonalWall rightUp(rightUpX, rightUpY, nbPointsRightUp);

    std::vector<double> leftUpX{17.4, 23.4, 23.4, 27.9, 27.9};
    std::vector<double> leftUpY{14.8, 14.8, 16, 16, 12};
    std::vector<int> nbPointsLeftUp{30, 10, 10, 10};
    polygonalWall leftUp(leftUpX, leftUpY, nbPointsLeftUp);

    std::vector<polygonalWall> wallVector {left, right, leftUp, rightUp, up};
    
    
    Arguments args {
        1.3,                                    //maxSpeedFactor
        350, 3000,                             // nbInidividual, nbTimeSteps -> nbIndividual must be multiple of 14
        0, 0,                               // iIndividual, jIndividual
        0,0,                                  // iwall, wallSize
        0,                                  // currentTimeStep
        2.0 / 100.0,                                // dt
        0.25,                               // lambda
        160,                                // F : maximalRepulsionForce 
        2,0.31,0.45,                              //k, D0, D1
        //2000,0.08,1.2e5,2.5e5,                            // Ai, Bi, ki, kappai
        2000,0.08,5.0e3,2.5e5,                            // Ai, Bi, ki, kappai
        nullptr, nullptr, nullptr, nullptr,                  //velocityX, velocityY,positionX, positionY
        nullptr, nullptr,                           //desiredSpeed, maxSpeed;
        nullptr,                            //relaxationTime;
        nullptr, nullptr,                   //normalX, normalY;
        nullptr, nullptr,                   //rijX, rijY;
        nullptr,                            //dij;
        nullptr, nullptr,                   //riWX, riWY;
        nullptr,                            //diW;
        nullptr,                            //Radius;
        nullptr,                             //mass;
        5.0 * M_PI / 12.0,                          // sight angle range
        //M_PI / 6,
        100,                                 // moving angle resolution
        10,                                  // sight distance range
        4                                    // maximum number of intermediate steps in the desired path
    };

    int nbRows = 2 * 7;
    int nbPersonPerRow = static_cast<int> (args.nbIndividual / nbRows);
    assert(nbPersonPerRow * nbRows == args.nbIndividual);
    DVector Y = DVector::LinSpaced(nbPersonPerRow, 1.2 + distanceToBorder, 14.9 - distanceToBorder);

    // positions relaxation in regard to the central row
    double downY = 1.2 + distanceToBorder;
    double upY = 14.9 - distanceToBorder;
    double centralDownY = 0.5*(upY + downY - centralRowWidth);
    double centralUpY = centralDownY + centralRowWidth;
    double middleY = 0.5 * (upY + downY);
    double scaleDown = (centralDownY - downY) / (middleY - downY);
    for (int i=0; i<nbPersonPerRow/2; i++) {
        Y[i] = downY + scaleDown * (Y[i] - downY);
    }
    double scaleUp = (centralUpY - upY) / (middleY - upY);
    for(int i=nbPersonPerRow/2; i<nbPersonPerRow; i++) {
        Y[i] = upY + scaleUp * (Y[i] - upY);
    }

    DVector X1 = DVector::LinSpaced(7, 17+distanceToBorder, 25.8-distanceToBorder);

    DVector initPositionX(args.nbIndividual), initPositionY(args.nbIndividual);
    for(int i=0; i<7; i++) {
        for(int j=0; j<nbPersonPerRow; j++) {
            initPositionX[i*nbPersonPerRow + j] = X1[i];
        }
    }
    for(int i=0; i<7; i++) {
        initPositionY.segment(i*nbPersonPerRow, nbPersonPerRow) = Y;
    }

    initPositionY.segment(args.nbIndividual / 2, args.nbIndividual/2) = initPositionY.segment(0, args.nbIndividual/2);
    for(int i=0; i<args.nbIndividual/2; i++) {
        initPositionX[i + args.nbIndividual/2] = initPositionX[i] - 10.5;
    }

    // creation of walls modelling the chairs rows
    std::vector<polygonalWall> rowWalls;
    rowWalls.reserve(nbRows * 2 - 2);
    std::vector<int> nbPointsPerRow {30};
    double dx = 0.5*(initPositionX[nbPersonPerRow] + initPositionX[0]) - initPositionX[0];
    for(int i = 0; i < nbRows; i++){
        double xRightDown = initPositionX[i*nbPersonPerRow] + dx;
        double XRightUp = xRightDown;
        double yRightDown = initPositionY[i*nbPersonPerRow];
        double yRightUp = initPositionY[i*nbPersonPerRow + nbPersonPerRow/2 - 1];
        std::vector<double> X{xRightDown, XRightUp};
        std::vector<double> Y{yRightDown, yRightUp};
        polygonalWall tmp{X, Y, nbPointsPerRow};
        rowWalls.push_back(tmp);

        double xLeftDown = xRightDown;
        double XLeftUp = xRightDown;
        double yLeftDown = initPositionY[i*nbPersonPerRow + nbPersonPerRow/2];
        double yLeftUp = initPositionY[(i+1)*nbPersonPerRow - 1];
        std::vector<double> Xl{xLeftDown, XLeftUp};
        std::vector<double> Yl{yLeftDown, yLeftUp};
        polygonalWall tmpl{Xl, Yl, nbPointsPerRow};
        rowWalls.push_back(tmpl);
    }
    wallVector.insert(wallVector.end(), rowWalls.begin(), rowWalls.end());

    /*
    std::string rwalls("cinemaRowWalls.dat");
    std::ofstream wallf(rwalls.c_str());  
    for(auto w : rowWalls) {
        for(int i=0; i < w.size(); i++){
            std::vector<double> point{w.point(i)};
            wallf << point[0]  << "\t" << point[1] << std::endl;
        }
    }
    wallf.close(); */

    Wall w(wallVector);

    std::vector<double> gatesX{27.9, 3.7, 22.3, 27.9, 22.3}, gatesY{4.75, 0.0, 0.0, 11.25, 16};
    std::vector<double> intermediateX{27.9, 3.7, 16.0, 27.9,16}, intermediateY{4.75, 0.0, 0.6, 11.25, 15.4};
    std::vector<double> finalX{35.0, 3.7, 22.3, 35, 22.3}, finalY{4.75, -10, -10, 11.25, 36};

    int nbIndividual = args.nbIndividual;

    DVector initVelocityX(nbIndividual), initVelocityY(nbIndividual);
    for(int i=0; i<nbIndividual; i++) initVelocityX[i] = 0.0;
    for(int i=0; i<nbIndividual; i++) initVelocityY[i] = 0.0;

    DMatrix finalPositionX(args.maxStepsInPath, nbIndividual), finalPositionY(args.maxStepsInPath, nbIndividual);

    for(int i=0; i<nbIndividual; i++){
        int ipath = closerGateIdx(initPositionX[i],initPositionY[i], gatesX, gatesY);
        double UPY = initPositionY[i] <= centralDownY ? centralDownY : upY;
        double DOWNY = initPositionY[i] <= centralDownY ? downY : centralUpY;
        finalPositionX(0,i) = initPositionX[i];
        finalPositionX(1,i) = intermediateX[ipath];
        finalPositionX(2,i) = gatesX[ipath];
        finalPositionX(3,i) = finalX[ipath];
        // ensures that every pedestrian lives the entire sitting row
        finalPositionY(0,i) = movesDown(initPositionY[i], gatesY[ipath], UPY, DOWNY) ? DOWNY - DistanceEps : UPY + DistanceEps;
        finalPositionY(1,i) = intermediateY[ipath];
        finalPositionY(2,i) = gatesY[ipath];
        finalPositionY(3,i) = finalY[ipath];
    }
    
    DVector Radius(nbIndividual);
    for(int i=0; i<nbIndividual; i++) Radius[i] = 0.25;
    args.Radius = &Radius;

    DVector relaxationTime(nbIndividual);
    args.relaxationTime = &relaxationTime;

    if(CAPKMODE) setArgumentsSFM(args, std::stod(argv[1]), std::stod(argv[2]), std::stod(argv[3]), 
        std::stod(argv[4]), std::stod(argv[5]), std::stod(argv[6]), std::stod(argv[7]), std::stod(argv[8]),
        std::stod(argv[9]), std::stod(argv[10]), std::stod(argv[11]));
    else
    {
        setArgumentsSFM(args, args.dt, args.lambda, args.Ai, args.Bi, 
            args.ki, args.kappai, args.maximalRepulsionForce, args.k, args.D0, args.D1, 0.25);
    }

    double meanVelocity { 1.3 };
    double stddeviation { 0.15 };
    DVector desiredSpeed(nbIndividual), maxSpeed(nbIndividual);
    std::normal_distribution<double> distributionSpeed(meanVelocity, stddeviation);
    std::default_random_engine generatorSpeed;
    for(int i=0; i < nbIndividual; i++) desiredSpeed[i] = distributionSpeed(generatorSpeed);
    args.desiredSpeed = &desiredSpeed;
    maxSpeed = args.maxSpeedFactor * desiredSpeed;
    args.maxSpeed = &maxSpeed;

    std::vector<std::function<std::vector<double>(Arguments &, Wall &)>> forces;
    forces.push_back(desiredSpeedMotivationSFM);
    forces.push_back(wallRepulsionSFM);
    forces.push_back(individualRepulsionCrowdTurbulence);
    totalForces TF(forces);

    DVector mass(nbIndividual);
    for(int i=0; i<nbIndividual; i++) mass[i] = 80.0;
    
    args.mass = &mass;
    args.wallSize = w.size();

    std::cout << "*** starting simulation for the movie theater case ***" << std::endl;
    std::cout << "*********  with  " << args.nbIndividual << "  points ************" << std::endl;

    //grid g = initializeGridCinema();

    /* Euler iter(TF, w, args); */
    double dxx = 0.15, dyy=0.1;
    EulerFMM<Cinema> iter(TF, w, args, dxx, dyy);
    int nbTseps = iter.run(initPositionX, initPositionY, initVelocityX, initVelocityY, finalPositionX, 
        finalPositionY);
    iter.output("movie350X.dat", "movie350Y.dat", "movie350VX.dat", "movie350VY.dat", "movie350Wall.dat");

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

    std::cout << "total elapsed time " << duration.count() << " seconds" << std::endl;
    std::cout << "Physical time (seconds) and finished (0 or 1) " << nbTseps * args.dt << " " << 
        (nbTseps == args.nbTimeSteps ? 0 : 1) << std::endl;
}