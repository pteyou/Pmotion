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
#include "placePeopleInZenith.h"
#include <omp.h>
#include <csignal>

#define maximalInitialCompression 0.1

bool regularizePositions(DVector &X, DVector &Y, DVector &Radius) {
    bool flag {false};
    DMatrix rijX, rijY, dij;
    std::vector<int> idxToErase;
    do {
        long nbPoints {X.size()};
        rijX.resize(nbPoints, nbPoints);
        rijY.resize(nbPoints, nbPoints);
        dij.resize(nbPoints, nbPoints);
        if (!idxToErase.empty()) idxToErase.clear();
        idxToErase.reserve(200);
        for(int i=0; i<nbPoints; i++) {
            for(int j=0; j<i; j++) {
                rijX(i,j) = X[i] - X[j];
                rijY(i,j) = Y[i] - Y[j];
                dij(i,j) = L2norm(rijX(i,j), rijY(i,j)) - maximalInitialCompression;
                // dij(i,j) = L2norm(rijX(i,j), rijY(i,j)) - Radius(i) - Radius(j);
                if (dij(i,j) <= 0.0) { idxToErase.push_back(i); break; }
            }
        }
        if(idxToErase.empty()) flag = true;
        else {
            std::sort(idxToErase.begin(), idxToErase.end(), std::greater<int>());
            idxToErase.erase( std::unique( idxToErase.begin(), idxToErase.end() ), idxToErase.end() );
            for(auto i : idxToErase) {
            //erase the rightmost points
                std::cout << "$$$ erasing point number " << i << std::endl;
                std::copy(X.data() + i + 1, X.data() + X.size(), X.data() + i);
                X.conservativeResize(X.size() -1);
                std::copy(Y.data() + i + 1, Y.data() + Y.size(), Y.data() + i);
                Y.conservativeResize(Y.size() -1);
                std::copy(Radius.data() + i + 1, Radius.data() + Radius.size(), Radius.data() + i);
                Radius.conservativeResize(Radius.size() -1);
            }
            flag = true;
        }
    }
    while (!flag);
    return flag;
}

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

EulerBase* EulerBase::me = nullptr;

int main() {
    //omp_set_num_threads(1);
    auto start = std::chrono::high_resolution_clock::now();
    
    // wall and gates definition
    std::vector<double> gatesX, gatesY, intermediateX, intermediateY;
    gatesX.reserve(12); gatesY.reserve(12); 
    intermediateX.reserve(12); intermediateY.reserve(12);
    std::vector<polygonalWall> wallVector;
    wallVector.reserve(20);
    std::vector<double> tmp {getVectorFromFile("../inputs/zenith_up_doors.txt")};
    std::vector<int> nbPointsPerDoorRows (3, 15);
    for(int i=0; i<tmp.size(); i+=8) {
        std::vector<double> gateX(4), gateY(4);
        for(int j=0; j<4; j++){
            gateX[j] = tmp[i + j*2];
            gateY[j] = tmp[i + j*2 + 1];
        }
        wallVector.push_back(polygonalWall(gateX, gateY, nbPointsPerDoorRows));
        if(i < 6 * 8) // Portes du haut 
        {   
            gatesX.push_back(0.5 * (gateX[0] + gateX[3]));
            gatesY.push_back(0.5 * (gateY[0] + gateY[3]));
            intermediateX.push_back(gatesX.back());
            intermediateY.push_back(gatesY.back());
        }
        else { // portes milieu
            gatesX.push_back(0.5 * (gateX[1] + gateX[2]));
            gatesY.push_back(0.5 * (gateY[1] + gateY[2]));
            intermediateX.push_back(0.5 * (gateX[0] + gateX[3]));
            intermediateY.push_back(0.5 * (gateY[0] + gateY[3]));
        }
        
    }

    tmp.clear();
    tmp = getVectorFromFile("../inputs/zenith_ext.txt");
    std::vector<double> verticesWallUpX(15), verticesWallUpY(15);
    int k=0;
    for(int i=0; i < 15*2; i+=2) {
        verticesWallUpX[k] = tmp[i];
        verticesWallUpY[k++] = tmp[i+1];
    }
    std::vector<int> wallUpNbPoints(14, 15);
    wallUpNbPoints[0] = 30;
    wallUpNbPoints.back() = 30;
    wallVector.push_back(polygonalWall(verticesWallUpX, verticesWallUpY, 
        wallUpNbPoints));

    k = 0;
    std::vector<double> verticesWallDownX(8), verticesWallDownY(8);
    for(int i=15*2; i<23*2; i+=2) {
        verticesWallDownX[k] = tmp[i];
        verticesWallDownY[k++] = tmp[i+1];
    }
    std::vector<int> wallDownNbPoints(7, 15);
    wallVector.push_back(polygonalWall(verticesWallDownX, verticesWallDownY, 
        wallDownNbPoints));

    for(int i=23*2; i<tmp.size(); i+=4) {
        std::vector<double> X, Y;
        X.push_back(tmp[i]);
        Y.push_back(tmp[i+1]);
        X.push_back(tmp[i+2]);
        Y.push_back(tmp[i+3]);
        std::vector<int> nbPoints{15};
        wallVector.push_back(polygonalWall(X, Y, nbPoints));
    }
    gatesX.push_back(0.5*(tmp[23*2] + tmp[22*2]));
    gatesY.push_back(0.5*(tmp[23*2 + 1] + tmp[22*2 + 1]));
    intermediateX.push_back(gatesX.back());
    intermediateY.push_back(gatesY.back());
    gatesX.push_back(-gatesX.back());
    gatesY.push_back(gatesY.back());
    intermediateX.push_back(gatesX.back());
    intermediateY.push_back(gatesY.back());

    // front walls
    tmp.clear();
    tmp = getVectorFromFile("../inputs/zenith_front_walls.txt");
    for(int wallIdx = 0; wallIdx<2; wallIdx++) {
        std::vector<double> X(4), Y(4);
        std::vector<int> nbPoints {20, 15, 10};
        int k=0;
        for(int i=wallIdx*8; i < (wallIdx+1) * 8; i+= 2) {
            X[k] = tmp[i];
            Y[k++] = tmp[i+1];
        }
        wallVector.push_back(polygonalWall(X,Y, nbPoints));
    }
    {
        std::vector<double> X(5), Y(5);
        std::vector<int> nbPoints(4, 100);
        int k=0;
        for(int i=8*2; i<tmp.size(); i+=2) {
            X[k] = tmp[i];
            Y[k++] = tmp[i+1];
        }
        X.back() = X[0];
        Y.back() = Y[0];
        wallVector.push_back(polygonalWall(X, Y, nbPoints));
    }

    Wall w(wallVector);

    Arguments args {
        1.3,                                    //maxSpeedFactor
        0, 6000,                             // nbInidividual, nbTimeSteps 
        0, 0,                               // iIndividual, jIndividual
        0,0,                                  // iwall, wallSize
        0,                                  // currentTimeStep
        1.0 / 100.0,                                // dt
        0.5,                               // lambda
        160,                                //maximalRepulsionForce
        2,0.31,0.45,                              //k, D0, D1
        2000,0.08,1.2e5,2.5e5,                            // Ai, Bi, ki, kappai
        //2000,0.08,5.0e3,2.5e5,                            // Ai, Bi, ki, kappai
        nullptr, nullptr, nullptr, nullptr,                  //velocityX, velocityY,positionX, positionY
        nullptr, nullptr,                           //desiredSpeed, maxSpeed;
        nullptr,                            //relaxationTime;
        nullptr, nullptr,                   //normalX, normalY;
        nullptr, nullptr,                   //rijX, rijY;
        nullptr,                            //dij;
        nullptr, nullptr,                   //riWX, riWY
        nullptr,                            //diW;
        nullptr,                            //Radius
        nullptr,                             //mass
        //5.0 * M_PI / 12.0,                          // sight angle range
        M_PI / 6,
        100,                                 // moving angle resolution
        10,                                  // sight distance range
        2                                    // maximum number of intermediate steps in the desired path
    };

    std::vector<double> initX, initY;
    placePeople(initX, initY);
    DVector initPositionX {Map<DVector> (initX.data(), initX.size())};
    DVector initPositionY {Map<DVector> (initY.data(), initY.size())};

    int nbIndividual = (args.nbIndividual = initPositionX.size());

    DVector initVelocityX(nbIndividual), initVelocityY(nbIndividual);
    for(int i=0; i<nbIndividual; i++) initVelocityX[i] = 0.0;
    for(int i=0; i<nbIndividual; i++) initVelocityY[i] = 0.0;


    DMatrix finalPositionX(args.maxStepsInPath, nbIndividual);
    DMatrix finalPositionY(args.maxStepsInPath, nbIndividual);

    for(int i=0; i<nbIndividual; i++){
        int ipath = closerGateIdx(initPositionX[i],initPositionY[i], 
            gatesX, gatesY);
        finalPositionX(0,i) = intermediateX[ipath];
        finalPositionX(1,i) = gatesX[ipath];
        finalPositionY(0,i) = intermediateY[ipath];
        finalPositionY(1,i) = gatesY[ipath];
    }

    DVector Radius(nbIndividual);
    for(int i=0; i<nbIndividual; i++) Radius[i] = 0.25;
    args.Radius = &Radius;

    DVector relaxationTime(nbIndividual);
    for(int i=0; i<nbIndividual; i++) relaxationTime[i] = 0.5;
    args.relaxationTime = &relaxationTime;

    double meanVelocity { 5.0 };
    double stddeviation { 0.2 };
    DVector desiredSpeed(nbIndividual), maxSpeed(nbIndividual);
    std::normal_distribution<double> distributionSpeed(meanVelocity, stddeviation);
    std::default_random_engine generatorSpeed;
    for(int i=0; i < nbIndividual; i++) 
        desiredSpeed[i] = distributionSpeed(generatorSpeed);
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

    std::cout << "*** starting simulation for the zenith case ***" << std::endl;
    std::cout << "*********  with  " << args.nbIndividual << 
        "  points ************" << std::endl;

    //Euler iter(TF, w, args);

    signal(SIGINT, Euler::signalHandler);

    /* iter.run(initPositionX, initPositionY, initVelocityX, initVelocityY, finalPositionX, 
        finalPositionY); */
    double dxx = 0.2, dyy=0.2;
    EulerFMM<Zenith> iter(TF, w, args, dxx, dyy);
    int nbTseps = iter.run(initPositionX, initPositionY, initVelocityX, 
        initVelocityY, finalPositionX, finalPositionY);
    iter.output("zenithXFMM.dat", "zenithYFMM.dat", "zenithVXFMM.dat", 
        "zenithVYFMM.dat", "zenithWall.dat");

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

    std::cout << "total elapsed time " << duration.count() << " seconds" << std::endl;
    std::cout << "Physical time (seconds) " << nbTseps * args.dt << std::endl;

}
