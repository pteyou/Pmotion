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
#include "utilities.h"

#define distanceFromTable 0.5
#define CAPKMODE 0

void defineInitialPositionsMeetingRoom(DVector &initX, DVector &initY, polygonalWall &table) {
    double Xcenter = std::accumulate(table.pointsX.begin(), table.pointsX.end(), 0.0) / table.size();
    double Ycenter = std::accumulate(table.pointsY.begin(), table.pointsY.end(), 0.0) / table.size();
    int idx=0;
    for (int pos = 0; pos < initX.size() * static_cast<int> (table.size() / initX.size())-1; 
        pos += static_cast<int> (table.size() / initX.size())){
        double dx = table.pointsX[pos] - Xcenter;
        double dy = table.pointsY[pos] - Ycenter;
        double dist = sqrt(dx*dx + dy*dy);
        initX[idx] = Xcenter + (dist + distanceFromTable) * dx / dist;
        initY[idx++] = Ycenter + (dist + distanceFromTable) * dy / dist;
    }
}

void setArgumentsSFM(Arguments &args, double dt, double lambda, double Ai, double Bi, 
    double ki, double kappai, double relaxationTime){
    int nbIndividual = args.nbIndividual;
    args.dt = dt;
    args.lambda = lambda;
    args.Ai = Ai;
    args.Bi = Bi;
    args.ki = ki;
    args.kappai = kappai;
    for(int i=0; i<nbIndividual; i++) args.relaxationTime->coeffRef(i) = relaxationTime;
}

EulerBase* EulerBase::me = nullptr;

int main(int argc, char* argv[]) {
    //omp_set_num_threads(1);
    
    auto start = std::chrono::high_resolution_clock::now();

    if (CAPKMODE) assert(argc == 8);
    
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
    
    Arguments args {
        1.3,                                    //maxSpeedFactor
        10, 1000,                             // nbInidividual, nbTimeSteps
        0, 0,                               // iIndividual, jIndividual
        0,0,                                  // iwall, wallSize
        0,                                  // currentTimeStep
        2.0 / 100.0,                                // dt
        0.25,                               // lambda
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
        nullptr, nullptr,                   //riWX, riWY;
        nullptr,                            //diW;
        nullptr,                            //Radius;
        nullptr,                             //mass;
        5.0 * M_PI / 12.0,                          // sight angle range
        //M_PI / 6.0,
        100,                                 // moving angle resolution
        10,                                  // sight distance range
        2                                    // maximum number of intermediate steps in the desired path
    };

    double dx = 0.15, dy = 0.1;

    int nbIndividual = args.nbIndividual;
    DVector relaxationTime(nbIndividual);
    args.relaxationTime = &relaxationTime;
    
    if(CAPKMODE) setArgumentsSFM(args, std::stod(argv[1]), std::stod(argv[2]), std::stod(argv[3]), 
        std::stod(argv[4]), std::stod(argv[5]), std::stod(argv[6]), std::stod(argv[7]));
    else
    {
        setArgumentsSFM(args, args.dt, args.lambda, args.Ai, args.Bi, 
            args.ki, args.kappai, 0.25);
    }
    

    DVector initPositionX(nbIndividual), initPositionY(nbIndividual);
    defineInitialPositionsMeetingRoom(initPositionX, initPositionY, table);
    

    DMatrix finalPositionX(args.maxStepsInPath, nbIndividual), finalPositionY(args.maxStepsInPath, nbIndividual);
    DVector initVelocityX(nbIndividual), initVelocityY(nbIndividual);

    for(int i=0; i<nbIndividual; i++) initVelocityX[i] = 0.0;
    for(int i=0; i<nbIndividual; i++) initVelocityY[i] = 0.0;
    
    for(int i=0; i<nbIndividual; i++) finalPositionX(0,i) = 0.0;
    for(int i=0; i<nbIndividual; i++) finalPositionX(1,i) = -0.5;
    for(int i=0; i<nbIndividual; i++) finalPositionY(0,i) = 0.96;
    for(int i=0; i<nbIndividual; i++) finalPositionY(1,i) = 0.96;

    DVector Radius(nbIndividual);
    for(int i=0; i<nbIndividual; i++) Radius[i] = 0.25;
    args.Radius = &Radius;

    

    double meanVelocity { 1.3 };
    double stddeviation { 0.2 };
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

    std::cout << "*** starting simulation for the meeting room case ***" << std::endl;
    std::cout << "*********  with  " << args.nbIndividual << "  points ************" << std::endl;

    EulerFMM<Reunion> iter(TF, w, args, dx, dy);
    int nbTseps = iter.run(initPositionX, initPositionY, initVelocityX, initVelocityY, finalPositionX, 
        finalPositionY);
    iter.output("room12X.dat", "room12Y.dat", "room12VX.dat", "room12VY.dat", "room12Wall.dat");

    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);

    std::cout << "total elapsed time " << duration.count() << " seconds" << std::endl;
    std::cout << "Physical time (seconds) " << nbTseps * args.dt << std::endl;
    return 0;
}