#pragma once
#include <Eigen/Dense>

using namespace Eigen;
typedef Matrix<double, Dynamic, Dynamic, RowMajor> DMatrix;
typedef Matrix<double, 1, Dynamic> DVector;

struct Arguments {
    double maxSpeedFactor;
    long nbIndividual;
    int nbTimeSteps;
    long iIndividual, jIndividual;
    int iWall, wallSize;
    int currentTimeStep;
    double dt;
    double lambda;
    double maximalRepulsionForce;
    double k, D0, D1;
    double Ai, Bi, ki, kappai;
    DMatrix *velocityX, *velocityY, *positionX, *positionY;
    DVector *desiredSpeed, *maxSpeed;
    DVector *relaxationTime;
    DMatrix *normalX, *normalY;
    DMatrix *rijX, *rijY;
    DMatrix *dij;
    DMatrix *riWX, *riWY;
    DMatrix *diW;
    DVector *Radius;
    DVector *mass;
    double sightAngleRange;
    int movingAngleResolution;
    double sightRangeDistance;
    int maxStepsInPath;
};
