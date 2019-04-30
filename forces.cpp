#include "forces.h"
#include <vector>
#include "math.h"
#include <cmath>
#include <iostream>
#include <algorithm>
#include "polygonalWall.h"
#include <fstream>
#include <sstream>
#include <random>
#include <chrono>

#define MAXIMAL_FORCE 10e80

double Theta(double cos_phi_max, double lambda, double nX, double nY, double rijX, double rijY) {
    double prodScal = nX * rijX + nY * rijY;
    double res = prodScal < cos_phi_max ? lambda+(1-lambda)*0.5*(1+prodScal) : 1;
    return res;
}

double g(double x) {
    return x < 0 ? 0 : x;
}

double collisionDistanceIndividual(double x0, double y0, double vx0, double vy0, double x1, double y1,
    double vx1, double vy1, double maxDistance, double Radius0, double Radius1) {
        double distance;
        double A = (vx1 - vx0) * (vx1 - vx0) + (vy1 - vy0) * (vy1 - vy0);
        if ( abs(A) < epss ) return maxDistance;
        double B = 2 * (vx1 - vx0) * (x1 -x0) + 2 * (vy1 - vy0) * (y1 -y0);
        double C = (x1 -x0) * (x1 - x0) + (y1 -y0) * (y1 - y0) - (Radius0 + Radius1) * (Radius0 + Radius1);
        double Delta = B * B - 4 * A * C;
        double sqrtDelta = sqrt(Delta);
        double t1 = (-B - sqrtDelta) / (2 * A);
        double t2 = (-B + sqrtDelta) / (2 * A);
        double deltaT = std::min(t1, t2) <= 0 ? std::max(t1, t2) : std::min(t1, t2);
        distance = deltaT > 0 ? sqrt(vx0 * vx0 + vy0 * vy0) * deltaT : maxDistance;
        distance = std::min(distance, maxDistance);
        return distance;
    }

double collisionDistanceWall(double x0, double y0, double vx0, double vy0, double Radius0, double maxDistance, Wall &wall) {
    double vNorm = sqrt(vx0*vx0 + vy0*vy0);
    double cosAlpha = vx0 / vNorm;
    double sinAlpha = vy0 / vNorm;
    std::vector<polygonalWall> *walls = &wall.walls_;
    double minCollisionDistance = maxDistance;
    for(int idx = 0; idx < walls->size(); idx++) {
        std::vector<double> verticesX = (*walls)[idx].getVerticesX();
        std::vector<double> verticesY = (*walls)[idx].getVerticesY();
        for(int j=0; j<verticesX.size() - 1; j++) {
            double xw1 = verticesX[j];
            double xw2 = verticesX[j+1];
            double yw1 = verticesY[j];
            double yw2 = verticesY[j+1];
            double a,b,c;
            if(abs(xw2 - xw1) > epss) {
                double slope = (yw2 - yw1) / (xw2 - xw1);
                a = slope;
                b = -1;
                c = -slope * xw1 + yw1;
            }
            else {
                a = 1;
                b = 0;
                c = -xw1; 
            }
            double distanceToWall;
            if (abs(a * cosAlpha + b * sinAlpha) < epss) distanceToWall = maxDistance;
            else {
                double delta = a * cosAlpha + b * sinAlpha;
                double xCol = (-c * cosAlpha - b * (y0 * cosAlpha - x0 * sinAlpha)) / delta;
                double ycol = (a * (y0 * cosAlpha - x0 * sinAlpha) - c * sinAlpha) / delta;

                if ((xCol >= std::min(xw1, xw2) - Radius0) && (xCol <= std::max(xw1, xw2) + Radius0) &&
                    (ycol >= std::min(yw1, yw2) - Radius0) && (ycol <= std::max(yw1, yw2) + Radius0) && 
                    (xCol - x0) * vx0 + (ycol - y0) * vy0 >= 0)
                    distanceToWall = sqrt ((xCol - x0) * (xCol - x0) + (ycol - y0) * (ycol - y0)) - Radius0;
                else distanceToWall = maxDistance;
            }

            distanceToWall = std::min(maxDistance, distanceToWall);
            minCollisionDistance = std::min(minCollisionDistance, distanceToWall);
        }
    }
    return minCollisionDistance;
}

std::vector<double> desiredSpeedMotivationSFM (Arguments &args, Wall &wall) {
    long i = args.iIndividual;
    double mi = args.mass->coeff(i);
    double resX = mi / args.relaxationTime->coeff(i) * (args.desiredSpeed->coeff(i) * 
        args.normalX->coeff(args.currentTimeStep, i) - args.velocityX->coeff(args.currentTimeStep, i));
    double resY = mi / args.relaxationTime->coeff(i) * (args.desiredSpeed->coeff(i) * 
        args.normalY->coeff(args.currentTimeStep, i) - args.velocityY->coeff(args.currentTimeStep, i));
    std::vector<double> res{resX, resY};
    return res;
}

std::vector<double> desiredSpeedMotivationHeuristics(Arguments &args, Wall &wall){
    long i = args.iIndividual;
    double mi = args.mass->coeff(i);
    double velocityNorm = sqrt(pow(args.velocityX->coeff(args.currentTimeStep, i), 2) + 
                        pow(args.velocityY->coeff(args.currentTimeStep, i), 2));
    while(velocityNorm < epss) {
        // singular case, when the pedestrian is idle
        // Impossible to compute is moving angle
        // we force him to have a slow motion in a random direction
        std::uniform_real_distribution<double> distributionV(-epss, epss);
        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
        std::default_random_engine generator (seed);
        args.velocityX->coeffRef(args.currentTimeStep, i) = distributionV(generator);
        args.velocityY->coeffRef(args.currentTimeStep, i) = distributionV(generator);
        velocityNorm = sqrt(pow(args.velocityX->coeff(args.currentTimeStep, i), 2) + 
                        pow(args.velocityY->coeff(args.currentTimeStep, i), 2));
    }
    double actualDirectionX = args.velocityX->coeff(args.currentTimeStep, i) / velocityNorm;
    double actualDirectionY = args.velocityY->coeff(args.currentTimeStep, i) / velocityNorm;
    double actualAngle = acos(actualDirectionX);
    
    if (actualDirectionY < 0.0) actualAngle *= -1;

    double targetAngle = acos(args.normalX->coeff(args.currentTimeStep, i));
    if (args.normalY->coeff(args.currentTimeStep, i) < 0 ) targetAngle *= -1;

    double minVisionAngle = actualAngle - args.sightAngleRange;
    double maxVisionAngle = actualAngle + args.sightAngleRange;

    int movingAngleResolution = args.movingAngleResolution;
    std::vector<double> potentialNextAngle(movingAngleResolution);
    std::vector<double> distanceToCollision(movingAngleResolution, args.sightRangeDistance);
    for(int iAngle=0; iAngle<movingAngleResolution; iAngle++){
        potentialNextAngle[iAngle] = minVisionAngle + (maxVisionAngle-minVisionAngle) * 
            iAngle / (movingAngleResolution - 1);
    }
    for(int idx=0; idx<movingAngleResolution; idx++) {
        // detect potential collision
        double minimalColDistance = args.sightRangeDistance;
        for(int jIndividual=0; jIndividual<args.nbIndividual; jIndividual++){
            if(jIndividual != i){
                double colDistanceIndividual = collisionDistanceIndividual(args.positionX->coeff(args.currentTimeStep, i),
                                                       args.positionY->coeff(args.currentTimeStep, i),
                                                       args.desiredSpeed->coeff(i) * cos(potentialNextAngle[idx]),
                                                       args.desiredSpeed->coeff(i) * sin(potentialNextAngle[idx]),
                                                       args.positionX->coeff(args.currentTimeStep, jIndividual),
                                                       args.positionY->coeff(args.currentTimeStep, jIndividual),
                                                       args.velocityX->coeff(args.currentTimeStep, jIndividual),
                                                       args.velocityY->coeff(args.currentTimeStep, jIndividual),
                                                       args.sightRangeDistance, 
                                                       args.Radius->coeff(i), 
                                                       args.Radius->coeff(jIndividual));
                minimalColDistance = std::min(minimalColDistance, colDistanceIndividual);
            }
        }
        double colDistanceWall = collisionDistanceWall(args.positionX->coeff(args.currentTimeStep, i),
                                                       args.positionY->coeff(args.currentTimeStep, i),
                                                       args.desiredSpeed->coeff(i) * cos(potentialNextAngle[idx]),
                                                       args.desiredSpeed->coeff(i) * sin(potentialNextAngle[idx]),
                                                       args.Radius->coeff(i),
                                                       args.sightRangeDistance, 
                                                       wall);
        minimalColDistance = std::min(minimalColDistance, colDistanceWall);                               
        distanceToCollision[idx] = minimalColDistance;
    }
    
    DVector distanceAvoidingObstacles(movingAngleResolution);
    for(int idx=0; idx<movingAngleResolution; idx++) {
        distanceAvoidingObstacles[idx] = args.sightRangeDistance * args.sightRangeDistance +
                                       distanceToCollision[idx] * distanceToCollision[idx] -
                                       2.0 * args.sightRangeDistance * distanceToCollision[idx] *
                                       cos(potentialNextAngle[idx] - targetAngle);
    }
    int desiredAngleIdx;
    double minObstacleDistance = distanceAvoidingObstacles.minCoeff(&desiredAngleIdx);
    double desiredAngle { potentialNextAngle[desiredAngleIdx] };
    double desiredDirectionX = cos(desiredAngle);
    double desiredDirectionY = sin(desiredAngle);
    
    double desiredSpeed = std::min(args.desiredSpeed->coeff(i), distanceToCollision[desiredAngleIdx] / args.relaxationTime->coeff(i));

    double resX = mi / args.relaxationTime->coeff(i) * (desiredSpeed * desiredDirectionX - 
        args.velocityX->coeff(args.currentTimeStep, i));
    double resY = mi / args.relaxationTime->coeff(i) * (desiredSpeed * desiredDirectionY -
        args.velocityY->coeff(args.currentTimeStep, i));
    std::vector<double> res{resX, resY};

    if(args.nbIndividual <= 2) {
        std::ostringstream oss;
        oss << "distance_" << args.iIndividual << "_" << args.currentTimeStep << ".dat";
        std::ofstream ofs(oss.str().c_str());
        for(int idx=0; idx<movingAngleResolution; idx++) {
            ofs << potentialNextAngle[idx] << "\t" << distanceToCollision[idx] << std::endl;
        }
    }
    
    return res;
}

std::vector<double> wallRepulsionSFM (Arguments &args, Wall &wall) {
    long i = args.iIndividual;
    double Ri = args.Radius->coeff(i);
    std::vector<double> res(2);
    for(int iWall=0; iWall <args.wallSize; iWall++) {
        args.iWall = iWall;
        double diW = args.diW->coeff(i, args.iWall);
        double normalForceNorm = args.Ai * exp((Ri - diW)/args.Bi) + args.ki * g(Ri - diW);
        double prodScal = -args.velocityX->coeff(args.currentTimeStep, i) * args.riWY->coeff(i, args.iWall) + 
                       args.velocityY->coeff(args.currentTimeStep, i) * args.riWX->coeff(i, args.iWall);
        double tangentialForceNorm = args.kappai * g(Ri - diW)* prodScal;
    
        res[0] += normalForceNorm * args.riWX->coeff(i, args.iWall);
        res[0] += tangentialForceNorm * args.riWY->coeff(i, args.iWall);
        res[1] += normalForceNorm * args.riWY->coeff(i, args.iWall);
        res[1] -= tangentialForceNorm * args.riWX->coeff(i, args.iWall);
    }
    return res;
}

std::vector<double> wallRepulsionHeuristics (Arguments &args, Wall &wall) {
    long i = args.iIndividual;
    double Ri = args.Radius->coeff(i);
    std::vector<double> res(2);
    for(int iWall=0; iWall <args.wallSize; iWall++) {
        args.iWall = iWall;
        double diW = args.diW->coeff(i, args.iWall);
        double normalForceNorm = args.ki * g(Ri - diW);
    
        res[0] += normalForceNorm * args.riWX->coeff(i, args.iWall);
        res[1] += normalForceNorm * args.riWY->coeff(i, args.iWall);
    }
    return res;
}

std::vector<double> individualRepulsionCrowdTurbulence (Arguments &args, Wall &wall) {
    std::vector<double> res(2);
    for(int jIndividual = 0; jIndividual<args.nbIndividual; jIndividual++) {
        if (jIndividual == args.iIndividual) continue;
        args.jIndividual = jIndividual;
        double rijX = args.rijX->coeff(args.iIndividual, args.jIndividual);
        double rijY = args.rijY->coeff(args.iIndividual, args.jIndividual);
        double theta = Theta(cos(args.sightAngleRange), args.lambda, args.normalX->coeff(args.currentTimeStep, args.iIndividual), 
            args.normalY->coeff(args.currentTimeStep, args.iIndividual), rijX, rijY);
        double dij = args.dij->coeff(args.iIndividual, args.jIndividual);
        double norm = args.maximalRepulsionForce * theta * exp(-dij/args.D0 + pow(args.D1 / dij, args.k));
        if (isinf(norm) || isnan(norm) || norm > MAXIMAL_FORCE) norm = MAXIMAL_FORCE;
        res[0] += norm * rijX;
        res[1] += norm * rijY;
    }
    return res;
}

std::vector<double> individualRepulsionSFM (Arguments &args, Wall &wall) {
    std::vector<double> res(2);
    double Ri = args.Radius->coeff(args.iIndividual);
    for(int jIndividual = 0; jIndividual < args.nbIndividual; jIndividual++) {
        if (jIndividual == args.iIndividual) continue;
        args.jIndividual = jIndividual;
        double Rj = args.Radius->coeff(jIndividual);
        double rijX = args.rijX->coeff(args.iIndividual, args.jIndividual);
        double rijY = args.rijY->coeff(args.iIndividual, args.jIndividual);
        double nX = args.normalX->coeff(args.currentTimeStep, args.iIndividual);
        double nY = args.normalY->coeff(args.currentTimeStep, args.iIndividual);
        double prodScal = nX * rijX + nY * rijY;
        double viewingFactor = args.lambda + (1 - args.lambda) * 0.5 * (1 + prodScal);
        double dij = args.dij->coeff(args.iIndividual, args.jIndividual);
        double normalForceNorm = args.Ai * exp((Ri + Rj - dij)/args.Bi) * viewingFactor + args.ki * g(Ri + Rj - dij);
        double DeltaVscalTij = (args.velocityX->coeff(args.currentTimeStep, args.iIndividual) - 
                                args.velocityX->coeff(args.currentTimeStep, args.jIndividual)) * rijY + 
                               (args.velocityY->coeff(args.currentTimeStep, args.jIndividual) - 
                                args.velocityY->coeff(args.currentTimeStep, args.iIndividual)) * rijX;
        double tangentialForceNorm = args.kappai * g(Ri + Rj - dij) * DeltaVscalTij;
        //double tangentialForceNorm = 0.0;
        res[0] += normalForceNorm * rijX;
        res[0] -= tangentialForceNorm * rijY;
        res[1] += normalForceNorm * rijY;
        res[1] += tangentialForceNorm * rijX;
    }
    return res;
}

std::vector<double> individualRepulsionHeuristics(Arguments &args, Wall &wall) {
    long i = args.iIndividual;
    double Ri = args.Radius->coeff(i);
    std::vector<double> res(2);
    for(int jIndividual=0; jIndividual < args.nbIndividual; jIndividual++) {
        if (jIndividual == i) continue;
        args.jIndividual = jIndividual;
        double Rj = args.Radius->coeff(jIndividual);
        double rijX = args.rijX->coeff(args.iIndividual, args.jIndividual);
        double rijY = args.rijY->coeff(args.iIndividual, args.jIndividual);
        double dij = args.dij->coeff(args.iIndividual, args.jIndividual);
        double normalForceNorm = args.ki * g(Ri + Rj - dij);
    
        res[0] += normalForceNorm * rijX;
        res[1] += normalForceNorm * rijY;
    }
    return res;
}