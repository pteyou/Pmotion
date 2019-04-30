#pragma once
#include <string>
#include <fstream>
#include "assert.h"
#include "totalForces.h"
#include "myTypes.h"
#include "math.h"
#include <cmath>
#include "polygonalWall.h"
#include <random>
#include <chrono>
#include <omp.h>
#include <sstream>
#include <fstream>
#include <csignal>
#include "FMM/FMM.h"
#include "utilities.h"
#include "cases.h"
#include <set>
#define loin 1e5
#define deb false

//#define DistanceEps 0.5
#define DistanceEps 0.3     // reunion
//#define DistanceEps 1     // cinema

double L2norm(double x, double y){
    return sqrt(x*x + y*y);
}

std::vector<int> crossedWallPointsRegul (const DVector &initX, const DVector &initY, double *finalX, 
    double *finalY, const DVector &initVX, const DVector &initVY, double *finalVX, double *finalVY, DVector *radius,
    Wall &wall) {
    std::vector<int> res;
    std::vector<polygonalWall> *walls = &wall.walls_;
    std::vector<std::vector<double>> wallVerticesX; wallVerticesX.reserve(wall.nbLines());
    std::vector<std::vector<double>> wallVerticesY; wallVerticesY.reserve(wall.nbLines());
    for(int idx = 0; idx < walls->size(); idx++) {
        wallVerticesX.push_back((*walls)[idx].getVerticesX());
        wallVerticesY.push_back((*walls)[idx].getVerticesY());
    }
    int nbPoint = initX.size();
    #pragma omp parallel for
    for(int iPoint=0; iPoint<nbPoint; iPoint++) {
        double Xinit = initX[iPoint];
        double Yinit = initY[iPoint];
        double Xfinal = finalX[iPoint];
        double Yfinal = finalY[iPoint];
        double VinitX = initVX[iPoint];
        double VinitY = initVY[iPoint];
        double vNorm = L2norm(VinitX, VinitY);
        while(vNorm < epss) {
            // singular case, when the pedestrian is idle
            // Impossible to compute is moving angle
            // we force him to have a slow motion in a random direction
            std::uniform_real_distribution<double> distributionV(-epss, epss);
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine generator (seed);
            VinitX = distributionV(generator);
            VinitY = distributionV(generator);
            vNorm = L2norm(VinitX, VinitY);
        }
        double cosAlpha = VinitX / vNorm;
        double sinAlpha = VinitY / vNorm; 
        double alpha = acos(cosAlpha);
        if(sinAlpha < 0) alpha *= -1; 
        double Radius0 = radius->coeff(iPoint);
        Xfinal += Radius0 * cosAlpha;
        Yfinal += Radius0 * sinAlpha;
        for(int jList=0; jList<wallVerticesX.size(); jList++) {
            for (int idx=0; idx < wallVerticesX[jList].size() - 1; idx++) {
                double xw1 = wallVerticesX[jList][idx];
                double xw2 = wallVerticesX[jList][idx + 1];
                double yw1 = wallVerticesY[jList][idx];
                double yw2 = wallVerticesY[jList][idx + 1];
                double delta = (xw2 - xw1) * (Yinit - Yfinal) - (Xinit - Xfinal) * (yw2 - yw1);
                double dist = L2norm(yw2 - yw1, xw2 - xw1);
                double alphaW = abs(xw2 - xw1) < epss ? M_PI * 0.5 : acos ((xw2 - xw1) / dist);
                if (yw2 - yw1 < 0) alphaW *= -1;
                if (abs(delta) < epss) continue; // colinear
                else {
                    double ta = ( (yw1 - yw2) * (Xinit - xw1) + (xw2 - xw1) * (Yinit - yw1) ) / delta;
                    double tb = ( (Yinit - Yfinal) * (Xinit - xw1) + (Xfinal - Xinit) * (Yinit - yw1) ) / delta;
                    double xCol = Xinit + ta * (Xfinal - Xinit);
                    double ycol = Yinit + ta * (Yfinal - Yinit);

                    if (ta >= 0 && ta <= 1 && tb >= 0 && tb <= 1 && (xCol - Xinit) * VinitX + (ycol - Yinit) * VinitY >= 0) {
                        finalX[iPoint] = abs(xCol - Xinit) < epss ? xCol - Radius0 * cosAlpha : 
                            Xinit + (ta - Radius0 * cosAlpha / (xCol - Xinit)) * (xCol - Xinit);
                        finalY[iPoint] =  abs(ycol - Yinit) < epss ? ycol - Radius0 * sinAlpha : 
                            Yinit + (ta - Radius0 * sinAlpha / (ycol - Yinit)) * (ycol - Yinit);
                        finalVX[iPoint] = vNorm * cos(-alpha + 2.0 * alphaW);
                        finalVY[iPoint] = vNorm * sin(-alpha + 2.0 * alphaW);
                        #pragma omp critical
                        res.push_back(iPoint);
                    }
                }
            }
        }
    }
    return res;
}


class EulerBase {
    friend struct Arguments;
public:
    EulerBase(totalForces &forces, Wall &wall, Arguments& args) : nbPoints_(args.nbIndividual),
     nbTimeSteps_(args.nbTimeSteps), timeStep_(args.dt), forces_(forces), wall_(wall),
     args_(args) {
        me = this;
        nbEjectedPoints_ = 0;
        nbCasualties_ = 0;
        velocityX_.resize(nbTimeSteps_ + 1, nbPoints_);
        velocityY_.resize(nbTimeSteps_ + 1, nbPoints_);
        positionX_.resize(nbTimeSteps_ + 1, nbPoints_);
        positionY_.resize(nbTimeSteps_ + 1, nbPoints_);
        normalX.resize(nbTimeSteps_ + 1, nbPoints_);
        normalY.resize(nbTimeSteps_ + 1, nbPoints_);
        currentTimeStep = 0;
        args_.positionX = &positionX_;
        args_.positionY = &positionY_;
        args_.velocityX = &velocityX_;
        args_.velocityY = &velocityY_;
        args_.normalX = &normalX;
        args_.normalY = &normalY;
        args_.rijX = &rijX_;
        args_.rijY = &rijY_;
        args_.dij = &dij_;
        args_.riWX = &riWX_;
        args_.riWY = &riWY_;
        args_.diW = &diW_;
    }
    void initialization(DVector &initX, DVector &initY, DVector &initVX, DVector &initVY, 
        DMatrix &finalX, DMatrix &finalY) {
        velocityX_.row(0) = initVX;
        velocityY_.row(0) = initVY;
        positionX_.row(0) = initX;
        positionY_.row(0) = initY;
        desiredPositionX_ = finalX;
        desiredPositionY_ = finalY;
        normalX.row(0) = desiredPositionX_.row(0) - initX;
        normalY.row(0) = desiredPositionY_.row(0) - initY;
        #pragma omp parallel for
        for(int i=0; i<nbPoints_; i++){
            double norm = L2norm(normalX(0, i), normalY(0, i));
            normalX(0,i) /= norm;
            normalY(0,i) /= norm;
        }
        
        rijX_.resize(nbPoints_, nbPoints_);
        rijY_.resize(nbPoints_, nbPoints_);
        dij_.resize(nbPoints_, nbPoints_);
        #pragma omp parallel for
        for(int i=0; i<nbPoints_; i++) {
            for(int j=0; j<nbPoints_; j++) {
                if(i != j) {
                    rijX_(i,j) = positionX_(0, i) - positionX_(0, j);
                    rijY_(i,j) = positionY_(0, i) - positionY_(0, j);
                    dij_(i,j) = L2norm(rijX_(i,j), rijY_(i, j));
                    #pragma omp critical
                    while(dij_(i, j) < epss) {
                        // avoid nan propagation when compression is huge
                        positionX_(0, i) += epss;
                        rijX_(i,j) = positionX_(0, i) - positionX_(0, j);
                        dij_(i,j) = L2norm(rijX_(i,j), rijY_(i, j));
                    }
                    rijX_(i,j) /= dij_(i,j);
                    rijY_(i,j) /= dij_(i,j);
                }
            }
        }
        riWX_.resize(nbPoints_, wall_.size());
        riWY_.resize(nbPoints_, wall_.size());
        diW_.resize(nbPoints_, wall_.size());
        wall_.storeCoordiantes();
        #pragma omp parallel for
        for(int i=0; i<nbPoints_; i++) {
            for(int j = 0; j<wall_.size(); j++) {
                riWX_(i, j) = positionX_(0, i) - wall_.wallX_[j];
                riWY_(i, j) = positionY_(0, i) - wall_.wallY_[j];
                diW_(i, j) = L2norm(riWX_(i, j), riWY_(i, j));
                riWX_(i,j) /= diW_(i, j);
                riWY_(i,j) /= diW_(i, j);
            }
        }
        currentStep.resize(nbPoints_);
    }
    
    virtual void step() = 0;

    void updateWallData(){
        #pragma omp parallel for
        for(int i=0; i<nbPoints_; i++) {
            for(int j = 0; j<wall_.size(); j++) {
                riWX_(i, j) = positionX_(currentTimeStep, i) - wall_.wallX_[j];
                riWY_(i, j) = positionY_(currentTimeStep, i) - wall_.wallY_[j];
                diW_(i, j) = L2norm(riWX_(i, j), riWY_(i, j));
                riWX_(i, j) /= diW_(i, j);
                riWY_(i, j) /= diW_(i, j);
            }
        }
    }

    void updateCasualties(const std::string CasualtiesFname = "casualtiesRecord.dat"){  
        std::ofstream ofs;
        if (currentTimeStep <= 1){
            ofs.open(CasualtiesFname.c_str(),  std::ios::trunc);
            ofs.close();
        }
        else if (!deathPeopleIdx_.empty()) {
            ofs.open(CasualtiesFname.c_str(), std::ios::app);
            std::vector<int> nbSubPoints(4, 5);
            std::uniform_real_distribution<double> distributionPX(loin-10000, loin+10000);
            unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
            std::default_random_engine generatorX (seed);
            std::default_random_engine generatorY(seed+10);
            for(auto d : deathPeopleIdx_){
                ofs << currentTimeStep << "\t" << positionX_(currentTimeStep, d) << "\t" <<
                    positionY_(currentTimeStep, d) << std::endl;
                
                double radius = args_.Radius->coeff(d);
                std::vector<double> X(5);
                X[0] = positionX_(currentTimeStep, d) - radius / 2;
                X[1] = X[0];
                X[2] = X[0] + radius;
                X[3] = X[2];
                X[4] = X[0];
                std::vector<double> Y(5);
                Y[0] = positionY_(currentTimeStep, d) - radius / 2;
                Y[1] = Y[0] + radius;
                Y[2] = Y[1];
                Y[3] = Y[0];
                Y[4] = Y[0];
                polygonalWall w(X,Y,nbSubPoints);
                wall_.addPolygonalWall(w);
                // release the point
                positionX_(currentTimeStep, d) = distributionPX(generatorX);
                positionY_(currentTimeStep, d) = distributionPX(generatorY);
                desiredPositionX_(0, d) = loin;
                desiredPositionY_(0, d) = loin;
                #pragma omp critical
                {
                    nbEjectedPoints_++;
                    nbCasualties_++;
                    if (nbEjectedPoints_ % 50 == 0) std::cout << "***** " << nbEjectedPoints_ << " points ejected **********" << std::endl;
                }
            }
            deathPeopleIdx_.clear();
            ofs.close();
            riWX_.resize(nbPoints_, wall_.size());
            riWY_.resize(nbPoints_, wall_.size());
            diW_.resize(nbPoints_, wall_.size());
            args_.wallSize = wall_.size();
            updateWallData();
        }
    }

    int run (DVector &initX, DVector &initY, DVector &initVX, DVector &initVY, 
        DMatrix &finalX, DMatrix &finalY) {
        initialization(initX, initY, initVX, initVY, finalX, finalY);
        int i = 1;
        while(i<=nbTimeSteps_ && !domainIsEmpty_()) {
            if (i % 50 == 0) std::cout << "********* computing time step " << i << " *************" << std::endl;
            step();

            updateCasualties();

            if (i % 100 == 0 && deb) { // partial output for debug purposes
                // partial output
                std::cout << "********* storing partial output at step " << i << " *************" << std::endl;
                std::ostringstream osX, osY, osVx, osVy, osW;
                osX << "X_" << i << ".dat";
                osY << "Y_" << i << ".dat";
                osVx << "VX_" << i << ".dat";
                osVy << "VY_" << i << ".dat";
                osW << "wall_" << i << ".dat";
                output(osX.str(), osY.str(), osVx.str(), osVy.str(), osW.str());
            }
            i++;
        }
        return --i;
    }

    void output (std::string positionFnameX, std::string positionFnameY, 
        std::string velocityFnameX, std::string velocityFnameY, std::string wallFname) {
        std::ofstream pfsX(positionFnameX.c_str());
        for(int i=0; i<=currentTimeStep; i++) {
            pfsX << timeStep_ * i << "\t" << positionX_.row(i) << std::endl;
        }
        pfsX.close();
        std::ofstream pfsY(positionFnameY.c_str());
        for(int i=0; i<=currentTimeStep; i++) {
            pfsY << timeStep_ * i << "\t" << positionY_.row(i) << std::endl;
        }
        pfsY.close();
        std::ofstream vfsX(velocityFnameX.c_str());
        for(int i=0; i<=currentTimeStep; i++) {
            vfsX << timeStep_ * i << "\t" << velocityX_.row(i) << std::endl;
        }
        vfsX.close();
        std::ofstream vfsY(velocityFnameY.c_str());
        for(int i=0; i<=currentTimeStep; i++) {
            vfsY << timeStep_ * i << "\t" << velocityY_.row(i) << std::endl;
        }
        vfsY.close();
        std::ofstream wallf(wallFname.c_str());  
        for(int i=0; i<wall_.size(); i++){
            std::vector<double> point{wall_.point(i)};
            wallf << point[0]  << "\t" << point[1] << std::endl;
        }
        wallf.close();
    }
    static void signalHandler(int x) {
        std::cout << "Interrupt signal (" << x << ") received.\n";
        std::cout << "saving current state before leaving" << std::endl;
        if (x == SIGINT) {
            me->output("IX.dat", "IY.dat", "IVX.dat", "IVY.dat", "IWall.dat");
        }
        exit(x);
    }

protected:
    long nbPoints_;
    volatile long nbEjectedPoints_;
    volatile long nbCasualties_;
    int nbTimeSteps_;
    double timeStep_;
    totalForces forces_;
    DMatrix velocityX_, velocityY_, positionX_, positionY_;
    DMatrix desiredPositionX_, desiredPositionY_;
    DMatrix normalX, normalY;
    DMatrix rijX_, rijY_, dij_, riWX_, riWY_, diW_;
    int currentTimeStep;
    Wall &wall_;
    Arguments args_;
    DVector currentStep;
    std::vector<int> deathPeopleIdx_;
    static EulerBase* me;
    virtual bool domainIsEmpty_() = 0; 
};

class Euler : public EulerBase {
protected:
    bool domainIsEmpty_() override {
        return nbEjectedPoints_ == nbPoints_;
    }
public:
    Euler(totalForces &forces, Wall &wall, Arguments& args) : 
        EulerBase(forces, wall, args) {}
    void step() {
        args_.currentTimeStep = currentTimeStep;
        assert(currentTimeStep <= nbTimeSteps_);
        #pragma omp parallel for
        for(int iPoint = 0; iPoint < nbPoints_; iPoint++) {
            args_.iIndividual = iPoint;
            std::vector<double> f { forces_.computeTotalForces(args_, wall_) };
            velocityX_(currentTimeStep + 1, iPoint) = velocityX_(currentTimeStep,  iPoint) + 
                timeStep_ / args_.mass->coeff(iPoint) * f[0];
            velocityY_(currentTimeStep + 1, iPoint) = velocityY_(currentTimeStep,  iPoint) + 
                timeStep_ / args_.mass->coeff(iPoint) * f[1];
        }
        currentTimeStep++;

        // Regularization of velocity
        // Must have a norm lower to Vmax
        #pragma omp parallel for
        for(int iPoint = 0; iPoint < nbPoints_; iPoint++) {
            double vNorm = L2norm(velocityX_(currentTimeStep, iPoint), velocityY_(currentTimeStep, iPoint));
            double vMax =  args_.maxSpeed->coeff(iPoint);
            if (vNorm > vMax){
                velocityX_(currentTimeStep, iPoint) *= vMax / vNorm;
                velocityY_(currentTimeStep, iPoint) *= vMax / vNorm;
            }
        }

        
        positionX_.row(currentTimeStep) = positionX_.row(currentTimeStep - 1) + 
                velocityX_.row(currentTimeStep) * timeStep_;
        positionY_.row(currentTimeStep) = positionY_.row(currentTimeStep - 1) + 
                velocityY_.row(currentTimeStep) * timeStep_;


        // check if the walll has been crossed and regularize if it is the case
        // application of Ping pong rule
        std::vector<int> reIntegratedPoints {crossedWallPointsRegul(
                                                positionX_.row(currentTimeStep-1),
                                                positionY_.row(currentTimeStep-1),
                                                positionX_.row(currentTimeStep).data(),
                                                positionY_.row(currentTimeStep).data(),
                                                velocityX_.row(currentTimeStep-1),
                                                velocityY_.row(currentTimeStep-1),
                                                velocityX_.row(currentTimeStep).data(),
                                                velocityY_.row(currentTimeStep).data(),
                                                args_.Radius,
                                                wall_                        )};
        if(!reIntegratedPoints.empty()) {
            std::cout << "**** Reintegrated " << reIntegratedPoints.size() << 
                " Points after wall crossing **** " << std::endl;
        }

        // check if one step passed
        if(args_.maxStepsInPath > 1) {
            #pragma omp parallel for
            for(int idx=0; idx < nbPoints_; idx++) {
                double distToObjective = L2norm(desiredPositionX_(0, idx) - positionX_(currentTimeStep, idx), 
                                                desiredPositionY_(0, idx) - positionY_(currentTimeStep, idx));
                if (distToObjective < DistanceEps) {
                    if(currentStep[idx] < args_.maxStepsInPath - 1) {
                        desiredPositionX_(0, idx) = desiredPositionX_(currentStep[idx] + 1, idx);
                        desiredPositionY_(0, idx) = desiredPositionY_(currentStep[idx] + 1, idx );
                        currentStep[idx] += 1;
                    }
                    else { // release the point
                        std::uniform_real_distribution<double> distributionPX(loin-10000, loin+10000);
                        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                        std::default_random_engine generatorX (seed);
                        std::default_random_engine generatorY(seed+10);
                        positionX_(currentTimeStep, idx) = distributionPX(generatorX);
                        positionY_(currentTimeStep, idx) = distributionPX(generatorY);
                        desiredPositionX_(0, idx) = loin;
                        desiredPositionY_(0, idx) = loin;
                        nbEjectedPoints_++;
                    }
                    
                }
            }
        }
        
        normalX.row(currentTimeStep) = desiredPositionX_.row(0) - positionX_.row(currentTimeStep);
        normalY.row(currentTimeStep) = desiredPositionY_.row(0) - positionY_.row(currentTimeStep);
        #pragma omp parallel for
        for(int i=0; i<nbPoints_; i++){
            double norm = L2norm(normalX(currentTimeStep, i), normalY(currentTimeStep, i));
            if (norm < epss) std::cout << "direction reg" << std::endl;
            while(norm < epss) {
                // avoid nan propagation when a point is arrived
                desiredPositionX_(0, i) += epss;
                normalX(currentTimeStep, i) = desiredPositionX_(0, i) - positionX_(currentTimeStep, i);
                norm = L2norm(normalX(currentTimeStep, i), normalY(currentTimeStep, i));
            }
            normalX(currentTimeStep,i) /= norm;
            normalY(currentTimeStep,i) /= norm;
        }

        #pragma omp parallel for
        for(int i=0; i<nbPoints_; i++) {
            for(int j=0; j<nbPoints_; j++) {
                if(i != j) {
                    rijX_(i,j) = positionX_(currentTimeStep, i) - positionX_(currentTimeStep, j);
                    rijY_(i,j) = positionY_(currentTimeStep, i) - positionY_(currentTimeStep, j);
                    dij_(i,j) = L2norm(rijX_(i,j), rijY_(i, j));
                    if (dij_(i,j) < epss) {
                        std::cout << "position reg" << std::endl;
                        //output("Xreprise.dat", "Yreprise.dat", "VXReprise.dat", "VYReprise.dat", "wall.dat");
                        //int cmp = 0;
                        #pragma omp critical
                        while(dij_(i, j) < epss) {
                            // avoid nan propagation when compression is huge
                            positionX_(currentTimeStep, i) += epss;
                            rijX_(i,j) = positionX_(currentTimeStep, i) - positionX_(currentTimeStep, j);
                            dij_(i, j) = L2norm(rijX_(i, j), rijY_(i, j));
                            //if(cmp++ < 20) std::cout << rijX_(i, j) << " " << dij_(i, j) << std::endl;
                        }
                    }
                    rijX_(i,j) /= dij_(i,j);
                    rijY_(i,j) /= dij_(i,j);
                }
            }
        }
        updateWallData();
    }
};

template<class Case>
class EulerFMM : public EulerBase {
private:
    grid g_;
protected:
    bool domainIsEmpty_() override {
        return Case::domainIsEmpty(nbEjectedPoints_, nbPoints_);
    }
public:
    EulerFMM(totalForces &forces, Wall &wall, Arguments& args, double dx, double dy) : 
        EulerBase(forces, wall, args) {
            g_ = Case::initializeGrid(wall_, dx, dy);
        }
    void step() {
        deathPeopleIdx_.clear();
        args_.currentTimeStep = currentTimeStep;
        assert(currentTimeStep <= nbTimeSteps_);
        #pragma omp parallel for
        for(int iPoint = 0; iPoint < nbPoints_; iPoint++) {
            // la ligne suivante est la cause du problem avec openMP -- à fixer
            args_.iIndividual = iPoint;
            std::vector<double> f { forces_.computeTotalForces(args_, wall_) };
            velocityX_(currentTimeStep + 1, iPoint) = velocityX_(currentTimeStep,  iPoint) + 
                timeStep_ / args_.mass->coeff(iPoint) * f[0];
            velocityY_(currentTimeStep + 1, iPoint) = velocityY_(currentTimeStep,  iPoint) + 
                timeStep_ / args_.mass->coeff(iPoint) * f[1];

            /*
            double normF = L2norm(f[0], f[1]);
            if (normF/2.0/M_PI/args_.Radius->coeff(iPoint) > 1600 && 
                positionX_(currentTimeStep, iPoint) <= g_.xMin_ + g_.dx_ * g_.iMax_)
                deathPeopleIdx_.push_back(iPoint); 
            */ 
        }

        if (!deathPeopleIdx_.empty()) std::cout << "-------- time step " << currentTimeStep << 
            "\t" << deathPeopleIdx_.size() << " people collapsed ---------" << std::endl;


        currentTimeStep++;

        // Regularization of velocity
        // Must have a norm lower to Vmax
        #pragma omp parallel for
        for(int iPoint = 0; iPoint < nbPoints_; iPoint++) {
            double vNorm = L2norm(velocityX_(currentTimeStep, iPoint), velocityY_(currentTimeStep, iPoint));
            double vMax =  args_.maxSpeed->coeff(iPoint);
            if (vNorm > vMax){
                velocityX_(currentTimeStep, iPoint) *= vMax / vNorm;
                velocityY_(currentTimeStep, iPoint) *= vMax / vNorm;
            }
        }

        
        positionX_.row(currentTimeStep) = positionX_.row(currentTimeStep - 1) + 
                velocityX_.row(currentTimeStep) * timeStep_;
        positionY_.row(currentTimeStep) = positionY_.row(currentTimeStep - 1) + 
                velocityY_.row(currentTimeStep) * timeStep_;


        // check if the walll has been crossed and regularize if it is the case
        // application of Ping pong rule
        std::vector<int> reIntegratedPoints {crossedWallPointsRegul(
                                                positionX_.row(currentTimeStep-1),
                                                positionY_.row(currentTimeStep-1),
                                                positionX_.row(currentTimeStep).data(),
                                                positionY_.row(currentTimeStep).data(),
                                                velocityX_.row(currentTimeStep-1),
                                                velocityY_.row(currentTimeStep-1),
                                                velocityX_.row(currentTimeStep).data(),
                                                velocityY_.row(currentTimeStep).data(),
                                                args_.Radius,
                                                wall_                        )};
        /*if(!reIntegratedPoints.empty()) {
            std::cout << "**** Reintegrated " << reIntegratedPoints.size() << 
                " Points after wall crossing **** " << std::endl;
        }*/

        // check arrival at doors
        #pragma omp parallel for shared(nbEjectedPoints_)
        for(int idx=0; idx < nbPoints_; idx++) {
            for(auto door : g_.doors_) {
                // check if the door is reached
                gridPoint currentPosition(positionX_(currentTimeStep, idx), positionY_(currentTimeStep, idx), 0);
                gridPoint pastPosition(positionX_(currentTimeStep - 1, idx), positionY_(currentTimeStep - 1, idx), 0);
                if (orientation(door.a_, door.b_, currentPosition) == 0 && onSegment(door.a_, currentPosition, door.b_) || 
                    doIntersect(door.a_, door.b_, pastPosition, currentPosition)) {
                        // release the point
                        std::uniform_real_distribution<double> distributionPX(loin-10000, loin+10000);
                        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                        std::default_random_engine generatorX (seed);
                        std::default_random_engine generatorY(seed+10);
                        positionX_(currentTimeStep, idx) = distributionPX(generatorX);
                        positionY_(currentTimeStep, idx) = distributionPX(generatorY);
                        #pragma omp critical
                        {
                            nbEjectedPoints_++;
                            if (nbEjectedPoints_ % 50 == 0) std::cout << "***** " << nbEjectedPoints_ << " points ejected **********" << std::endl;
                        }
                        break;
                }
            }
        }
        // compute the current normal direction to take, according to the grid distribution
        #pragma omp parallel for
        for(int i=0; i<nbPoints_; i++) {
            double x {positionX_(currentTimeStep, i)};
            double y {positionY_(currentTimeStep, i)};
            int kx = (int) std::floor((x-g_.xMin_) / g_.dx_);
            int ky = (int) std::floor((y-g_.yMin_) / g_.dy_);
            if (kx <= g_.iMax_ && ky <= g_.jMax_) {
                if (abs(x - (double) kx * g_.dx_) < epsi && abs(y - (double) ky * g_.dy_) < epsi) {
                    // the pedestrian is on a grid point
                    assert(kx != 0 && kx != g_.iMax_ && ky != 0 && ky != g_.jMax_);
                    normalX(currentTimeStep, i) = g_.get(kx, ky).gradX;
                    normalY(currentTimeStep, i) = g_.get(kx, ky).gradY;
                    double norm = L2norm(normalX(currentTimeStep, i), normalY(currentTimeStep, i));
                    while(norm == 0) {
                        // give a random direction to the point
                        std::uniform_real_distribution<double> direstionDist(-1.0, 1.0);
                        unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
                        std::default_random_engine generatorX (seed);
                        std::default_random_engine generatorY(seed+10);
                        normalX(currentTimeStep, i) = direstionDist(generatorX);
                        normalY(currentTimeStep, i) = direstionDist(generatorY);
                        norm = L2norm(normalX(currentTimeStep, i), normalY(currentTimeStep, i));
                    }
                }
                else if(abs(x - (double) kx * g_.dx_)< epsi) {
                    // the pedestrian is on a vertical edge
                    //assert(kx != 0 && kx != g_.iMax_ && ky + 1 <= g_.jMax_);
                    if(kx == g_.iMax_) kx--;
                    if(ky == g_.jMax_) ky--;
                    normalX(currentTimeStep, i) = 0.5 * (g_.get(kx, ky).gradX + g_.get(kx, ky+1).gradX);
                    normalY(currentTimeStep, i) = 0.5 * (g_.get(kx, ky).gradY + g_.get(kx, ky+1).gradY);
                    double norm = L2norm(normalX(currentTimeStep, i), normalY(currentTimeStep, i));
                    if(norm == 0) {
                        if(g_.get(kx, ky).T_ < g_.get(kx, ky+1).T_) {
                            normalX(currentTimeStep, i) = g_.get(kx, ky).gradX;
                            normalY(currentTimeStep, i) = g_.get(kx, ky).gradY;
                        }
                        else {
                            normalX(currentTimeStep, i) = g_.get(kx, ky+1).gradX;
                            normalY(currentTimeStep, i) = g_.get(kx, ky+1).gradY;
                        }
                    }
                }
                else if(abs(y - (double) ky * g_.dy_)< epsi) {
                    // the pedestrian is on a horizontal edge
                    //assert(ky != 0 && ky != g_.jMax_ && kx + 1 <= g_.iMax_);
                    if(kx == g_.iMax_) kx--;
                    if(ky == g_.jMax_) ky--;
                    normalX(currentTimeStep, i) = 0.5 * (g_.get(kx, ky).gradX + g_.get(kx + 1, ky).gradX);
                    normalY(currentTimeStep, i) = 0.5 * (g_.get(kx, ky).gradY + g_.get(kx + 1, ky).gradY);
                    double norm = L2norm(normalX(currentTimeStep, i), normalY(currentTimeStep, i));
                    if(norm == 0) {
                        if(g_.get(kx, ky).T_ < g_.get(kx+1, ky).T_) {
                            normalX(currentTimeStep, i) = g_.get(kx, ky).gradX;
                            normalY(currentTimeStep, i) = g_.get(kx, ky).gradY;
                        }
                        else{
                            normalX(currentTimeStep, i) = g_.get(kx+1, ky).gradX;
                            normalY(currentTimeStep, i) = g_.get(kx+1, ky).gradY;
                        }
                    }
                }
                else {
                    if(kx + 1 <= g_.iMax_ && ky + 1 <= g_.jMax_) { // inside the domain
                        normalX(currentTimeStep, i) = 0.25 * (g_.get(kx, ky).gradX + g_.get(kx + 1, ky).gradX + 
                                                            g_.get(kx, ky+1).gradX + g_.get(kx + 1, ky+1).gradX);
                        normalY(currentTimeStep, i) = 0.25 * (g_.get(kx, ky).gradY + g_.get(kx + 1, ky).gradY + 
                                                            g_.get(kx, ky+1).gradY + g_.get(kx + 1, ky+1).gradY);
                    }
                    else {
                        // point outside the domain
                        normalX(currentTimeStep, i) = 1.0;
                        normalY(currentTimeStep, i) = 1.0;
                    }
                    double norm = L2norm(normalX(currentTimeStep, i), normalY(currentTimeStep, i));
                    if(norm == 0) {
                        std::vector<double> dist {g_.get(kx, ky).T_, g_.get(kx+1, ky).T_, g_.get(kx, ky+1).T_, g_.get(kx+1, ky+1).T_};
                        int idx = std::distance(dist.begin(), std::min_element(dist.begin(), dist.end()));
                        switch (idx)
                        {
                            case 0:
                                normalX(currentTimeStep, i) = g_.get(kx, ky).gradX;
                                normalY(currentTimeStep, i) = g_.get(kx, ky).gradY;
                                break;
                            case 1:
                                normalX(currentTimeStep, i) = g_.get(kx+1, ky).gradX;
                                normalY(currentTimeStep, i) = g_.get(kx+1, ky).gradY;
                                break;
                            case 2:
                                normalX(currentTimeStep, i) = g_.get(kx, ky+1).gradX;
                                normalY(currentTimeStep, i) = g_.get(kx, ky+1).gradY;
                                break;
                            case 3:
                                normalX(currentTimeStep, i) = g_.get(kx+1, ky+1).gradX;
                                normalY(currentTimeStep, i) = g_.get(kx+1, ky+1).gradY;
                                break;
                            default:
                                break;
                        }
                    }
                }
            }
            else {
                    // point outside the domain
                    normalX(currentTimeStep, i) = 1.0;
                    normalY(currentTimeStep, i) = 1.0;
                }
        }
        #pragma omp parallel for
        for(int i=0; i<nbPoints_; i++){
            double norm = L2norm(normalX(currentTimeStep, i), normalY(currentTimeStep, i));
            if(norm == 0) {
                #pragma omp critical
                if(deb) std::cout << "need to adopt the nearest normal at x = " << 
                    positionX_(currentTimeStep, i) << " y = " << positionY_(currentTimeStep, i) << 
                    std::endl;
                std::vector<polygonalWall> *walls = &wall_.walls_;
                std::vector<std::vector<double>> wallVerticesX; wallVerticesX.reserve(wall_.nbLines());
                std::vector<std::vector<double>> wallVerticesY; wallVerticesY.reserve(wall_.nbLines());
                for(int idx = 0; idx < walls->size(); idx++) {
                    wallVerticesX.push_back((*walls)[idx].getVerticesX());
                    wallVerticesY.push_back((*walls)[idx].getVerticesY());
                }
                double x {positionX_(currentTimeStep, i)};
                double y {positionY_(currentTimeStep, i)};
                int kx = (int) std::floor((x-g_.xMin_) / g_.dx_);
                int ky = (int) std::floor((y-g_.yMin_) / g_.dy_);
                std::vector<std::pair<double, int>> tmp;
                std::set<int> notTakenIntoAccount;
                for(int jList=0; jList<wallVerticesX.size(); jList++) {
                    for (int idx=0; idx < wallVerticesX[jList].size() - 1; idx++) {
                        double xw1 = wallVerticesX[jList][idx];
                        double xw2 = wallVerticesX[jList][idx + 1];
                        double yw1 = wallVerticesY[jList][idx];
                        double yw2 = wallVerticesY[jList][idx + 1];
                        gridPoint wa(xw1, yw1, 0.0);
                        gridPoint wb(xw2, yw2, 0.0);
                        gridPoint pa(positionX_(currentTimeStep, i), positionY_(currentTimeStep, i), 0.0);
                        int kk = 0;
                        for(int i=-4; i<=4; i++){
                            for(int j=-4; j<=4; j++){
                                if(kx+i >= 0 && kx+i <= g_.iMax_ && ky+j>=0 && ky+j<=g_.jMax_){
                                    gridPoint pb(g_.get(kx+i, ky+j).x, 
                                                 g_.get(kx+i, ky+j).y, 0.0);
                                    if(doIntersect(wa, wb, pa, pb)){
                                        notTakenIntoAccount.insert(kk);
                                    }
                                }
                                kk++;
                            }
                        }
                    }
                }
                int kk=0;
                for(int i=-4; i<=4; i++){
                    for(int j=-4; j<=4; j++){
                        if(kx+i >= 0 && kx+i <= g_.iMax_ && ky+j>=0 && ky+j<=g_.jMax_){
                            if (std::find(notTakenIntoAccount.begin(), notTakenIntoAccount.end(), kk) == 
                                notTakenIntoAccount.end())
                                tmp.push_back(std::make_pair(g_.get(kx+i, ky+j).T_, 
                                                                g_.pointsMap_[std::make_pair(kx+i, ky+j)]));
                        }
                        kk++;
                    }
                }
                std::sort(tmp.begin(), tmp.end());

                

                int idx=0;
                while(isnan(norm) || norm == 0){
                    assert(idx < tmp.size());
                    normalX(currentTimeStep, i) = g_.get(tmp[idx].second).gradX;
                    normalY(currentTimeStep, i) = g_.get(tmp[idx++].second).gradY;
                    norm = L2norm(normalX(currentTimeStep, i), normalY(currentTimeStep, i));
                    if(isnan(norm) && deb) for(auto t : tmp)
                            std::cout << t.first << " " << g_.get(t.second).x << " " << 
                            g_.get(t.second).y << " " << g_.get(t.second).gradX << " " << 
                            g_.get(t.second).gradY << std::endl;
                }
            }
            assert(norm > 0);
            normalX(currentTimeStep,i) /= norm;
            normalY(currentTimeStep,i) /= norm;
        }

        #pragma omp parallel for
        for(int i=0; i<nbPoints_; i++) {
            for(int j=0; j<nbPoints_; j++) {
                if(i != j) {
                    rijX_(i,j) = positionX_(currentTimeStep, i) - positionX_(currentTimeStep, j);
                    rijY_(i,j) = positionY_(currentTimeStep, i) - positionY_(currentTimeStep, j);
                    dij_(i,j) = L2norm(rijX_(i,j), rijY_(i, j));
                    if (dij_(i,j) < epss) {
                        std::cout << "position regularization :: forcing two close pedestrians to step away" << std::endl;
                        //output("Xreprise.dat", "Yreprise.dat", "VXReprise.dat", "VYReprise.dat", "wall.dat");
                        //int cmp = 0;
                        #pragma omp critical
                        while(dij_(i, j) < epss) {
                            // avoid nan propagation when compression is huge
                            positionX_(currentTimeStep, i) += epss;
                            rijX_(i,j) = positionX_(currentTimeStep, i) - positionX_(currentTimeStep, j);
                            dij_(i, j) = L2norm(rijX_(i, j), rijY_(i, j));
                            //if(cmp++ < 20) std::cout << rijX_(i, j) << " " << dij_(i, j) << std::endl;
                        }
                    }
                    rijX_(i,j) /= dij_(i,j);
                    rijY_(i,j) /= dij_(i,j);
                }
            }
        }
        updateWallData();
    }
};