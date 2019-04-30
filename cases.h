#pragma once
#include "FMM/grid.h"
class Reunion{
    Reunion(){}
    ~Reunion(){}
public:
    static grid initializeGrid(Wall &w, double dx = 0.15, double dy = 0.1);
    static bool domainIsEmpty(volatile long nbEjectedPoints, long nbPoints) {
        return nbEjectedPoints == nbPoints;
    }
};

class Cinema{
    Cinema(){}
    ~Cinema(){}
public:
    static grid initializeGrid(Wall &w, double dx = 0.15, double dy = 0.1);
    //static grid initializeGrid(Wall &w, double dx = 0.5, double dy = 0.3);
    static bool domainIsEmpty(volatile long nbEjectedPoints, long nbPoints) {
        return nbEjectedPoints == nbPoints;
    }
};


class Zenith{
    Zenith(){}
    ~Zenith(){}
public:
    static grid initializeGrid(Wall &w, double dx = 0.15, double dy = 0.1);
    static bool domainIsEmpty(volatile long nbEjectedPoints, long nbPoints) {
        return nbEjectedPoints >= (nbPoints - 10);
    }
};