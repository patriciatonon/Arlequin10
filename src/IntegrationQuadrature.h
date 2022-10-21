//------------------------------------------------------------------------------
// 
//        Jeferson W D Fernandes,Rodolfo A K Sanches and Patrícia Tonon
//                             University of Sao Paulo
//                           (C) 2019 All Rights Reserved
//
// <LicenseText>
//----------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------QUADRATURE POINTS-------------------------------
//------------------------------------------------------------------------------

#ifndef INTEG_QUADRATURE_H
#define INTEG_QUADRATURE_H

#include "QuadraticShapeFunction.h"

#include<cstdlib>
#include<fstream>
#include<iostream>
#include <math.h>       /* atan */

// Defines the domain integration Hammer quadrature

template<int DIM>
class IntegQuadrature{

public:

    // Isogeometric parameters
    typedef IsogeometricParameters<DIM>                         IParameters;

    // Returns the index of the first integration point
    double* beginFem() {
        return std::begin(pointWeightFem);
    }

    double* beginIso() {
        return std::begin(pointWeightIso);
    }

    // Returns the index of the last integration point
    double* endFem() {
        return std::end(pointWeightFem);
    }

    double* endIso() {
        return std::end(pointWeightIso);
    }

    // Returns the integration point coordinate
    double PointListFem(int i, int j); 
    double PointListIso(int i, int j); 

    //Retuns the integration point weight
    double WeightListFem(int i);
    double WeightListIso(int i);
  
    double interpolateQuadraticVariableFem(double *nValues, int &point);
    double interpolateQuadraticVariableIso(double *nValues, int &point, double *wpc, int *INC_, std::vector<IParameters *> &iparam, int &patch);

private:
    std::vector<IParameters *> *iparameters;

    //List of integration points coordinates
    double pointCoordFem[8*DIM-9][DIM];
    double pointCoordIso[18*DIM-27][DIM];

    //List of integration points weights
    double pointWeightFem[8*DIM-9];
    double pointWeightIso[18*DIM-27];

};


#endif
