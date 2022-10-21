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

#ifndef SPECIAL_INTEG_QUADRATURE_H
#define SPECIAL_INTEG_QUADRATURE_H

#include "QuadraticShapeFunction.h"

#include<cstdlib>
#include<fstream>
#include<iostream>
#include <math.h>       /* atan */

// Defines the domain integration Hammer quadrature

template<int DIM>
class SpecialIntegQuadrature{

public:

    // Isogeometric parameters
    typedef IsogeometricParameters<DIM>                         IParameters;

    // Returns the index of the first integration point
    double* beginIso(){
    	return std::begin(pointWeightSpecialIso);
    }

    // Returns the index of the last integration point
    double* endIso() {
        return std::end(pointWeightSpecialIso);
    }

      // Returns the index of the first integration point
    double* beginFem(){
        return std::begin(pointWeightSpecialFem);
    }

    // Returns the index of the last integration point
    double* endFem() {
        return std::end(pointWeightSpecialFem);
    }


    // Returns the integration point coordinate
    double PointListIso(int i, int j);

    //Retuns the integration point weight
    double WeightListIso(int i);

    // Returns the integration point coordinate
    double PointListFem(int i, int j);

    //Retuns the integration point weight
    double WeightListFem(int i);
    
    //Interpolates variables
    double interpolateQuadraticVariableIso( double *nValues, int &point, double *wpc, int *INC_, std::vector<IParameters *> &iparam, int &patch);
    double interpolateQuadraticVariableFem( double *nValues, int &point);

private:
    std::vector<IParameters *> *iparameters;

    //List of integration points coordinates
    // double pointCoordSpecial[64][DIM];
    // double pointCoordSpecial[25][DIM]; //Reduced
    double pointCoordSpecialIso[18*DIM-27][DIM];

    //List of integration points weights
    // double pointWeightSpecial[64];
    // double pointWeightSpecial[25]; //Reduced
    double pointWeightSpecialIso[18*DIM-27]; 

    double pointCoordSpecialFem[8*DIM-9][DIM];
    double pointWeightSpecialFem[8*DIM-9]; //Reduced

    // double pointCoordSpecialFem[12][DIM];
    // double pointWeightSpecialFem[12]; 
};

#endif
