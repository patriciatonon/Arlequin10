//------------------------------------------------------------------------------
// 
//        Jeferson W D Fernandes,Rodolfo A K Sanches and Patrícia Tonon
//                             University of Sao Paulo
//                           (C) 2019 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------QUADRATURE POINTS-------------------------------
//------------------------------------------------------------------------------

#ifndef BOUND_INTEG_QUADRATURE_ISO_H
#define BOUND_INTEG_QUADRATURE_ISO_H

#include "BoundaryShapeFunction.h"

#include<cstdlib>
#include<fstream>
#include<iostream>
#include <math.h>       /* atan */

// Defines the quadrature rule for the boundary integration
template<int DIM>

class BoundaryIntegQuadratureIso{
private:
    double pi = M_PI;  //Pi

public:
    // Returns the index of the first integration point
    double* begin() {
        return std::begin(pointWeight);
    }

    // Returns the index of the last integration point
    double* end() {
        return std::end(pointWeight);
    }

    // Returns the integration point coordinate
    double PointList(int i, int j) {return pointCoord[i][j];}; 

    // Retuns the integration point weight
    double WeightList(int i) {return pointWeight[i];};

    // Returns the integration point adimensional coordinate and weight
    void setQuadrature();

    /// Constructor of the domain integration quadrature
    BoundaryIntegQuadratureIso(){
        setQuadrature();
    }

private:
    ///List of integration points coordinates
    double pointCoord[7*DIM-12][DIM-1];

    ///List of integration points weights
    double pointWeight[7*DIM-12];
};


#endif
