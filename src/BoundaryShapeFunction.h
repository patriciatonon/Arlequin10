//------------------------------------------------------------------------------
// 
//        Jeferson W D Fernandes,Rodolfo A K Sanches and Patrícia Tonon
//                             University of Sao Paulo
//                           (C) 2019 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------QUADRATIC SHAPE FUNCTION----------------------------
//------------------------------------------------------------------------------

#ifndef BOUND_SHAPEFUNCTION_H
#define BOUND_SHAPEFUNCTION_H

// using namespace boost::numeric;

#include "IsogeometricParameters.h"

#include<cstdlib>
#include<fstream>
#include<iostream>
#include <math.h>       /* atan */
#include <vector>



/// Defines the fluid boundary shape functions
template<int DIM>
class BoundShapeFunction {
public:
   
    // Defines the class Isogeometric Parameters locally
    typedef IsogeometricParameters<DIM>                              IParameters_;
  

public:

void evaluateBoundaryFem(double *Xsi, double *phiB_);

void evaluateBoundaryIso(double *Xsi, double *phiB_, double *wpcs, int side, int *inc, std::vector<IParameters_ *> &iparameters, int Npatch);

void evaluateGradientBoundaryFem(double *Xsi, double **dphiB_);

void evaluateGradientBoundaryIso(double *Xsi, double **dphiB_, double *wpcs, int side, int *inc, std::vector<IParameters_ *> &iparameters, int Npatch);

void evaluateHessianBoundaryFem(double ***ddphiB_);
};


#endif
