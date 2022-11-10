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

#ifndef QUADSHAPEFUNCTION_H
#define QUADSHAPEFUNCTION_H

// using namespace boost::numeric;

#include "IsogeometricParameters.h"

#include<cstdlib>
#include<fstream>
#include<iostream>
#include <math.h>       /* atan */
#include <vector>


// Defines the quadratic shape functions and its derivatives

template<int DIM>

class QuadShapeFunction {
public:

    // Defines the class Isogeometric Parameters locally
    typedef IsogeometricParameters<DIM>         IParameters_;
 
public:
    
    // Evaluates the shape function value
    void evaluateFem(double *xsi, double *phi) const;
    void evaluateIso(double *xsi, double *phi, double* wpcs, int *inc, std::vector<IParameters_ *> &iparameters ,int Npatch);

    // Evaluates the values of the shape funtion derivatives
    void evaluateGradientFem(double *xsi, double **dphi) const;   
    void evaluateGradientIso(double *xsi, double **dphi, double* wpcs, int *inc, std::vector<IParameters_ *> &iparameters ,int Npatch);

    void evaluateHessianFem(double ***ddphi);
    void evaluateHessianIso(double *xsi, double ***ddphi, double* wpcs, int *inc, std::vector<IParameters_ *> &iparameters ,int Npatch);

};


#endif
