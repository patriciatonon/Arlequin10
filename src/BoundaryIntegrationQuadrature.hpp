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

#ifndef BOUND_INTEG_QUADRATURE_H
#define BOUND_INTEG_QUADRATURE_H

#include "BoundaryShapeFunction.hpp"

// Defines the quadrature rule for the boundary integration
template<int DIM>

class BoundaryIntegQuadrature{

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
    BoundaryIntegQuadrature(){
        setQuadrature();
    }

private:
    ///List of integration points coordinates
    double pointCoord[5*DIM-8][DIM-1];

    ///List of integration points weights
    double pointWeight[5*DIM-8];
};

 //------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------   

template<>
void BoundaryIntegQuadrature<2>::setQuadrature(){

    int numPoints = 2;
    double xmga, xlga, zga, p1ga, p2ga, p3ga, ppga, z1ga;
    int mga;
    int nga = numPoints;

    mga=(nga+1.)/2.;
    xmga=0.0;
    xlga=1.0;
    
    for (int iga=1; iga<=mga; iga++) {
        zga = std::cos(pi*(double(iga)-0.25)/(double(nga)+0.5));
    g1:
        p1ga = 1.0;
        p2ga = 0.0;
        for (int jga=1; jga <= nga; jga++) {
            p3ga = p2ga;
            p2ga = p1ga;
            p1ga = ((2.0*double(jga)-1.0)*zga*p2ga-(double(jga)-1.0)*p3ga)/(double(jga));
        };
     
        ppga = nga*(zga*p1ga-p2ga)/(zga*zga-1.0);
        z1ga = zga;
        zga = z1ga-p1ga/ppga;
        
        if (std::fabs(zga-z1ga) > 1.0e-15) goto g1;
        
        pointCoord[iga-1][0] = xmga-xlga*zga;
        pointCoord[nga-iga][0] = xmga + xlga*zga;
        pointWeight[iga-1] = 2.0*xlga/((1.0-zga*zga)*ppga*ppga);
        pointWeight[nga-iga] = pointWeight[iga-1];
    };
    return;

}; 

    // Returns the integration point coordinate to a surface boundary 
template<>
void BoundaryIntegQuadrature<3>::setQuadrature(){

    pointCoord[0][0] = 1. / 3.;
    pointCoord[0][1] = 1. / 3.;
    
    pointCoord[1][0] = (9. + 2. * sqrt(15.)) / 21.;
    pointCoord[1][1] = (6. - sqrt(15.)) / 21.;
  
    pointCoord[2][0] = (6. - sqrt(15.)) / 21.;
    pointCoord[2][1] = (9. + 2. * sqrt(15.)) / 21.;
  
    pointCoord[3][0] = (6. - sqrt(15.)) / 21.;
    pointCoord[3][1] = (6. - sqrt(15.)) / 21.;
  
    pointCoord[4][0] = (6. + sqrt(15.)) / 21.;
    pointCoord[4][1] = (6. + sqrt(15.)) / 21.;
  
    pointCoord[5][0] = (9. - 2. * sqrt(15.)) / 21.;
    pointCoord[5][1] = (6. + sqrt(15.)) / 21.;
  
    pointCoord[6][0] = (6. + sqrt(15.)) / 21.;
    pointCoord[6][1] = (9. - 2. * sqrt(15.)) / 21.;

    pointWeight[0] = 0.11250;
    pointWeight[1] = (155. - sqrt(15.)) / 2400.;
    pointWeight[2] = (155. - sqrt(15.)) / 2400.;
    pointWeight[3] = (155. - sqrt(15.)) / 2400.;
    pointWeight[4] = (155. + sqrt(15.)) / 2400.;
    pointWeight[5] = (155. + sqrt(15.)) / 2400.;
    pointWeight[6] = (155. + sqrt(15.)) / 2400.; 

    return;
};


#endif
