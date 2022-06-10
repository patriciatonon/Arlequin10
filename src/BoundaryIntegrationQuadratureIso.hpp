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

#include "BoundaryShapeFunction.hpp"

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

template<>
void BoundaryIntegQuadratureIso<2>::setQuadrature(){

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

template<>
void BoundaryIntegQuadratureIso<3>::setQuadrature(){

    pointCoord[0][0] = 0.;
    pointCoord[0][1] = 0.;
        
    pointCoord[1][0] = 0.774596669241483;
    pointCoord[1][1] = 0.000000000000000;
      
    pointCoord[2][0] = -0.774596669241483;
    pointCoord[2][1] = 0.;
      
    pointCoord[3][0] = 0.;
    pointCoord[3][1] = 0.774596669241483;
     
    pointCoord[4][0] =  0.774596669241483;
    pointCoord[4][1] =  0.774596669241483;
      
    pointCoord[5][0] = -0.774596669241483;
    pointCoord[5][1] =  0.774596669241483;

    pointCoord[6][0] =  0.;
    pointCoord[6][1] = -0.774596669241483;

    pointCoord[7][0] =  0.774596669241483;
    pointCoord[7][1] = -0.774596669241483;
      
    pointCoord[8][0] = -0.774596669241483;
    pointCoord[8][1] = -0.774596669241483;
   
   	pointWeight[0] = 0.790123456790124;
    pointWeight[1] = 0.493827160493828;
    pointWeight[2] = 0.493827160493828;
    pointWeight[3] = 0.493827160493828;
    pointWeight[4] = 0.308641975308642;
    pointWeight[5] = 0.308641975308642;
    pointWeight[6] = 0.493827160493828;
    pointWeight[7] = 0.308641975308642;
    pointWeight[8] = 0.308641975308642;
    return;
};


#endif
