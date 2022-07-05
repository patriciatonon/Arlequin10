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

#include "QuadraticShapeFunction.hpp"

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

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------QUADRATURE POINTS - COORDINATES------------------------
//------------------------------------------------------------------------------
template<>
double IntegQuadrature<2>::PointListFem(int i, int j){

    pointCoordFem[0][0] = 1. / 3.;
    pointCoordFem[0][1] = 1. / 3.;
        
    pointCoordFem[1][0] = (9. + 2. * sqrt(15.)) / 21.;
    pointCoordFem[1][1] = (6. - sqrt(15.)) / 21.;
      
    pointCoordFem[2][0]= (6. - sqrt(15.)) / 21.;
    pointCoordFem[2][1] = (9. + 2. * sqrt(15.)) / 21.;
      
    pointCoordFem[3][0] = (6. - sqrt(15.)) / 21.;
    pointCoordFem[3][1] = (6. - sqrt(15.)) / 21.;
      
    pointCoordFem[4][0] = (6. + sqrt(15.)) / 21.;
    pointCoordFem[4][1] = (6. + sqrt(15.)) / 21.;
      
    pointCoordFem[5][0] = (9. - 2. * sqrt(15.)) / 21.;
    pointCoordFem[5][1] = (6. + sqrt(15.)) / 21.;
      
    pointCoordFem[6][0] = (6. + sqrt(15.)) / 21.;
    pointCoordFem[6][1] = (9. - 2. * sqrt(15.)) / 21.;

    return pointCoordFem[i][j];
};

template<>
double IntegQuadrature<2>::PointListIso(int i, int j){

    pointCoordIso[0][0] = 0.;
    pointCoordIso[0][1] = 0.;
        
    pointCoordIso[1][0] = 0.774596669241483;
    pointCoordIso[1][1] = 0.000000000000000;
      
    pointCoordIso[2][0] = -0.774596669241483;
    pointCoordIso[2][1] = 0.;
      
    pointCoordIso[3][0] = 0.;
    pointCoordIso[3][1] = 0.774596669241483;
     
    pointCoordIso[4][0] =  0.774596669241483;
    pointCoordIso[4][1] =  0.774596669241483;
      
    pointCoordIso[5][0] = -0.774596669241483;
    pointCoordIso[5][1] =  0.774596669241483;

    pointCoordIso[6][0] =  0.;
    pointCoordIso[6][1] = -0.774596669241483;

    pointCoordIso[7][0] =  0.774596669241483;
    pointCoordIso[7][1] = -0.774596669241483;
      
    pointCoordIso[8][0] = -0.774596669241483;
    pointCoordIso[8][1] = -0.774596669241483;
    
    return pointCoordIso[i][j];
};


// template<>
// double IntegQuadrature<3>::PointList(int i, int j){
    
//     pointCoord[0][0] = 1. / 4.;
//     pointCoord[0][1] = 1. / 4.;
//     pointCoord[0][2] = 1. / 4.;

//     pointCoord[1][0] = 0.;
//     pointCoord[1][1] = 1./3.;
//     pointCoord[1][2] = 1./3.;

//     pointCoord[2][0] = 1./3.;
//     pointCoord[2][1] = 1./3.;
//     pointCoord[2][2] = 1./3.;

//     pointCoord[3][0] = 1./3.;
//     pointCoord[3][1] = 1./3.;
//     pointCoord[3][2] = 0.;

//     pointCoord[4][0] = 1./3.;
//     pointCoord[4][1] = 0.;
//     pointCoord[4][2] = 1./3.;

//     pointCoord[5][0] = 8./11.;
//     pointCoord[5][1] = 1./11.;
//     pointCoord[5][2] = 1./11.;

//     pointCoord[6][0] = 1./11.;
//     pointCoord[6][1] = 1./11.;
//     pointCoord[6][2] = 1./11.;

//     pointCoord[7][0] = 1./11.;
//     pointCoord[7][1] = 1./11.;
//     pointCoord[7][2] = 8./11.;

//     pointCoord[8][0] = 1./11.;
//     pointCoord[8][1] = 8./11.;
//     pointCoord[8][2] = 1./11.;

//     pointCoord[9][0] = 0.4334498464263357;
//     pointCoord[9][1] = 0.0665501535736643;
//     pointCoord[9][2] = 0.0665501535736643;

//     pointCoord[10][0] = 0.0665501535736643;
//     pointCoord[10][1] = 0.4334498464263357;
//     pointCoord[10][2] = 0.0665501535736643;

//     pointCoord[11][0] = 0.0665501535736643;
//     pointCoord[11][1] = 0.0665501535736643;
//     pointCoord[11][2] = 0.4334498464263357;

//     pointCoord[12][0] = 0.0665501535736643;
//     pointCoord[12][1] = 0.4334498464263357;
//     pointCoord[12][2] = 0.4334498464263357;

//     pointCoord[13][0] = 0.4334498464263357;
//     pointCoord[13][1] = 0.0665501535736643;
//     pointCoord[13][2] = 0.4334498464263357;

//     pointCoord[14][0] = 0.4334498464263357;
//     pointCoord[14][1] = 0.4334498464263357;
//     pointCoord[14][2] = 0.0665501535736643;

//     return pointCoord[i][j];
// }

// template<>
// double IntegQuadrature<3>::PointListIso(int i, int j){
    
//     pointCoordIso[0][0] =   0.;
//     pointCoordIso[0][1] =   0.;
//     pointCoordIso[0][2] =   0.;   
//     pointCoordIso[1][0] =   0.774596669241483;
//     pointCoordIso[1][1] =   0.000000000000000;
//     pointCoordIso[1][2] =   0.;
//     pointCoordIso[2][0] =  -0.774596669241483;
//     pointCoordIso[2][1] =   0.;
//     pointCoordIso[2][2] =   0.; 
//     pointCoordIso[3][0] =   0.;
//     pointCoordIso[3][1] =   0.774596669241483;
//     pointCoordIso[3][2] =   0.;
//     pointCoordIso[4][0] =   0.774596669241483;
//     pointCoordIso[4][1] =   0.774596669241483;
//     pointCoordIso[4][2] =   0.;
//     pointCoordIso[5][0] =  -0.774596669241483;
//     pointCoordIso[5][1] =   0.774596669241483;
//     pointCoordIso[5][2] =   0.;
//     pointCoordIso[6][0] =   0.;
//     pointCoordIso[6][1] =  -0.774596669241483;
//     pointCoordIso[6][2] =   0.;
//     pointCoordIso[7][0] =   0.774596669241483;
//     pointCoordIso[7][1] =   -0.774596669241483;
//     pointCoordIso[7][2] =   0.;  
//     pointCoordIso[8][0] =  -0.774596669241483;
//     pointCoordIso[8][1] =  -0.774596669241483;
//     pointCoordIso[8][2] =   0.;
//     pointCoordIso[9][0] =   0.;
//     pointCoordIso[9][1] =   0.;
//     pointCoordIso[9][2] =   0.774596669241483;  
//     pointCoordIso[10][0] =  0.774596669241483;
//     pointCoordIso[10][1] =  0.000000000000000;
//     pointCoordIso[10][2] =  0.774596669241483;
//     pointCoordIso[11][0] = -0.774596669241483;
//     pointCoordIso[11][1] =  0.;
//     pointCoordIso[11][2] =  0.774596669241483;
//     pointCoordIso[12][0] =  0.;
//     pointCoordIso[12][1] =  0.774596669241483;
//     pointCoordIso[12][2] =  0.774596669241483;
//     pointCoordIso[13][0] =  0.774596669241483;
//     pointCoordIso[13][1] =  0.774596669241483;
//     pointCoordIso[13][2] =  0.774596669241483;
//     pointCoordIso[14][0] = -0.774596669241483;
//     pointCoordIso[14][1] =  0.774596669241483;
//     pointCoordIso[14][2] =  0.774596669241483;
//     pointCoordIso[15][0] =  0.;
//     pointCoordIso[15][1] = -0.774596669241483;
//     pointCoordIso[15][2] =  0.774596669241483;
//     pointCoordIso[16][0] =  0.774596669241483;
//     pointCoordIso[16][1] = -0.774596669241483;
//     pointCoordIso[16][2] =  0.774596669241483;
//     pointCoordIso[17][0] = -0.774596669241483;
//     pointCoordIso[17][1] = -0.774596669241483;
//     pointCoordIso[17][2] =  0.774596669241483;
//     pointCoordIso[18][0] =  0.;
//     pointCoordIso[18][1] =  0.;
//     pointCoordIso[18][2] = -0.774596669241483; 
//     pointCoordIso[19][0] =  0.774596669241483;
//     pointCoordIso[19][1] =  0.000000000000000;
//     pointCoordIso[19][2] = -0.774596669241483;
//     pointCoordIso[20][0] = -0.774596669241483;
//     pointCoordIso[20][1] =  0.;
//     pointCoordIso[20][2] = -0.774596669241483;
//     pointCoordIso[21][0] =  0.;
//     pointCoordIso[21][1] =  0.774596669241483;
//     pointCoordIso[21][2] = -0.774596669241483;
//     pointCoordIso[22][0] =  0.774596669241483;
//     pointCoordIso[22][1] =  0.774596669241483;
//     pointCoordIso[22][2] = -0.774596669241483;  
//     pointCoordIso[23][0] = -0.774596669241483;
//     pointCoordIso[23][1] =  0.774596669241483;
//     pointCoordIso[23][2] = -0.774596669241483;
//     pointCoordIso[24][0] =  0.;
//     pointCoordIso[24][1] = -0.774596669241483;
//     pointCoordIso[24][2] = -0.774596669241483;
//     pointCoordIso[25][0] =  0.774596669241483;
//     pointCoordIso[25][1] = -0.774596669241483;
//     pointCoordIso[25][2] = -0.774596669241483;    
//     pointCoordIso[26][0] = -0.774596669241483;
//     pointCoordIso[26][1] = -0.774596669241483;
//     pointCoordIso[26][2] = -0.774596669241483;

//     return pointCoordIso[i][j];
// }

//------------------------------------------------------------------------------
//-------------------------QUADRATURE POINTS - WEIGHTS--------------------------
//------------------------------------------------------------------------------
template<>
double IntegQuadrature<2>::WeightListFem(int i){
    
    pointWeightFem[0] = 0.11250;
    pointWeightFem[1] = (155. - sqrt(15.)) / 2400.;
    pointWeightFem[2] = (155. - sqrt(15.)) / 2400.;
    pointWeightFem[3] = (155. - sqrt(15.)) / 2400.;
    pointWeightFem[4] = (155. + sqrt(15.)) / 2400.;
    pointWeightFem[5] = (155. + sqrt(15.)) / 2400.;
    pointWeightFem[6] = (155. + sqrt(15.)) / 2400.; 

    return pointWeightFem[i];
};

template<>
double IntegQuadrature<2>::WeightListIso(int i){
    
    
    pointWeightIso[0] = 0.790123456790124;
    pointWeightIso[1] = 0.493827160493828;
    pointWeightIso[2] = 0.493827160493828;
    pointWeightIso[3] = 0.493827160493828;
    pointWeightIso[4] = 0.308641975308642;
    pointWeightIso[5] = 0.308641975308642;
    pointWeightIso[6] = 0.493827160493828;
    pointWeightIso[7] = 0.308641975308642;
    pointWeightIso[8] = 0.308641975308642;


    return pointWeightIso[i];
};


// template<>
// double IntegQuadrature<3>::WeightList(int i){
    
//     pointWeight[0] = 0.1817020685825351/6.;
//     pointWeight[1] = 0.0361607142857143/6.;
//     pointWeight[2] = 0.0361607142857143/6.;
//     pointWeight[3] = 0.0361607142857143/6.;
//     pointWeight[4] = 0.0361607142857143/6.;
//     pointWeight[5] = 0.0698714945161738/6.;
//     pointWeight[6] = 0.0698714945161738/6.; 
//     pointWeight[7] = 0.0698714945161738/6.;
//     pointWeight[8] = 0.0698714945161738/6.; 
//     pointWeight[9] = 0.0656948493683187/6.;
//     pointWeight[10] = 0.0656948493683187/6.; 
//     pointWeight[11] = 0.0656948493683187/6.;
//     pointWeight[12] = 0.0656948493683187/6.;
//     pointWeight[13] = 0.0656948493683187/6.;
//     pointWeight[14] = 0.0656948493683187/6.;

//     return pointWeight[i];
// };

// template<>
// double IntegQuadrature<3>::WeightListIso(int i){
    
//     pointWeightIso[0] = 0.702331961591221;
//     pointWeightIso[1] = 0.438957475994514;
//     pointWeightIso[2] = 0.438957475994514;
//     pointWeightIso[3] = 0.438957475994514;
//     pointWeightIso[4] = 0.274348422496571;
//     pointWeightIso[5] = 0.274348422496571;
//     pointWeightIso[6] = 0.438957475994514;
//     pointWeightIso[7] = 0.274348422496571;
//     pointWeightIso[8] = 0.274348422496571;
//     pointWeightIso[9] = 0.438957475994514;
//     pointWeightIso[10] = 0.274348422496571;
//     pointWeightIso[11] = 0.274348422496571;
//     pointWeightIso[12] = 0.274348422496571;
//     pointWeightIso[13] = 0.171467764060357;
//     pointWeightIso[14] = 0.171467764060357;
//     pointWeightIso[15] = 0.274348422496571;
//     pointWeightIso[16] = 0.171467764060357;
//     pointWeightIso[17] = 0.171467764060357;
//     pointWeightIso[18] = 0.438957475994514;
//     pointWeightIso[19] = 0.274348422496571;
//     pointWeightIso[20] = 0.274348422496571;
//     pointWeightIso[21] = 0.274348422496571;
//     pointWeightIso[22] = 0.171467764060357;
//     pointWeightIso[23] = 0.171467764060357;
//     pointWeightIso[24] = 0.274348422496571;
//     pointWeightIso[25] = 0.171467764060357;
//     pointWeightIso[26] = 0.171467764060357;

//     return pointWeightIso[i];
// };

//------------------------------------------------------------------------------
//-----------COMPUTES THE VALUE INTERPOLATED IN THE INTEGRATION POINT-----------
//------------------------------------------------------------------------------
template<>
double IntegQuadrature<2>::interpolateQuadraticVariableFem(double *nValues, int &point) {
    
    int dim = 2;
    QuadShapeFunction<2>    shapeQuad;

    double xsi[dim];

    for (int i = 0; i < dim; i++) xsi[i] = PointListFem(point,i);
    
    double phi_[6];
    shapeQuad.evaluateFem(xsi,phi_);
    
    double int_value = 0.;
    for (int i = 0; i < 6; i++){
        int_value += nValues[i] * phi_[i];
    };

    return int_value;
};


template<>
double IntegQuadrature<2>::interpolateQuadraticVariableIso(double *nValues, int &point, double *wpc, int *INC_, std::vector<IParameters *> &iparam, int &patch) {
    
    int dim = 2;
    iparameters = &iparam;
    QuadShapeFunction<2>    shapeQuad;

    double xsi[dim];
    for (int i = 0; i < dim; i++) xsi[i] = PointListIso(point,i);

    double phi_[9];
    shapeQuad.evaluateIso(xsi,phi_,wpc,INC_,(*iparameters),patch);
    
    double int_value = 0.;
    for (int i = 0; i < 9; i++){
        int_value += nValues[i] * phi_[i];
    };

    return int_value;
};

#endif
