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

#include "QuadraticShapeFunction.hpp"

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

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------QUADRATURE POINTS - COORDINATES------------------------
//------------------------------------------------------------------------------
template<>
double SpecialIntegQuadrature<2>::PointListFem(int i, int j){


    // pointCoordSpecialFem[0][0] = 0.873821971016996;
    // pointCoordSpecialFem[0][1] = 0.063089014491502;

    // pointCoordSpecialFem[1][0] = 0.063089014491502;
    // pointCoordSpecialFem[1][1] = 0.873821971016996;
      
    // pointCoordSpecialFem[2][0] = 0.063089014491502;
    // pointCoordSpecialFem[2][1] = 0.063089014491502;
      
    // pointCoordSpecialFem[3][0] = 0.501426509658179;
    // pointCoordSpecialFem[3][1] = 0.249286745170910;
      
    // pointCoordSpecialFem[4][0] = 0.249286745170910;
    // pointCoordSpecialFem[4][1] = 0.501426509658179;
      
    // pointCoordSpecialFem[5][0] = 0.249286745170910;
    // pointCoordSpecialFem[5][1] = 0.249286745170910;
      
    // pointCoordSpecialFem[6][0] = 0.636502499121399;
    // pointCoordSpecialFem[6][1] = 0.310352451033785;

    // pointCoordSpecialFem[7][0] = 0.636502499121399;
    // pointCoordSpecialFem[7][1] = 0.053145049844816;

    // pointCoordSpecialFem[8][0] = 0.310352451033785;
    // pointCoordSpecialFem[8][1] = 0.636502499121399;

    // pointCoordSpecialFem[9][0] = 0.310352451033785;
    // pointCoordSpecialFem[9][1] = 0.053145049844816;

    // pointCoordSpecialFem[10][0] = 0.053145049844816;
    // pointCoordSpecialFem[10][1] = 0.636502499121399;

    // pointCoordSpecialFem[11][0] = 0.053145049844816;
    // pointCoordSpecialFem[11][1] = 0.310352451033785;
    pointCoordSpecialFem[0][0] = 1. / 3.;
    pointCoordSpecialFem[0][1] = 1. / 3.;
        
    pointCoordSpecialFem[1][0] = (9. + 2. * sqrt(15.)) / 21.;
    pointCoordSpecialFem[1][1] = (6. - sqrt(15.)) / 21.;
      
    pointCoordSpecialFem[2][0]= (6. - sqrt(15.)) / 21.;
    pointCoordSpecialFem[2][1] = (9. + 2. * sqrt(15.)) / 21.;
      
    pointCoordSpecialFem[3][0] = (6. - sqrt(15.)) / 21.;
    pointCoordSpecialFem[3][1] = (6. - sqrt(15.)) / 21.;
      
    pointCoordSpecialFem[4][0] = (6. + sqrt(15.)) / 21.;
    pointCoordSpecialFem[4][1] = (6. + sqrt(15.)) / 21.;
      
    pointCoordSpecialFem[5][0] = (9. - 2. * sqrt(15.)) / 21.;
    pointCoordSpecialFem[5][1] = (6. + sqrt(15.)) / 21.;
      
    pointCoordSpecialFem[6][0] = (6. + sqrt(15.)) / 21.;
    pointCoordSpecialFem[6][1] = (9. - 2. * sqrt(15.)) / 21.;

    return pointCoordSpecialFem[i][j];
};


template<>
double SpecialIntegQuadrature<2>::PointListIso(int i, int j){

    pointCoordSpecialIso[0][0] = 0.;
    pointCoordSpecialIso[0][1] = 0.;
        
    pointCoordSpecialIso[1][0] = 0.774596669241483;
    pointCoordSpecialIso[1][1] = 0.000000000000000;
      
    pointCoordSpecialIso[2][0] = -0.774596669241483;
    pointCoordSpecialIso[2][1] = 0.;
      
    pointCoordSpecialIso[3][0] = 0.;
    pointCoordSpecialIso[3][1] = 0.774596669241483;
     
    pointCoordSpecialIso[4][0] =  0.774596669241483;
    pointCoordSpecialIso[4][1] =  0.774596669241483;
      
    pointCoordSpecialIso[5][0] = -0.774596669241483;
    pointCoordSpecialIso[5][1] =  0.774596669241483;

    pointCoordSpecialIso[6][0] =  0.;
    pointCoordSpecialIso[6][1] = -0.774596669241483;

    pointCoordSpecialIso[7][0] =  0.774596669241483;
    pointCoordSpecialIso[7][1] = -0.774596669241483;
      
    pointCoordSpecialIso[8][0] = -0.774596669241483;
    pointCoordSpecialIso[8][1] = -0.774596669241483;
        
    // pointCoordSpecial[0][0] = 0.0000000000;  pointCoordSpecial[0][1] = 0.0000000000;
    // pointCoordSpecial[1][0] = 0.0000000000;  pointCoordSpecial[1][1] = 0.5384693101;
    // pointCoordSpecial[2][0] = 0.0000000000;  pointCoordSpecial[2][1] = -0.5384693101;
    // pointCoordSpecial[3][0] = 0.0000000000;  pointCoordSpecial[3][1] = 0.9061798459;
    // pointCoordSpecial[4][0] = 0.0000000000;  pointCoordSpecial[4][1] = -0.9061798459;
    // pointCoordSpecial[5][0] = 0.5384693101;  pointCoordSpecial[5][1] = 0.0000000000;
    // pointCoordSpecial[6][0] = 0.5384693101;  pointCoordSpecial[6][1] = 0.5384693101;
    // pointCoordSpecial[7][0] = 0.5384693101;  pointCoordSpecial[7][1] = -0.5384693101;
    // pointCoordSpecial[8][0] = 0.5384693101;  pointCoordSpecial[8][1] = 0.9061798459;
    // pointCoordSpecial[9][0] = 0.5384693101;  pointCoordSpecial[9][1] = -0.9061798459;
    // pointCoordSpecial[10][0] = -0.5384693101; pointCoordSpecial[10][1] = 0.0000000000;
    // pointCoordSpecial[11][0] = -0.5384693101; pointCoordSpecial[11][1] = 0.5384693101;
    // pointCoordSpecial[12][0] = -0.5384693101; pointCoordSpecial[12][1] = -0.5384693101;
    // pointCoordSpecial[13][0] = -0.5384693101; pointCoordSpecial[13][1] = 0.9061798459;
    // pointCoordSpecial[14][0] = -0.5384693101; pointCoordSpecial[14][1] = -0.9061798459;
    // pointCoordSpecial[15][0] = 0.9061798459; pointCoordSpecial[15][1] = 0.0000000000;
    // pointCoordSpecial[16][0] = 0.9061798459; pointCoordSpecial[16][1] = 0.5384693101;
    // pointCoordSpecial[17][0] = 0.9061798459; pointCoordSpecial[17][1] = -0.5384693101;
    // pointCoordSpecial[18][0] = 0.9061798459; pointCoordSpecial[18][1] = 0.9061798459;
    // pointCoordSpecial[19][0] = 0.9061798459; pointCoordSpecial[19][1] = -0.9061798459;
    // pointCoordSpecial[20][0] = -0.9061798459; pointCoordSpecial[20][1] = 0.0000000000;
    // pointCoordSpecial[21][0] = -0.9061798459; pointCoordSpecial[21][1] = 0.5384693101;
    // pointCoordSpecial[22][0] = -0.9061798459; pointCoordSpecial[22][1] = -0.5384693101;
    // pointCoordSpecial[23][0] = -0.9061798459; pointCoordSpecial[23][1] = 0.9061798459;
    // pointCoordSpecial[24][0] = -0.9061798459; pointCoordSpecial[24][1] = -0.9061798459;
   
    return pointCoordSpecialIso[i][j];
};


//------------------------------------------------------------------------------
//-------------------------QUADRATURE POINTS - WEIGHTS--------------------------
//------------------------------------------------------------------------------
// template<>
// double SpecialIntegQuadrature<2>::WeightList(int i){
    
// 	pointWeightSpecial[0] = 0.010247216561;
// 	pointWeightSpecial[1] = 0.022511306623;
// 	pointWeightSpecial[2] = 0.031756064592;
// 	pointWeightSpecial[3] = 0.036713948533;
// 	pointWeightSpecial[4] = 0.036713948533;
// 	pointWeightSpecial[5] = 0.031756064592;
// 	pointWeightSpecial[6] = 0.022511306623;
// 	pointWeightSpecial[7] = 0.010247216561;
// 	pointWeightSpecial[8] = 0.022511306623;
// 	pointWeightSpecial[9] = 0.049453324505;
// 	pointWeightSpecial[10] = 0.069762408445;
// 	pointWeightSpecial[11] = 0.080653994949;
// 	pointWeightSpecial[12] = 0.080653994949;
// 	pointWeightSpecial[13] = 0.069762408445;
// 	pointWeightSpecial[14] = 0.049453324505;
// 	pointWeightSpecial[15] = 0.022511306623;
// 	pointWeightSpecial[16] = 0.031756064592;
// 	pointWeightSpecial[17] = 0.069762408445;
// 	pointWeightSpecial[18] = 0.098411859682;
// 	pointWeightSpecial[19] = 0.113776313213;
// 	pointWeightSpecial[20] = 0.113776313213;
// 	pointWeightSpecial[21] = 0.098411859682;
// 	pointWeightSpecial[22] = 0.069762408445;
// 	pointWeightSpecial[23] = 0.031756064592;
// 	pointWeightSpecial[24] = 0.036713948533;
// 	pointWeightSpecial[25] = 0.080653994949;
// 	pointWeightSpecial[26] = 0.113776313213;
// 	pointWeightSpecial[27] = 0.131539526741;
// 	pointWeightSpecial[28] = 0.131539526741;
// 	pointWeightSpecial[29] = 0.113776313213;
// 	pointWeightSpecial[30] = 0.080653994949;
// 	pointWeightSpecial[31] = 0.036713948533;
// 	pointWeightSpecial[32] = 0.036713948533;
// 	pointWeightSpecial[33] = 0.080653994949;
// 	pointWeightSpecial[34] = 0.113776313213;
// 	pointWeightSpecial[35] = 0.131539526741;
// 	pointWeightSpecial[36] = 0.131539526741;
// 	pointWeightSpecial[37] = 0.113776313213;
// 	pointWeightSpecial[38] = 0.080653994949;
// 	pointWeightSpecial[39] = 0.036713948533;
// 	pointWeightSpecial[40] = 0.031756064592;
// 	pointWeightSpecial[41] = 0.069762408445;
// 	pointWeightSpecial[42] = 0.098411859682;
// 	pointWeightSpecial[43] = 0.113776313213;
// 	pointWeightSpecial[44] = 0.113776313213;
// 	pointWeightSpecial[45] = 0.098411859682;
// 	pointWeightSpecial[46] = 0.069762408445;
// 	pointWeightSpecial[47] = 0.031756064592;
// 	pointWeightSpecial[48] = 0.022511306623;
// 	pointWeightSpecial[49] = 0.049453324505;
// 	pointWeightSpecial[50] = 0.069762408445;
// 	pointWeightSpecial[51] = 0.080653994949;
// 	pointWeightSpecial[52] = 0.080653994949;
// 	pointWeightSpecial[53] = 0.069762408445;
// 	pointWeightSpecial[54] = 0.049453324505;
// 	pointWeightSpecial[55] = 0.022511306623;
// 	pointWeightSpecial[56] = 0.010247216561;
// 	pointWeightSpecial[57] = 0.022511306623;
// 	pointWeightSpecial[58] = 0.031756064592;
// 	pointWeightSpecial[59] = 0.036713948533;
// 	pointWeightSpecial[60] = 0.036713948533;
// 	pointWeightSpecial[61] = 0.031756064592;
// 	pointWeightSpecial[62] = 0.022511306623;
// 	pointWeightSpecial[63] = 0.010247216561;

//     return pointWeightSpecial[i];
// };

template<>
double SpecialIntegQuadrature<2>::WeightListFem(int i){

    // pointWeightSpecialFem[0] = 0.050844906370207;
    // pointWeightSpecialFem[1] = 0.050844906370207;
    // pointWeightSpecialFem[2] = 0.050844906370207;
    // pointWeightSpecialFem[3] = 0.116786275726379;
    // pointWeightSpecialFem[4] = 0.116786275726379;
    // pointWeightSpecialFem[5] = 0.116786275726379;
    // pointWeightSpecialFem[6] = 0.082851075618374;
    // pointWeightSpecialFem[7] = 0.082851075618374;
    // pointWeightSpecialFem[8] = 0.082851075618374;
    // pointWeightSpecialFem[9] = 0.082851075618374;
    // pointWeightSpecialFem[10] = 0.082851075618374;
    // pointWeightSpecialFem[11] = 0.082851075618374; 
    pointWeightSpecialFem[0] = 0.11250;
    pointWeightSpecialFem[1] = (155. - sqrt(15.)) / 2400.;
    pointWeightSpecialFem[2] = (155. - sqrt(15.)) / 2400.;
    pointWeightSpecialFem[3] = (155. - sqrt(15.)) / 2400.;
    pointWeightSpecialFem[4] = (155. + sqrt(15.)) / 2400.;
    pointWeightSpecialFem[5] = (155. + sqrt(15.)) / 2400.;
    pointWeightSpecialFem[6] = (155. + sqrt(15.)) / 2400.; 

    return pointWeightSpecialFem[i];

};

template<>
double SpecialIntegQuadrature<2>::WeightListIso(int i){

     pointWeightSpecialIso[0] = 0.790123456790124;
    pointWeightSpecialIso[1] = 0.493827160493828;
    pointWeightSpecialIso[2] = 0.493827160493828;
    pointWeightSpecialIso[3] = 0.493827160493828;
    pointWeightSpecialIso[4] = 0.308641975308642;
    pointWeightSpecialIso[5] = 0.308641975308642;
    pointWeightSpecialIso[6] = 0.493827160493828;
    pointWeightSpecialIso[7] = 0.308641975308642;
    pointWeightSpecialIso[8] = 0.308641975308642;

    // pointWeightSpecial[0] = 0.3236345679;
    // pointWeightSpecial[1] = 0.2722865326;
    // pointWeightSpecial[2] = 0.2722865326;
    // pointWeightSpecial[3] = 0.1347850724;
    // pointWeightSpecial[4] = 0.1347850724;
    // pointWeightSpecial[5] = 0.2722865326;
    // pointWeightSpecial[6] = 0.2290854042;
    // pointWeightSpecial[7] = 0.2290854042;
    // pointWeightSpecial[8] = 0.1134000000;
    // pointWeightSpecial[9] = 0.1134000000;
    // pointWeightSpecial[10] = 0.2722865326;
    // pointWeightSpecial[11] = 0.2290854042;
    // pointWeightSpecial[12] = 0.2290854042;
    // pointWeightSpecial[13] = 0.1134000000;
    // pointWeightSpecial[14] = 0.1134000000;
    // pointWeightSpecial[15] = 0.1347850724;
    // pointWeightSpecial[16] = 0.1134000000;
    // pointWeightSpecial[17] = 0.1134000000;
    // pointWeightSpecial[18] = 0.0561343489;
    // pointWeightSpecial[19] = 0.0561343489;
    // pointWeightSpecial[20] = 0.1347850724;
    // pointWeightSpecial[21] = 0.1134000000;
    // pointWeightSpecial[22] = 0.1134000000;
    // pointWeightSpecial[23] = 0.0561343489;
    // pointWeightSpecial[24] = 0.0561343489;

    return pointWeightSpecialIso[i];
};

//------------------------------------------------------------------------------
//-----------COMPUTES THE VALUE INTERPOLATED IN THE INTEGRATION POINT-----------
//------------------------------------------------------------------------------
template<>
double SpecialIntegQuadrature<2>::interpolateQuadraticVariableFem(double *nValues, int &point) {
    
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
double SpecialIntegQuadrature<2>::interpolateQuadraticVariableIso(double *nValues, int &point, double *wpc, int *INC_, std::vector<IParameters *> &iparam, int &patch) {
    
    iparameters = &iparam;
    QuadShapeFunction<2>    shapeQuad;
    int dim = 2;

    double xsi[dim];
    for (int i = 0; i <dim ;i++) xsi[i] = PointListIso(point,i);

    double phi_[9];
    
    shapeQuad.evaluateIso(xsi,phi_,wpc,INC_,(*iparameters),patch);
    
    double int_value = 0.;
    for (int i = 0; i < 9; i++){
        int_value += nValues[i] * phi_[i];
    };

    return int_value;
};




#endif
