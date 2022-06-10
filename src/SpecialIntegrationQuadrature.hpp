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
    double* begin(){
    	return std::begin(pointWeightSpecial);
    }

    // Returns the index of the last integration point
    double* end() {
        return std::end(pointWeightSpecial);
    }

    // Returns the integration point coordinate
    double PointList(int i, int j);

    //Retuns the integration point weight
    double WeightList(int i);
  
    double interpolateQuadraticVariable( double *nValues, int &point, double *wpc, int *INC_, std::vector<IParameters *> &iparam, int &patch);

private:
    std::vector<IParameters *> *iparameters;

    //List of integration points coordinates
    // double pointCoordSpecial[64][DIM];
    // double pointCoordSpecial[25][DIM]; //Reduced
    double pointCoordSpecial[9][DIM];

    //List of integration points weights
    // double pointWeightSpecial[64];
    // double pointWeightSpecial[25]; //Reduced
    double pointWeightSpecial[9]; //Reduced

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------QUADRATURE POINTS - COORDINATES------------------------
//------------------------------------------------------------------------------
// template<>
// double SpecialIntegQuadrature<2>::PointList(int i, int j){

        
//   	pointCoordSpecial[0][0] = 0.96029856500;   pointCoordSpecial[0][1] = 0.960298565000;
//     pointCoordSpecial[1][0] = 0.96029856500;   pointCoordSpecial[1][1] = 0.796666477400;
// 	pointCoordSpecial[2][0] = 0.96029856500;   pointCoordSpecial[2][1] = 0.525532409900;
// 	pointCoordSpecial[3][0] = 0.96029856500;   pointCoordSpecial[3][1] = 0.183434642500;
//     pointCoordSpecial[4][0] = 0.96029856500;   pointCoordSpecial[4][1] =-0.183434642500;
//     pointCoordSpecial[5][0] = 0.96029856500;   pointCoordSpecial[5][1] =	-0.525532409900;
//     pointCoordSpecial[6][0] = 0.96029856500;   pointCoordSpecial[6][1] =	-0.796666477400;
//     pointCoordSpecial[7][0] = 0.96029856500;   pointCoordSpecial[7][1] =	-0.960298565000;
//     pointCoordSpecial[8][0] = 0.79666647740;   pointCoordSpecial[8][1] =	0.960298565000;
//     pointCoordSpecial[9][0] = 0.79666647740;   pointCoordSpecial[9][1] = 0.796666477400;
//     pointCoordSpecial[10][0] = 0.79666647740;  pointCoordSpecial[10][1] = 0.525532409900;
//     pointCoordSpecial[11][0] = 0.79666647740;  pointCoordSpecial[11][1] = 0.183434642500;
//     pointCoordSpecial[12][0] = 0.79666647740;  pointCoordSpecial[12][1] =-0.183434642500;
//     pointCoordSpecial[13][0] = 0.79666647740;  pointCoordSpecial[13][1] =-0.525532409900;
//     pointCoordSpecial[14][0] = 0.79666647740;  pointCoordSpecial[14][1] = -0.796666477400;
// 	pointCoordSpecial[15][0] = 0.79666647740;  pointCoordSpecial[15][1] =-0.960298565000;
// 	pointCoordSpecial[16][0] = 0.52553240990;  pointCoordSpecial[16][1] = 0.960298565000;
// 	pointCoordSpecial[17][0] = 0.52553240990;  pointCoordSpecial[17][1] = 0.796666477400;
// 	pointCoordSpecial[18][0] = 0.52553240990;  pointCoordSpecial[18][1] = 0.525532409900;
// 	pointCoordSpecial[19][0] = 0.52553240990;  pointCoordSpecial[19][1] = 0.183434642500;
// 	pointCoordSpecial[20][0] = 0.52553240990;  pointCoordSpecial[20][1] = -0.183434642500;
// 	pointCoordSpecial[21][0] = 0.52553240990;  pointCoordSpecial[21][1] = -0.525532409900;
// 	pointCoordSpecial[22][0] = 0.52553240990;  pointCoordSpecial[22][1] = -0.796666477400;
// 	pointCoordSpecial[23][0] = 0.52553240990;  pointCoordSpecial[23][1] = -0.960298565000;
// 	pointCoordSpecial[24][0] = 0.18343464250;  pointCoordSpecial[24][1] = 0.960298565000;
// 	pointCoordSpecial[25][0] = 0.18343464250;  pointCoordSpecial[25][1] = 0.796666477400;
// 	pointCoordSpecial[26][0] = 0.18343464250;  pointCoordSpecial[26][1] = 0.525532409900;
// 	pointCoordSpecial[27][0] = 0.18343464250;  pointCoordSpecial[27][1] = 0.183434642500;
// 	pointCoordSpecial[28][0] = 0.18343464250;  pointCoordSpecial[28][1] = -0.183434642500;
// 	pointCoordSpecial[29][0] = 0.18343464250;  pointCoordSpecial[29][1] = -0.525532409900;
// 	pointCoordSpecial[30][0] = 0.18343464250;  pointCoordSpecial[30][1] = -0.796666477400;
// 	pointCoordSpecial[31][0] = 0.18343464250;  pointCoordSpecial[31][1] = -0.960298565000;
// 	pointCoordSpecial[32][0] = -0.18343464250; pointCoordSpecial[32][1] = 0.960298565000;
// 	pointCoordSpecial[33][0] = -0.18343464250; pointCoordSpecial[33][1] = 0.796666477400;
// 	pointCoordSpecial[34][0] = -0.18343464250; pointCoordSpecial[34][1] = 0.525532409900;
// 	pointCoordSpecial[35][0] = -0.18343464250; pointCoordSpecial[35][1] = 0.183434642500;
// 	pointCoordSpecial[36][0] = -0.18343464250; pointCoordSpecial[36][1] = -0.183434642500;
// 	pointCoordSpecial[37][0] = -0.18343464250; pointCoordSpecial[37][1] = -0.525532409900;
// 	pointCoordSpecial[38][0] = -0.18343464250; pointCoordSpecial[38][1] = -0.796666477400;
// 	pointCoordSpecial[39][0] = -0.18343464250; pointCoordSpecial[39][1] = -0.960298565000;
// 	pointCoordSpecial[40][0] = -0.52553240990; pointCoordSpecial[40][1] = 0.960298565000;
// 	pointCoordSpecial[41][0] = -0.52553240990; pointCoordSpecial[41][1] = 0.796666477400;
// 	pointCoordSpecial[42][0] = -0.52553240990; pointCoordSpecial[42][1] = 0.525532409900;
// 	pointCoordSpecial[43][0] = -0.52553240990; pointCoordSpecial[43][1] = 0.183434642500;
// 	pointCoordSpecial[44][0] = -0.52553240990; pointCoordSpecial[44][1] = -0.183434642500;
// 	pointCoordSpecial[45][0] = -0.52553240990; pointCoordSpecial[45][1] = -0.525532409900;
// 	pointCoordSpecial[46][0] = -0.52553240990; pointCoordSpecial[46][1] = -0.796666477400;
// 	pointCoordSpecial[47][0] = -0.52553240990; pointCoordSpecial[47][1] = -0.960298565000;
// 	pointCoordSpecial[48][0] = -0.79666647740; pointCoordSpecial[48][1] = 0.960298565000;
// 	pointCoordSpecial[49][0] = -0.79666647740; pointCoordSpecial[49][1] = 0.796666477400;
// 	pointCoordSpecial[50][0] = -0.79666647740; pointCoordSpecial[50][1] = 0.525532409900;
// 	pointCoordSpecial[51][0] = -0.79666647740; pointCoordSpecial[51][1] = 0.183434642500;
// 	pointCoordSpecial[52][0] = -0.79666647740; pointCoordSpecial[52][1] = -0.183434642500;
// 	pointCoordSpecial[53][0] = -0.79666647740; pointCoordSpecial[53][1] = -0.525532409900;
// 	pointCoordSpecial[54][0] = -0.79666647740; pointCoordSpecial[54][1] = -0.796666477400;
// 	pointCoordSpecial[55][0] = -0.79666647740; pointCoordSpecial[55][1] = -0.960298565000;
// 	pointCoordSpecial[56][0] = -0.96029856500; pointCoordSpecial[56][1] = 0.960298565000;
// 	pointCoordSpecial[57][0] = -0.96029856500; pointCoordSpecial[57][1] = 0.796666477400;
// 	pointCoordSpecial[58][0] = -0.96029856500; pointCoordSpecial[58][1] = 0.525532409900;
// 	pointCoordSpecial[59][0] = -0.96029856500; pointCoordSpecial[59][1] = 0.183434642500;
// 	pointCoordSpecial[60][0] = -0.96029856500; pointCoordSpecial[60][1] = -0.183434642500;
// 	pointCoordSpecial[61][0] = -0.96029856500; pointCoordSpecial[61][1] = -0.525532409900;
// 	pointCoordSpecial[62][0] = -0.96029856500; pointCoordSpecial[62][1] = -0.796666477400;
// 	pointCoordSpecial[63][0] = -0.96029856500; pointCoordSpecial[63][1] = -0.960298565000;
   
//     return pointCoordSpecial[i][j];
// };


template<>
double SpecialIntegQuadrature<2>::PointList(int i, int j){

    pointCoordSpecial[0][0] = 0.;
    pointCoordSpecial[0][1] = 0.;
        
    pointCoordSpecial[1][0] = 0.774596669241483;
    pointCoordSpecial[1][1] = 0.000000000000000;
      
    pointCoordSpecial[2][0] = -0.774596669241483;
    pointCoordSpecial[2][1] = 0.;
      
    pointCoordSpecial[3][0] = 0.;
    pointCoordSpecial[3][1] = 0.774596669241483;
     
    pointCoordSpecial[4][0] =  0.774596669241483;
    pointCoordSpecial[4][1] =  0.774596669241483;
      
    pointCoordSpecial[5][0] = -0.774596669241483;
    pointCoordSpecial[5][1] =  0.774596669241483;

    pointCoordSpecial[6][0] =  0.;
    pointCoordSpecial[6][1] = -0.774596669241483;

    pointCoordSpecial[7][0] =  0.774596669241483;
    pointCoordSpecial[7][1] = -0.774596669241483;
      
    pointCoordSpecial[8][0] = -0.774596669241483;
    pointCoordSpecial[8][1] = -0.774596669241483;
        
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
   
    return pointCoordSpecial[i][j];
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
double SpecialIntegQuadrature<2>::WeightList(int i){

     pointWeightSpecial[0] = 0.790123456790124;
    pointWeightSpecial[1] = 0.493827160493828;
    pointWeightSpecial[2] = 0.493827160493828;
    pointWeightSpecial[3] = 0.493827160493828;
    pointWeightSpecial[4] = 0.308641975308642;
    pointWeightSpecial[5] = 0.308641975308642;
    pointWeightSpecial[6] = 0.493827160493828;
    pointWeightSpecial[7] = 0.308641975308642;
    pointWeightSpecial[8] = 0.308641975308642;

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

    return pointWeightSpecial[i];
};

//------------------------------------------------------------------------------
//-----------COMPUTES THE VALUE INTERPOLATED IN THE INTEGRATION POINT-----------
//------------------------------------------------------------------------------

template<>
double SpecialIntegQuadrature<2>::interpolateQuadraticVariable(double *nValues, int &point, double *wpc, int *INC_, std::vector<IParameters *> &iparam, int &patch) {
    
    iparameters = &iparam;
    QuadShapeFunction<2>    shapeQuad;
    int dim = 2;

    double xsi[dim];
    for (int i = 0; i <dim ;i++) xsi[i] = PointList(point,i);

    double phi_[9];
    
    shapeQuad.evaluateIso(xsi,phi_,wpc,INC_,(*iparameters),patch);
    
    double int_value = 0.;
    for (int i = 0; i < 9; i++){
        int_value += nValues[i] * phi_[i];
    };

    return int_value;
};




#endif
