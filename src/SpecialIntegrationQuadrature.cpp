
#include "SpecialIntegrationQuadrature.h"


//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------QUADRATURE POINTS - COORDINATES------------------------
//------------------------------------------------------------------------------
template<>
double SpecialIntegQuadrature<2>::PointListFem(int i, int j){
   
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

    return pointCoordSpecialFem[i][j];
};

template<>
double SpecialIntegQuadrature<3>::PointListFem(int i, int j){

    pointCoordSpecialFem[0][0] = 1. / 4.;
    pointCoordSpecialFem[0][1] = 1. / 4.;
    pointCoordSpecialFem[0][2] = 1. / 4.;

    pointCoordSpecialFem[1][0] = 0.;
    pointCoordSpecialFem[1][1] = 1./3.;
    pointCoordSpecialFem[1][2] = 1./3.;

    pointCoordSpecialFem[2][0] = 1./3.;
    pointCoordSpecialFem[2][1] = 1./3.;
    pointCoordSpecialFem[2][2] = 1./3.;

    pointCoordSpecialFem[3][0] = 1./3.;
    pointCoordSpecialFem[3][1] = 1./3.;
    pointCoordSpecialFem[3][2] = 0.;

    pointCoordSpecialFem[4][0] = 1./3.;
    pointCoordSpecialFem[4][1] = 0.;
    pointCoordSpecialFem[4][2] = 1./3.;

    pointCoordSpecialFem[5][0] = 8./11.;
    pointCoordSpecialFem[5][1] = 1./11.;
    pointCoordSpecialFem[5][2] = 1./11.;

    pointCoordSpecialFem[6][0] = 1./11.;
    pointCoordSpecialFem[6][1] = 1./11.;
    pointCoordSpecialFem[6][2] = 1./11.;

    pointCoordSpecialFem[7][0] = 1./11.;
    pointCoordSpecialFem[7][1] = 1./11.;
    pointCoordSpecialFem[7][2] = 8./11.;

    pointCoordSpecialFem[8][0] = 1./11.;
    pointCoordSpecialFem[8][1] = 8./11.;
    pointCoordSpecialFem[8][2] = 1./11.;

    pointCoordSpecialFem[9][0] = 0.4334498464263357;
    pointCoordSpecialFem[9][1] = 0.0665501535736643;
    pointCoordSpecialFem[9][2] = 0.0665501535736643;

    pointCoordSpecialFem[10][0] = 0.0665501535736643;
    pointCoordSpecialFem[10][1] = 0.4334498464263357;
    pointCoordSpecialFem[10][2] = 0.0665501535736643;

    pointCoordSpecialFem[11][0] = 0.0665501535736643;
    pointCoordSpecialFem[11][1] = 0.0665501535736643;
    pointCoordSpecialFem[11][2] = 0.4334498464263357;

    pointCoordSpecialFem[12][0] = 0.0665501535736643;
    pointCoordSpecialFem[12][1] = 0.4334498464263357;
    pointCoordSpecialFem[12][2] = 0.4334498464263357;

    pointCoordSpecialFem[13][0] = 0.4334498464263357;
    pointCoordSpecialFem[13][1] = 0.0665501535736643;
    pointCoordSpecialFem[13][2] = 0.4334498464263357;

    pointCoordSpecialFem[14][0] = 0.4334498464263357;
    pointCoordSpecialFem[14][1] = 0.4334498464263357;
    pointCoordSpecialFem[14][2] = 0.0665501535736643;

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
        
    // pointCoordSpecialIso[0][0] = 0.0000000000;  pointCoordSpecialIso[0][1] = 0.0000000000;
    // pointCoordSpecialIso[1][0] = 0.0000000000;  pointCoordSpecialIso[1][1] = 0.5384693101;
    // pointCoordSpecialIso[2][0] = 0.0000000000;  pointCoordSpecialIso[2][1] = -0.5384693101;
    // pointCoordSpecialIso[3][0] = 0.0000000000;  pointCoordSpecialIso[3][1] = 0.9061798459;
    // pointCoordSpecialIso[4][0] = 0.0000000000;  pointCoordSpecialIso[4][1] = -0.9061798459;
    // pointCoordSpecialIso[5][0] = 0.5384693101;  pointCoordSpecialIso[5][1] = 0.0000000000;
    // pointCoordSpecialIso[6][0] = 0.5384693101;  pointCoordSpecialIso[6][1] = 0.5384693101;
    // pointCoordSpecialIso[7][0] = 0.5384693101;  pointCoordSpecialIso[7][1] = -0.5384693101;
    // pointCoordSpecialIso[8][0] = 0.5384693101;  pointCoordSpecialIso[8][1] = 0.9061798459;
    // pointCoordSpecialIso[9][0] = 0.5384693101;  pointCoordSpecialIso[9][1] = -0.9061798459;
    // pointCoordSpecialIso[10][0] = -0.5384693101; pointCoordSpecialIso[10][1] = 0.0000000000;
    // pointCoordSpecialIso[11][0] = -0.5384693101; pointCoordSpecialIso[11][1] = 0.5384693101;
    // pointCoordSpecialIso[12][0] = -0.5384693101; pointCoordSpecialIso[12][1] = -0.5384693101;
    // pointCoordSpecialIso[13][0] = -0.5384693101; pointCoordSpecialIso[13][1] = 0.9061798459;
    // pointCoordSpecialIso[14][0] = -0.5384693101; pointCoordSpecialIso[14][1] = -0.9061798459;
    // pointCoordSpecialIso[15][0] = 0.9061798459; pointCoordSpecialIso[15][1] = 0.0000000000;
    // pointCoordSpecialIso[16][0] = 0.9061798459; pointCoordSpecialIso[16][1] = 0.5384693101;
    // pointCoordSpecialIso[17][0] = 0.9061798459; pointCoordSpecialIso[17][1] = -0.5384693101;
    // pointCoordSpecialIso[18][0] = 0.9061798459; pointCoordSpecialIso[18][1] = 0.9061798459;
    // pointCoordSpecialIso[19][0] = 0.9061798459; pointCoordSpecialIso[19][1] = -0.9061798459;
    // pointCoordSpecialIso[20][0] = -0.9061798459; pointCoordSpecialIso[20][1] = 0.0000000000;
    // pointCoordSpecialIso[21][0] = -0.9061798459; pointCoordSpecialIso[21][1] = 0.5384693101;
    // pointCoordSpecialIso[22][0] = -0.9061798459; pointCoordSpecialIso[22][1] = -0.5384693101;
    // pointCoordSpecialIso[23][0] = -0.9061798459; pointCoordSpecialIso[23][1] = 0.9061798459;
    // pointCoordSpecialIso[24][0] = -0.9061798459; pointCoordSpecialIso[24][1] = -0.9061798459;
   
    return pointCoordSpecialIso[i][j];
};

template<>
double SpecialIntegQuadrature<3>::PointListIso(int i, int j){

    pointCoordSpecialIso[0][0] =   0.;
    pointCoordSpecialIso[0][1] =   0.;
    pointCoordSpecialIso[0][2] =   0.;   
   
    pointCoordSpecialIso[1][0] =   0.774596669241483;
    pointCoordSpecialIso[1][1] =   0.000000000000000;
    pointCoordSpecialIso[1][2] =   0.;
    
    pointCoordSpecialIso[2][0] =  -0.774596669241483;
    pointCoordSpecialIso[2][1] =   0.;
    pointCoordSpecialIso[2][2] =   0.; 
    
    pointCoordSpecialIso[3][0] =   0.;
    pointCoordSpecialIso[3][1] =   0.774596669241483;
    pointCoordSpecialIso[3][2] =   0.;
    
    pointCoordSpecialIso[4][0] =   0.774596669241483;
    pointCoordSpecialIso[4][1] =   0.774596669241483;
    pointCoordSpecialIso[4][2] =   0.;
    
    pointCoordSpecialIso[5][0] =  -0.774596669241483;
    pointCoordSpecialIso[5][1] =   0.774596669241483;
    pointCoordSpecialIso[5][2] =   0.;
    
    pointCoordSpecialIso[6][0] =   0.;
    pointCoordSpecialIso[6][1] =  -0.774596669241483;
    pointCoordSpecialIso[6][2] =   0.;
    
    pointCoordSpecialIso[7][0] =   0.774596669241483;
    pointCoordSpecialIso[7][1] =   -0.774596669241483;
    pointCoordSpecialIso[7][2] =   0.;  
    
    pointCoordSpecialIso[8][0] =  -0.774596669241483;
    pointCoordSpecialIso[8][1] =  -0.774596669241483;
    pointCoordSpecialIso[8][2] =   0.;
    
    pointCoordSpecialIso[9][0] =   0.;
    pointCoordSpecialIso[9][1] =   0.;
    pointCoordSpecialIso[9][2] =   0.774596669241483;  
    
    pointCoordSpecialIso[10][0] =  0.774596669241483;
    pointCoordSpecialIso[10][1] =  0.000000000000000;
    pointCoordSpecialIso[10][2] =  0.774596669241483;
    
    pointCoordSpecialIso[11][0] = -0.774596669241483;
    pointCoordSpecialIso[11][1] =  0.;
    pointCoordSpecialIso[11][2] =  0.774596669241483;
    
    pointCoordSpecialIso[12][0] =  0.;
    pointCoordSpecialIso[12][1] =  0.774596669241483;
    pointCoordSpecialIso[12][2] =  0.774596669241483;
    
    pointCoordSpecialIso[13][0] =  0.774596669241483;
    pointCoordSpecialIso[13][1] =  0.774596669241483;
    pointCoordSpecialIso[13][2] =  0.774596669241483;
    
    pointCoordSpecialIso[14][0] = -0.774596669241483;
    pointCoordSpecialIso[14][1] =  0.774596669241483;
    pointCoordSpecialIso[14][2] =  0.774596669241483;
    
    pointCoordSpecialIso[15][0] =  0.;
    pointCoordSpecialIso[15][1] = -0.774596669241483;
    pointCoordSpecialIso[15][2] =  0.774596669241483;
    
    pointCoordSpecialIso[16][0] =  0.774596669241483;
    pointCoordSpecialIso[16][1] = -0.774596669241483;
    pointCoordSpecialIso[16][2] =  0.774596669241483;
    
    pointCoordSpecialIso[17][0] = -0.774596669241483;
    pointCoordSpecialIso[17][1] = -0.774596669241483;
    pointCoordSpecialIso[17][2] =  0.774596669241483;
    
    pointCoordSpecialIso[18][0] =  0.;
    pointCoordSpecialIso[18][1] =  0.;
    pointCoordSpecialIso[18][2] = -0.774596669241483; 
    
    pointCoordSpecialIso[19][0] =  0.774596669241483;
    pointCoordSpecialIso[19][1] =  0.000000000000000;
    pointCoordSpecialIso[19][2] = -0.774596669241483;
   
    pointCoordSpecialIso[20][0] = -0.774596669241483;
    pointCoordSpecialIso[20][1] =  0.;
    pointCoordSpecialIso[20][2] = -0.774596669241483;
    
    pointCoordSpecialIso[21][0] =  0.;
    pointCoordSpecialIso[21][1] =  0.774596669241483;
    pointCoordSpecialIso[21][2] = -0.774596669241483;
    
    pointCoordSpecialIso[22][0] =  0.774596669241483;
    pointCoordSpecialIso[22][1] =  0.774596669241483;
    pointCoordSpecialIso[22][2] = -0.774596669241483;  
    
    pointCoordSpecialIso[23][0] = -0.774596669241483;
    pointCoordSpecialIso[23][1] =  0.774596669241483;
    pointCoordSpecialIso[23][2] = -0.774596669241483;
    
    pointCoordSpecialIso[24][0] =  0.;
    pointCoordSpecialIso[24][1] = -0.774596669241483;
    pointCoordSpecialIso[24][2] = -0.774596669241483;
    
    pointCoordSpecialIso[25][0] =  0.774596669241483;
    pointCoordSpecialIso[25][1] = -0.774596669241483;
    pointCoordSpecialIso[25][2] = -0.774596669241483;    
    
    pointCoordSpecialIso[26][0] = -0.774596669241483;
    pointCoordSpecialIso[26][1] = -0.774596669241483;
    pointCoordSpecialIso[26][2] = -0.774596669241483;

    return pointCoordSpecialIso[i][j];
};


//------------------------------------------------------------------------------
//-------------------------QUADRATURE POINTS - WEIGHTS--------------------------
//------------------------------------------------------------------------------
template<>
double SpecialIntegQuadrature<2>::WeightListFem(int i){

    pointWeightSpecialFem[0] = 0.11250;
    pointWeightSpecialFem[1] = (155. - sqrt(15.)) / 2400.;
    pointWeightSpecialFem[2] = (155. - sqrt(15.)) / 2400.;
    pointWeightSpecialFem[3] = (155. - sqrt(15.)) / 2400.;
    pointWeightSpecialFem[4] = (155. + sqrt(15.)) / 2400.;
    pointWeightSpecialFem[5] = (155. + sqrt(15.)) / 2400.;
    pointWeightSpecialFem[6] = (155. + sqrt(15.)) / 2400.; 

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

    return pointWeightSpecialFem[i];

};

template<>
double SpecialIntegQuadrature<3>::WeightListFem(int i){
    
    pointWeightSpecialFem[0] = 0.1817020685825351/6.;
    pointWeightSpecialFem[1] = 0.0361607142857143/6.;
    pointWeightSpecialFem[2] = 0.0361607142857143/6.;
    pointWeightSpecialFem[3] = 0.0361607142857143/6.;
    pointWeightSpecialFem[4] = 0.0361607142857143/6.;
    pointWeightSpecialFem[5] = 0.0698714945161738/6.;
    pointWeightSpecialFem[6] = 0.0698714945161738/6.; 
    pointWeightSpecialFem[7] = 0.0698714945161738/6.;
    pointWeightSpecialFem[8] = 0.0698714945161738/6.; 
    pointWeightSpecialFem[9] = 0.0656948493683187/6.;
    pointWeightSpecialFem[10] = 0.0656948493683187/6.; 
    pointWeightSpecialFem[11] = 0.0656948493683187/6.;
    pointWeightSpecialFem[12] = 0.0656948493683187/6.;
    pointWeightSpecialFem[13] = 0.0656948493683187/6.;
    pointWeightSpecialFem[14] = 0.0656948493683187/6.;

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

    // pointWeightSpecialIso[0] = 0.3236345679;
    // pointWeightSpecialIso[1] = 0.2722865326;
    // pointWeightSpecialIso[2] = 0.2722865326;
    // pointWeightSpecialIso[3] = 0.1347850724;
    // pointWeightSpecialIso[4] = 0.1347850724;
    // pointWeightSpecialIso[5] = 0.2722865326;
    // pointWeightSpecialIso[6] = 0.2290854042;
    // pointWeightSpecialIso[7] = 0.2290854042;
    // pointWeightSpecialIso[8] = 0.1134000000;
    // pointWeightSpecialIso[9] = 0.1134000000;
    // pointWeightSpecialIso[10] = 0.2722865326;
    // pointWeightSpecialIso[11] = 0.2290854042;
    // pointWeightSpecialIso[12] = 0.2290854042;
    // pointWeightSpecialIso[13] = 0.1134000000;
    // pointWeightSpecialIso[14] = 0.1134000000;
    // pointWeightSpecialIso[15] = 0.1347850724;
    // pointWeightSpecialIso[16] = 0.1134000000;
    // pointWeightSpecialIso[17] = 0.1134000000;
    // pointWeightSpecialIso[18] = 0.0561343489;
    // pointWeightSpecialIso[19] = 0.0561343489;
    // pointWeightSpecialIso[20] = 0.1347850724;
    // pointWeightSpecialIso[21] = 0.1134000000;
    // pointWeightSpecialIso[22] = 0.1134000000;
    // pointWeightSpecialIso[23] = 0.0561343489;
    // pointWeightSpecialIso[24] = 0.0561343489;

    return pointWeightSpecialIso[i];
};

template<>
double SpecialIntegQuadrature<3>::WeightListIso(int i){
    
    pointWeightSpecialIso[0] = 0.702331961591221;
    pointWeightSpecialIso[1] = 0.438957475994514;
    pointWeightSpecialIso[2] = 0.438957475994514;
    pointWeightSpecialIso[3] = 0.438957475994514;
    pointWeightSpecialIso[4] = 0.274348422496571;
    pointWeightSpecialIso[5] = 0.274348422496571;
    pointWeightSpecialIso[6] = 0.438957475994514;
    pointWeightSpecialIso[7] = 0.274348422496571;
    pointWeightSpecialIso[8] = 0.274348422496571;
    pointWeightSpecialIso[9] = 0.438957475994514;
    pointWeightSpecialIso[10] = 0.274348422496571;
    pointWeightSpecialIso[11] = 0.274348422496571;
    pointWeightSpecialIso[12] = 0.274348422496571;
    pointWeightSpecialIso[13] = 0.171467764060357;
    pointWeightSpecialIso[14] = 0.171467764060357;
    pointWeightSpecialIso[15] = 0.274348422496571;
    pointWeightSpecialIso[16] = 0.171467764060357;
    pointWeightSpecialIso[17] = 0.171467764060357;
    pointWeightSpecialIso[18] = 0.438957475994514;
    pointWeightSpecialIso[19] = 0.274348422496571;
    pointWeightSpecialIso[20] = 0.274348422496571;
    pointWeightSpecialIso[21] = 0.274348422496571;
    pointWeightSpecialIso[22] = 0.171467764060357;
    pointWeightSpecialIso[23] = 0.171467764060357;
    pointWeightSpecialIso[24] = 0.274348422496571;
    pointWeightSpecialIso[25] = 0.171467764060357;
    pointWeightSpecialIso[26] = 0.171467764060357;

    return pointWeightSpecialIso[i];
};

//------------------------------------------------------------------------------
//-----------COMPUTES THE VALUE INTERPOLATED IN THE INTEGRATION POINT-----------
//------------------------------------------------------------------------------
template<int DIM>
double SpecialIntegQuadrature<DIM>::interpolateQuadraticVariableFem(double *nValues, int &point) {
    
    QuadShapeFunction<DIM>    shapeQuad;

    double xsi[DIM];
    for (int i = 0; i < DIM; i++) xsi[i] = PointListFem(point,i);
    
    int LNN = 4*DIM-2;;
    double phi_[LNN];
    shapeQuad.evaluateFem(xsi,phi_);
    
    double int_value = 0.;
    for (int i = 0; i < LNN; i++){
        int_value += nValues[i] * phi_[i];
    };
    return int_value;
};


template<int DIM>
double SpecialIntegQuadrature<DIM>::interpolateQuadraticVariableIso(double *nValues, int &point, double *wpc, int *INC_, std::vector<IParameters *> &iparam, int &patch) {
    
    iparameters = &iparam;
    QuadShapeFunction<DIM>    shapeQuad;

    double xsi[DIM];
    for (int i = 0; i <DIM ;i++) xsi[i] = PointListIso(point,i);

    int LNN = 18*DIM -27;
    double phi_[LNN];
    
    shapeQuad.evaluateIso(xsi,phi_,wpc,INC_,(*iparameters),patch);
    
    double int_value = 0.;
    for (int i = 0; i < LNN; i++){
        int_value += nValues[i] * phi_[i];
    };

    return int_value;
};

template class SpecialIntegQuadrature<2>;
template class SpecialIntegQuadrature<3>;

