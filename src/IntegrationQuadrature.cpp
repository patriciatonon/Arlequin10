#include "IntegrationQuadrature.h"

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
double IntegQuadrature<3>::PointListFem(int i, int j){



        pointCoordFem[0][0] = 0.310885919263300609797345733763458e+0;
    pointCoordFem[0][1] = 0.310885919263300609797345733763458e+0;
    pointCoordFem[0][2] = 0.673422422100981706079627987096265e-1;

    pointCoordFem[1][0] = 0.310885919263300609797345733763458e+0;
    pointCoordFem[1][1] = 0.673422422100981706079627987096265e-1;
    pointCoordFem[1][2] = 0.310885919263300609797345733763458e+0;

    pointCoordFem[2][0] = 0.673422422100981706079627987096265e-1;
    pointCoordFem[2][1] = 0.310885919263300609797345733763458e+0;
    pointCoordFem[2][2] = 0.310885919263300609797345733763458e+0;

    pointCoordFem[3][0] = 0.310885919263300609797345733763458e+0;
    pointCoordFem[3][1] = 0.310885919263300609797345733763458e+0;
    pointCoordFem[3][2] = 0.310885919263300609797345733763458e+0;

    pointCoordFem[4][0] = 0.927352503108912264023239137370306e-1;
    pointCoordFem[4][1] = 0.927352503108912264023239137370306e-1;
    pointCoordFem[4][2] = 0.721794249067326320793028258788908e+0;

    pointCoordFem[5][0] = 0.927352503108912264023239137370306e-1;
    pointCoordFem[5][1] = 0.721794249067326320793028258788908e+0;
    pointCoordFem[5][2] = 0.927352503108912264023239137370306e-1;

    pointCoordFem[6][0] = 0.721794249067326320793028258788908e+0;
    pointCoordFem[6][1] = 0.927352503108912264023239137370306e-1;
    pointCoordFem[6][2] = 0.927352503108912264023239137370306e-1;

    pointCoordFem[7][0] = 0.927352503108912264023239137370306e-1;
    pointCoordFem[7][1] = 0.927352503108912264023239137370306e-1;
    pointCoordFem[7][2] = 0.927352503108912264023239137370306e-1;

    pointCoordFem[8][0] = 0.455037041256496494918805262793394e-1;
    pointCoordFem[8][1] = 0.454496295874350350508119473720661e+0;
    pointCoordFem[8][2] = 0.454496295874350350508119473720661e+0;

    pointCoordFem[9][0] = 0.454496295874350350508119473720661e+0;
    pointCoordFem[9][1] = 0.455037041256496494918805262793394e-1;
    pointCoordFem[9][2] = 0.454496295874350350508119473720661e+0;

    pointCoordFem[10][0] = 0.455037041256496494918805262793394e-1;
    pointCoordFem[10][1] = 0.455037041256496494918805262793394e-1;
    pointCoordFem[10][2] = 0.454496295874350350508119473720661e+0;

    pointCoordFem[11][0] = 0.455037041256496494918805262793394e-1;
    pointCoordFem[11][1] = 0.454496295874350350508119473720661e+0;
    pointCoordFem[11][2] = 0.455037041256496494918805262793394e-1;

    pointCoordFem[12][0] = 0.454496295874350350508119473720661e+0;
    pointCoordFem[12][1] = 0.455037041256496494918805262793394e-1;
    pointCoordFem[12][2] = 0.455037041256496494918805262793394e-1;

    pointCoordFem[13][0] = 0.454496295874350350508119473720661e+0;
    pointCoordFem[13][1] = 0.454496295874350350508119473720661e+0;
    pointCoordFem[13][2] = 0.455037041256496494918805262793394e-1;

    //     const double a = (1. + sqrt(5. / 14.)) / 4.;
    // const double b = (1. - sqrt(5. / 14.)) / 4.;

    // pointCoordFem[0][0] = 1. / 4.;
    // pointCoordFem[0][1] = 1. / 4.;
    // pointCoordFem[0][2] = 1. / 4.;

    // pointCoordFem[1][0] = 11. / 14.;
    // pointCoordFem[1][1] = 1. / 14.;
    // pointCoordFem[1][2] = 1. / 14.;

    // pointCoordFem[2][0] = 1. / 14.;
    // pointCoordFem[2][1] = 11. / 14.;
    // pointCoordFem[2][2] = 1. / 14.;

    // pointCoordFem[3][0] = 1. / 14.;
    // pointCoordFem[3][1] = 1. / 14.;
    // pointCoordFem[3][2] = 11. / 14.;

    // pointCoordFem[4][0] = 1. / 14.;
    // pointCoordFem[4][1] = 1. / 14.;
    // pointCoordFem[4][2] = 1. / 14.;

    // pointCoordFem[5][0] = a;
    // pointCoordFem[5][1] = a;
    // pointCoordFem[5][2] = b;

    // pointCoordFem[6][0] = a;
    // pointCoordFem[6][1] = b;
    // pointCoordFem[6][2] = a;

    // pointCoordFem[7][0] = a;
    // pointCoordFem[7][1] = b;
    // pointCoordFem[7][2] = b;

    // pointCoordFem[8][0] = b;
    // pointCoordFem[8][1] = a;
    // pointCoordFem[8][2] = a;

    // pointCoordFem[9][0] = b;
    // pointCoordFem[9][1] = a;
    // pointCoordFem[9][2] = b;

    // pointCoordFem[10][0] = b;
    // pointCoordFem[10][1] = b;
    // pointCoordFem[10][2] = a;
    
    // pointCoordFem[0][0] = 1. / 4.;
    // pointCoordFem[0][1] = 1. / 4.;
    // pointCoordFem[0][2] = 1. / 4.;

    // pointCoordFem[1][0] = 0.;
    // pointCoordFem[1][1] = 1./3.;
    // pointCoordFem[1][2] = 1./3.;

    // pointCoordFem[2][0] = 1./3.;
    // pointCoordFem[2][1] = 1./3.;
    // pointCoordFem[2][2] = 1./3.;

    // pointCoordFem[3][0] = 1./3.;
    // pointCoordFem[3][1] = 1./3.;
    // pointCoordFem[3][2] = 0.;

    // pointCoordFem[4][0] = 1./3.;
    // pointCoordFem[4][1] = 0.;
    // pointCoordFem[4][2] = 1./3.;

    // pointCoordFem[5][0] = 8./11.;
    // pointCoordFem[5][1] = 1./11.;
    // pointCoordFem[5][2] = 1./11.;

    // pointCoordFem[6][0] = 1./11.;
    // pointCoordFem[6][1] = 1./11.;
    // pointCoordFem[6][2] = 1./11.;

    // pointCoordFem[7][0] = 1./11.;
    // pointCoordFem[7][1] = 1./11.;
    // pointCoordFem[7][2] = 8./11.;

    // pointCoordFem[8][0] = 1./11.;
    // pointCoordFem[8][1] = 8./11.;
    // pointCoordFem[8][2] = 1./11.;

    // pointCoordFem[9][0] = 0.4334498464263357;
    // pointCoordFem[9][1] = 0.0665501535736643;
    // pointCoordFem[9][2] = 0.0665501535736643;

    // pointCoordFem[10][0] = 0.0665501535736643;
    // pointCoordFem[10][1] = 0.4334498464263357;
    // pointCoordFem[10][2] = 0.0665501535736643;

    // pointCoordFem[11][0] = 0.0665501535736643;
    // pointCoordFem[11][1] = 0.0665501535736643;
    // pointCoordFem[11][2] = 0.4334498464263357;

    // pointCoordFem[12][0] = 0.0665501535736643;
    // pointCoordFem[12][1] = 0.4334498464263357;
    // pointCoordFem[12][2] = 0.4334498464263357;

    // pointCoordFem[13][0] = 0.4334498464263357;
    // pointCoordFem[13][1] = 0.0665501535736643;
    // pointCoordFem[13][2] = 0.4334498464263357;

    // pointCoordFem[14][0] = 0.4334498464263357;
    // pointCoordFem[14][1] = 0.4334498464263357;
    // pointCoordFem[14][2] = 0.0665501535736643;

    return pointCoordFem[i][j];
}

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

template<>
double IntegQuadrature<3>::PointListIso(int i, int j){
    

    // pointCoordIso[0][0] = -0.33998104358485597976;
    // pointCoordIso[0][1] = -0.33998104358485597976;
    // pointCoordIso[0][2] = -0.33998104358485597976;

    // pointCoordIso[1][0] = -0.33998104358485597976;
    // pointCoordIso[1][1] = -0.33998104358485597976;
    // pointCoordIso[1][2] = 0.33998104358485597976;

    // pointCoordIso[2][0] = -0.33998104358485597976;
    // pointCoordIso[2][1] = -0.33998104358485597976;
    // pointCoordIso[2][2] = -0.86113631159405301663;

    // pointCoordIso[3][0] = -0.33998104358485597976;
    // pointCoordIso[3][1] = -0.33998104358485597976;
    // pointCoordIso[3][2] = 0.86113631159405301663;

    // pointCoordIso[4][0] = -0.33998104358485597976;
    // pointCoordIso[4][1] = 0.33998104358485597976;
    // pointCoordIso[4][2] = -0.33998104358485597976;

    // pointCoordIso[5][0] = -0.33998104358485597976;
    // pointCoordIso[5][1] = 0.33998104358485597976;
    // pointCoordIso[5][2] = 0.33998104358485597976;

    // pointCoordIso[6][0] = -0.33998104358485597976;
    // pointCoordIso[6][1] = 0.33998104358485597976;
    // pointCoordIso[6][2] = -0.86113631159405301663;

    // pointCoordIso[7][0] = -0.33998104358485597976;
    // pointCoordIso[7][1] = 0.33998104358485597976;
    // pointCoordIso[7][2] = 0.86113631159405301663;

    // pointCoordIso[8][0] = -0.33998104358485597976;
    // pointCoordIso[8][1] = -0.86113631159405301663;
    // pointCoordIso[8][2] = -0.33998104358485597976;

    // pointCoordIso[9][0] = -0.33998104358485597976;
    // pointCoordIso[9][1] = -0.86113631159405301663;
    // pointCoordIso[9][2] = 0.33998104358485597976;

    // pointCoordIso[10][0] = -0.33998104358485597976;
    // pointCoordIso[10][1] = -0.86113631159405301663;
    // pointCoordIso[10][2] = -0.86113631159405301663;

    // pointCoordIso[11][0] = -0.33998104358485597976;
    // pointCoordIso[11][1] = -0.86113631159405301663;
    // pointCoordIso[11][2] = 0.86113631159405301663;

    // pointCoordIso[12][0] = -0.33998104358485597976;
    // pointCoordIso[12][1] = 0.86113631159405301663;
    // pointCoordIso[12][2] = -0.33998104358485597976;

    // pointCoordIso[13][0] = -0.33998104358485597976;
    // pointCoordIso[13][1] = 0.86113631159405301663;
    // pointCoordIso[13][2] = 0.33998104358485597976;

    // pointCoordIso[14][0] = -0.33998104358485597976;
    // pointCoordIso[14][1] = 0.86113631159405301663;
    // pointCoordIso[14][2] = -0.86113631159405301663;

    // pointCoordIso[15][0] = -0.33998104358485597976;
    // pointCoordIso[15][1] = 0.86113631159405301663;
    // pointCoordIso[15][2] = 0.86113631159405301663;

    // pointCoordIso[16][0] = 0.33998104358485597976;
    // pointCoordIso[16][1] = -0.33998104358485597976;
    // pointCoordIso[16][2] = -0.33998104358485597976;

    // pointCoordIso[17][0] = 0.33998104358485597976;
    // pointCoordIso[17][1] = -0.33998104358485597976;
    // pointCoordIso[17][2] = 0.33998104358485597976;

    // pointCoordIso[18][0] = 0.33998104358485597976;
    // pointCoordIso[18][1] = -0.33998104358485597976;
    // pointCoordIso[18][2] = -0.86113631159405301663;

    // pointCoordIso[19][0] = 0.33998104358485597976;
    // pointCoordIso[19][1] = -0.33998104358485597976;
    // pointCoordIso[19][2] = 0.86113631159405301663;

    // pointCoordIso[20][0] = 0.33998104358485597976;
    // pointCoordIso[20][1] = 0.33998104358485597976;
    // pointCoordIso[20][2] = -0.33998104358485597976;

    // pointCoordIso[21][0] = 0.33998104358485597976;
    // pointCoordIso[21][1] = 0.33998104358485597976;
    // pointCoordIso[21][2] = 0.33998104358485597976;

    // pointCoordIso[22][0] = 0.33998104358485597976;
    // pointCoordIso[22][1] = 0.33998104358485597976;
    // pointCoordIso[22][2] = -0.86113631159405301663;

    // pointCoordIso[23][0] = 0.33998104358485597976;
    // pointCoordIso[23][1] = 0.33998104358485597976;
    // pointCoordIso[23][2] = 0.86113631159405301663;

    // pointCoordIso[24][0] = 0.33998104358485597976;
    // pointCoordIso[24][1] = -0.86113631159405301663;
    // pointCoordIso[24][2] = -0.33998104358485597976;

    // pointCoordIso[25][0] = 0.33998104358485597976;
    // pointCoordIso[25][1] = -0.86113631159405301663;
    // pointCoordIso[25][2] = 0.33998104358485597976;

    // pointCoordIso[26][0] = 0.33998104358485597976;
    // pointCoordIso[26][1] = -0.86113631159405301663;
    // pointCoordIso[26][2] = -0.86113631159405301663;

    // pointCoordIso[27][0] = 0.33998104358485597976;
    // pointCoordIso[27][1] = -0.86113631159405301663;
    // pointCoordIso[27][2] = 0.86113631159405301663;

    // pointCoordIso[28][0] = 0.33998104358485597976;
    // pointCoordIso[28][1] = 0.86113631159405301663;
    // pointCoordIso[28][2] = -0.33998104358485597976;

    // pointCoordIso[29][0] = 0.33998104358485597976;
    // pointCoordIso[29][1] = 0.86113631159405301663;
    // pointCoordIso[29][2] = 0.33998104358485597976;

    // pointCoordIso[30][0] = 0.33998104358485597976;
    // pointCoordIso[30][1] = 0.86113631159405301663;
    // pointCoordIso[30][2] = -0.86113631159405301663;

    // pointCoordIso[31][0] = 0.33998104358485597976;
    // pointCoordIso[31][1] = 0.86113631159405301663;
    // pointCoordIso[31][2] = 0.86113631159405301663;

    // pointCoordIso[32][0] = -0.86113631159405301663;
    // pointCoordIso[32][1] = -0.33998104358485597976;
    // pointCoordIso[32][2] = -0.33998104358485597976;

    // pointCoordIso[33][0] = -0.86113631159405301663;
    // pointCoordIso[33][1] = -0.33998104358485597976;
    // pointCoordIso[33][2] = 0.33998104358485597976;

    // pointCoordIso[34][0] = -0.86113631159405301663;
    // pointCoordIso[34][1] = -0.33998104358485597976;
    // pointCoordIso[34][2] = -0.86113631159405301663;

    // pointCoordIso[35][0] = -0.86113631159405301663;
    // pointCoordIso[35][1] = -0.33998104358485597976;
    // pointCoordIso[35][2] = 0.86113631159405301663;

    // pointCoordIso[36][0] = -0.86113631159405301663;
    // pointCoordIso[36][1] = 0.33998104358485597976;
    // pointCoordIso[36][2] = -0.33998104358485597976;

    // pointCoordIso[37][0] = -0.86113631159405301663;
    // pointCoordIso[37][1] = 0.33998104358485597976;
    // pointCoordIso[37][2] = 0.33998104358485597976;

    // pointCoordIso[38][0] = -0.86113631159405301663;
    // pointCoordIso[38][1] = 0.33998104358485597976;
    // pointCoordIso[38][2] = -0.86113631159405301663;

    // pointCoordIso[39][0] = -0.86113631159405301663;   
    // pointCoordIso[39][1] = 0.33998104358485597976;
    // pointCoordIso[39][2] = 0.86113631159405301663;

    // pointCoordIso[40][0] = -0.86113631159405301663;
    // pointCoordIso[40][1] = -0.86113631159405301663;
    // pointCoordIso[40][2] = -0.33998104358485597976;

    // pointCoordIso[41][0] = -0.86113631159405301663;
    // pointCoordIso[41][1] = -0.86113631159405301663;
    // pointCoordIso[41][2] = 0.33998104358485597976;

    // pointCoordIso[42][0] = -0.86113631159405301663;
    // pointCoordIso[42][1] = -0.86113631159405301663;
    // pointCoordIso[42][2] = -0.86113631159405301663;

    // pointCoordIso[43][0] = -0.86113631159405301663;
    // pointCoordIso[43][1] = -0.86113631159405301663;
    // pointCoordIso[43][2] = 0.86113631159405301663;

    // pointCoordIso[44][0] = -0.86113631159405301663;
    // pointCoordIso[44][1] = 0.86113631159405301663;
    // pointCoordIso[44][2] = -0.33998104358485597976;

    // pointCoordIso[45][0] = -0.86113631159405301663;
    // pointCoordIso[45][1] = 0.86113631159405301663;
    // pointCoordIso[45][2] = 0.33998104358485597976;

    // pointCoordIso[46][0] = -0.86113631159405301663;
    // pointCoordIso[46][1] = 0.86113631159405301663;
    // pointCoordIso[46][2] = -0.86113631159405301663;

    // pointCoordIso[47][0] = -0.86113631159405301663;
    // pointCoordIso[47][1] = 0.86113631159405301663;
    // pointCoordIso[47][2] = 0.86113631159405301663;

    // pointCoordIso[48][0] = 0.86113631159405301663;
    // pointCoordIso[48][1] = -0.33998104358485597976;
    // pointCoordIso[48][2] = -0.33998104358485597976;

    // pointCoordIso[49][0] = 0.86113631159405301663;
    // pointCoordIso[49][1] = -0.33998104358485597976;
    // pointCoordIso[49][2] = 0.33998104358485597976;

    // pointCoordIso[50][0] = 0.86113631159405301663;
    // pointCoordIso[50][1] = -0.33998104358485597976;
    // pointCoordIso[50][2] = -0.86113631159405301663;

    // pointCoordIso[51][0] = 0.86113631159405301663;
    // pointCoordIso[51][1] = -0.33998104358485597976;
    // pointCoordIso[51][2] = 0.86113631159405301663;

    // pointCoordIso[52][0] = 0.86113631159405301663;
    // pointCoordIso[52][1] = 0.33998104358485597976;
    // pointCoordIso[52][2] = -0.33998104358485597976;

    // pointCoordIso[53][0] = 0.86113631159405301663;
    // pointCoordIso[53][1] = 0.33998104358485597976;
    // pointCoordIso[53][2] = 0.33998104358485597976;

    // pointCoordIso[54][0] = 0.86113631159405301663;
    // pointCoordIso[54][1] = 0.33998104358485597976;
    // pointCoordIso[54][2] = -0.86113631159405301663;

    // pointCoordIso[55][0] = 0.86113631159405301663;
    // pointCoordIso[55][1] = 0.33998104358485597976;
    // pointCoordIso[55][2] = 0.86113631159405301663;

    // pointCoordIso[56][0] = 0.86113631159405301663;
    // pointCoordIso[56][1] = -0.86113631159405301663;
    // pointCoordIso[56][2] = -0.33998104358485597976;

    // pointCoordIso[57][0] = 0.86113631159405301663;
    // pointCoordIso[57][1] = -0.86113631159405301663;
    // pointCoordIso[57][2] = 0.33998104358485597976;

    // pointCoordIso[58][0] = 0.86113631159405301663;
    // pointCoordIso[58][1] = -0.86113631159405301663;
    // pointCoordIso[58][2] = -0.86113631159405301663;

    // pointCoordIso[59][0] = 0.86113631159405301663;
    // pointCoordIso[59][1] = -0.86113631159405301663;
    // pointCoordIso[59][2] = 0.86113631159405301663;

    // pointCoordIso[60][0] = 0.86113631159405301663;
    // pointCoordIso[60][1] = 0.86113631159405301663;
    // pointCoordIso[60][2] = -0.33998104358485597976;

    // pointCoordIso[61][0] = 0.86113631159405301663;
    // pointCoordIso[61][1] = 0.86113631159405301663;
    // pointCoordIso[61][2] = 0.33998104358485597976;

    // pointCoordIso[62][0] = 0.86113631159405301663;
    // pointCoordIso[62][1] = 0.86113631159405301663;
    // pointCoordIso[62][2] = -0.86113631159405301663;

    // pointCoordIso[63][0] = 0.86113631159405301663;
    // pointCoordIso[63][1] = 0.86113631159405301663;
    // pointCoordIso[63][2] = 0.86113631159405301663;
   
    pointCoordIso[0][0] =   0.;
    pointCoordIso[0][1] =   0.;
    pointCoordIso[0][2] =   0.;   
   
    pointCoordIso[1][0] =   0.774596669241483;
    pointCoordIso[1][1] =   0.000000000000000;
    pointCoordIso[1][2] =   0.;
    
    pointCoordIso[2][0] =  -0.774596669241483;
    pointCoordIso[2][1] =   0.;
    pointCoordIso[2][2] =   0.; 
    
    pointCoordIso[3][0] =   0.;
    pointCoordIso[3][1] =   0.774596669241483;
    pointCoordIso[3][2] =   0.;
    
    pointCoordIso[4][0] =   0.774596669241483;
    pointCoordIso[4][1] =   0.774596669241483;
    pointCoordIso[4][2] =   0.;
    
    pointCoordIso[5][0] =  -0.774596669241483;
    pointCoordIso[5][1] =   0.774596669241483;
    pointCoordIso[5][2] =   0.;
    
    pointCoordIso[6][0] =   0.;
    pointCoordIso[6][1] =  -0.774596669241483;
    pointCoordIso[6][2] =   0.;
    
    pointCoordIso[7][0] =   0.774596669241483;
    pointCoordIso[7][1] =   -0.774596669241483;
    pointCoordIso[7][2] =   0.;  
    
    pointCoordIso[8][0] =  -0.774596669241483;
    pointCoordIso[8][1] =  -0.774596669241483;
    pointCoordIso[8][2] =   0.;
    
    pointCoordIso[9][0] =   0.;
    pointCoordIso[9][1] =   0.;
    pointCoordIso[9][2] =   0.774596669241483;  
    
    pointCoordIso[10][0] =  0.774596669241483;
    pointCoordIso[10][1] =  0.000000000000000;
    pointCoordIso[10][2] =  0.774596669241483;
    
    pointCoordIso[11][0] = -0.774596669241483;
    pointCoordIso[11][1] =  0.;
    pointCoordIso[11][2] =  0.774596669241483;
    
    pointCoordIso[12][0] =  0.;
    pointCoordIso[12][1] =  0.774596669241483;
    pointCoordIso[12][2] =  0.774596669241483;
    
    pointCoordIso[13][0] =  0.774596669241483;
    pointCoordIso[13][1] =  0.774596669241483;
    pointCoordIso[13][2] =  0.774596669241483;
    
    pointCoordIso[14][0] = -0.774596669241483;
    pointCoordIso[14][1] =  0.774596669241483;
    pointCoordIso[14][2] =  0.774596669241483;
    
    pointCoordIso[15][0] =  0.;
    pointCoordIso[15][1] = -0.774596669241483;
    pointCoordIso[15][2] =  0.774596669241483;
    
    pointCoordIso[16][0] =  0.774596669241483;
    pointCoordIso[16][1] = -0.774596669241483;
    pointCoordIso[16][2] =  0.774596669241483;
    
    pointCoordIso[17][0] = -0.774596669241483;
    pointCoordIso[17][1] = -0.774596669241483;
    pointCoordIso[17][2] =  0.774596669241483;
    
    pointCoordIso[18][0] =  0.;
    pointCoordIso[18][1] =  0.;
    pointCoordIso[18][2] = -0.774596669241483; 
    
    pointCoordIso[19][0] =  0.774596669241483;
    pointCoordIso[19][1] =  0.000000000000000;
    pointCoordIso[19][2] = -0.774596669241483;
   
    pointCoordIso[20][0] = -0.774596669241483;
    pointCoordIso[20][1] =  0.;
    pointCoordIso[20][2] = -0.774596669241483;
    
    pointCoordIso[21][0] =  0.;
    pointCoordIso[21][1] =  0.774596669241483;
    pointCoordIso[21][2] = -0.774596669241483;
    
    pointCoordIso[22][0] =  0.774596669241483;
    pointCoordIso[22][1] =  0.774596669241483;
    pointCoordIso[22][2] = -0.774596669241483;  
    
    pointCoordIso[23][0] = -0.774596669241483;
    pointCoordIso[23][1] =  0.774596669241483;
    pointCoordIso[23][2] = -0.774596669241483;
    
    pointCoordIso[24][0] =  0.;
    pointCoordIso[24][1] = -0.774596669241483;
    pointCoordIso[24][2] = -0.774596669241483;
    
    pointCoordIso[25][0] =  0.774596669241483;
    pointCoordIso[25][1] = -0.774596669241483;
    pointCoordIso[25][2] = -0.774596669241483;    
    
    pointCoordIso[26][0] = -0.774596669241483;
    pointCoordIso[26][1] = -0.774596669241483;
    pointCoordIso[26][2] = -0.774596669241483;

    return pointCoordIso[i][j];
}

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
double IntegQuadrature<3>::WeightListFem(int i){


    pointWeightFem[0] = 0.187813209530026417998642753888811e-1;
    pointWeightFem[1] = 0.187813209530026417998642753888811e-1;
    pointWeightFem[2] = 0.187813209530026417998642753888811e-1;
    pointWeightFem[3] = 0.187813209530026417998642753888811e-1;
    pointWeightFem[4] = 0.122488405193936582572850342477213e-1;
    pointWeightFem[5] = 0.122488405193936582572850342477213e-1;
    pointWeightFem[6] = 0.122488405193936582572850342477213e-1;
    pointWeightFem[7] = 0.122488405193936582572850342477213e-1;
    pointWeightFem[8] = 0.709100346284691107301157135337624e-2;
    pointWeightFem[9] = 0.709100346284691107301157135337624e-2;
    pointWeightFem[10] = 0.709100346284691107301157135337624e-2;
    pointWeightFem[11] = 0.709100346284691107301157135337624e-2;
    pointWeightFem[12] = 0.709100346284691107301157135337624e-2;
    pointWeightFem[13] = 0.709100346284691107301157135337624e-2;

    // pointWeightFem[0] = -74. / 5625.;
    // pointWeightFem[1] = 343. / 45000.;
    // pointWeightFem[2] = 343. / 45000.;
    // pointWeightFem[3] = 343. / 45000.;
    // pointWeightFem[4] = 343. / 45000.;
    // pointWeightFem[5] = 56. / 2250.;
    // pointWeightFem[6] = 56. / 2250.; 
    // pointWeightFem[7] = 56. / 2250.;
    // pointWeightFem[8] = 56. / 2250.; 
    // pointWeightFem[9] = 56. / 2250.;
    // pointWeightFem[10] = 56. / 2250.; 
    
    // pointWeightFem[0] = 0.1817020685825351/6.;
    // pointWeightFem[1] = 0.0361607142857143/6.;
    // pointWeightFem[2] = 0.0361607142857143/6.;
    // pointWeightFem[3] = 0.0361607142857143/6.;
    // pointWeightFem[4] = 0.0361607142857143/6.;
    // pointWeightFem[5] = 0.0698714945161738/6.;
    // pointWeightFem[6] = 0.0698714945161738/6.; 
    // pointWeightFem[7] = 0.0698714945161738/6.;
    // pointWeightFem[8] = 0.0698714945161738/6.; 
    // pointWeightFem[9] = 0.0656948493683187/6.;
    // pointWeightFem[10] = 0.0656948493683187/6.; 
    // pointWeightFem[11] = 0.0656948493683187/6.;
    // pointWeightFem[12] = 0.0656948493683187/6.;
    // pointWeightFem[13] = 0.0656948493683187/6.;
    // pointWeightFem[14] = 0.0656948493683187/6.;

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

template<>
double IntegQuadrature<3>::WeightListIso(int i){

// pointWeightIso[0] = 0.27735296695391276067;
// pointWeightIso[1] = 0.27735296695391276067;
// pointWeightIso[2] = 0.27735296695391276067;
// pointWeightIso[3] = 0.27735296695391276067;
// pointWeightIso[4] = 0.27735296695391276067;
// pointWeightIso[5] = 0.27735296695391276067;
// pointWeightIso[6] = 0.27735296695391276067;
// pointWeightIso[7] = 0.27735296695391276067;
// pointWeightIso[8] = 0.27735296695391276067;
// pointWeightIso[9] = 0.27735296695391276067;
// pointWeightIso[10] = 0.27735296695391276067;
// pointWeightIso[11] = 0.27735296695391276067;
// pointWeightIso[12] = 0.27735296695391276067;
// pointWeightIso[13] = 0.27735296695391276067;
// pointWeightIso[14] = 0.27735296695391276067;
// pointWeightIso[15] = 0.27735296695391276067;
// pointWeightIso[16] = 0.27735296695391276067;
// pointWeightIso[17] = 0.27735296695391276067;
// pointWeightIso[18] = 0.27735296695391276067;
// pointWeightIso[19] = 0.27735296695391276067;
// pointWeightIso[20] = 0.27735296695391276067;
// pointWeightIso[21] = 0.27735296695391276067;
// pointWeightIso[22] = 0.27735296695391276067;
// pointWeightIso[23] = 0.27735296695391276067;
// pointWeightIso[24] = 0.27735296695391276067;
// pointWeightIso[25] = 0.27735296695391276067;
// pointWeightIso[26] = 0.27735296695391276067;
// pointWeightIso[27] = 0.27735296695391276067;
// pointWeightIso[28] = 0.27735296695391276067;
// pointWeightIso[29] = 0.27735296695391276067;
// pointWeightIso[30] = 0.27735296695391276067;
// pointWeightIso[31] = 0.27735296695391276067;
// pointWeightIso[32] = 0.14794033605678127974;
// pointWeightIso[33] = 0.14794033605678127974;
// pointWeightIso[34] = 0.14794033605678127974;
// pointWeightIso[35] = 0.14794033605678127974;
// pointWeightIso[36] = 0.14794033605678127974;
// pointWeightIso[37] = 0.14794033605678127974;
// pointWeightIso[38] = 0.14794033605678127974;
// pointWeightIso[39] = 0.14794033605678127974;
// pointWeightIso[40] = 0.14794033605678127974;
// pointWeightIso[41] = 0.14794033605678127974;
// pointWeightIso[42] = 0.14794033605678127974;
// pointWeightIso[43] = 0.14794033605678127974;
// pointWeightIso[44] = 0.14794033605678127974;
// pointWeightIso[45] = 0.14794033605678127974;
// pointWeightIso[46] = 0.14794033605678127974;
// pointWeightIso[47] = 0.14794033605678127974;
// pointWeightIso[48] = 0.14794033605678127974;
// pointWeightIso[49] = 0.14794033605678127974;
// pointWeightIso[50] = 0.14794033605678127974;
// pointWeightIso[51] = 0.14794033605678127974;
// pointWeightIso[52] = 0.14794033605678127974;
// pointWeightIso[53] = 0.14794033605678127974;
// pointWeightIso[54] = 0.14794033605678127974;
// pointWeightIso[55] = 0.14794033605678127974;
// pointWeightIso[56] = 0.14794033605678127974;
// pointWeightIso[57] = 0.14794033605678127974;
// pointWeightIso[58] = 0.14794033605678127974;
// pointWeightIso[59] = 0.14794033605678127974;
// pointWeightIso[60] = 0.14794033605678127974;
// pointWeightIso[61] = 0.14794033605678127974;
// pointWeightIso[62] = 0.14794033605678127974;
// pointWeightIso[63] = 0.14794033605678127974;
    
    pointWeightIso[0] = 0.702331961591221;
    pointWeightIso[1] = 0.438957475994514;
    pointWeightIso[2] = 0.438957475994514;
    pointWeightIso[3] = 0.438957475994514;
    pointWeightIso[4] = 0.274348422496571;
    pointWeightIso[5] = 0.274348422496571;
    pointWeightIso[6] = 0.438957475994514;
    pointWeightIso[7] = 0.274348422496571;
    pointWeightIso[8] = 0.274348422496571;
    pointWeightIso[9] = 0.438957475994514;
    pointWeightIso[10] = 0.274348422496571;
    pointWeightIso[11] = 0.274348422496571;
    pointWeightIso[12] = 0.274348422496571;
    pointWeightIso[13] = 0.171467764060357;
    pointWeightIso[14] = 0.171467764060357;
    pointWeightIso[15] = 0.274348422496571;
    pointWeightIso[16] = 0.171467764060357;
    pointWeightIso[17] = 0.171467764060357;
    pointWeightIso[18] = 0.438957475994514;
    pointWeightIso[19] = 0.274348422496571;
    pointWeightIso[20] = 0.274348422496571;
    pointWeightIso[21] = 0.274348422496571;
    pointWeightIso[22] = 0.171467764060357;
    pointWeightIso[23] = 0.171467764060357;
    pointWeightIso[24] = 0.274348422496571;
    pointWeightIso[25] = 0.171467764060357;
    pointWeightIso[26] = 0.171467764060357;

    return pointWeightIso[i];
};

//------------------------------------------------------------------------------
//-----------COMPUTES THE VALUE INTERPOLATED IN THE INTEGRATION POINT-----------
//------------------------------------------------------------------------------
template<int DIM>
double IntegQuadrature<DIM>::interpolateQuadraticVariableFem(double *nValues, int &point) {
    
    QuadShapeFunction<DIM>    shapeQuad;

    double xsi[DIM];
    for (int i = 0; i < DIM; i++) xsi[i] = PointListFem(point,i);
    
    //local number of nodes
    int LNN = 4*DIM-2;
    double phi_[LNN];
    shapeQuad.evaluateFem(xsi,phi_);
    
    double int_value = 0.;
    for (int i = 0; i < LNN; i++){
        int_value += nValues[i] * phi_[i];
    };

    return int_value;
};


template<int DIM>
double IntegQuadrature<DIM>::interpolateQuadraticVariableIso(double *nValues, int &point, double *wpc, int *INC_, std::vector<IParameters *> &iparam, int &patch) {
    
    iparameters = &iparam;
    QuadShapeFunction<DIM>    shapeQuad;

    double xsi[DIM];
    for (int i = 0; i < DIM; i++) xsi[i] = PointListIso(point,i);

    int LNN = 18*DIM -27;

    double phi_[LNN];
    shapeQuad.evaluateIso(xsi,phi_,wpc,INC_,(*iparameters),patch);
    
    double int_value = 0.;
    for (int i = 0; i < LNN; i++){
        int_value += nValues[i] * phi_[i];
    };

    return int_value;
};

template class IntegQuadrature<2>;
template class IntegQuadrature<3>;

