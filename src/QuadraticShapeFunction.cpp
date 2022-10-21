#include "QuadraticShapeFunction.h"

//------------------------------------------------------------------------------
//-------------------------COMPUTE SHAPE FUNCTION VALUE-------------------------
//------------------------------------------------------------------------------
// Defines quadratic shape functions and its derivatives 
// for triangles and tetrahedrons
template<>
void QuadShapeFunction<2>::evaluateFem(double *xsi, double *phi) const {
    
     const double xsi1 = xsi[0];
     const double xsi2 = xsi[1];
     const double xsi3 = 1. - xsi1 - xsi2;

     phi[0] = xsi3 * (2.0 * xsi3 - 1.0);
     phi[1] = xsi1 * (2.0 * xsi1 - 1.0);
     phi[2] = xsi2 * (2.0 * xsi2 - 1.0);
     phi[3] = 4.0 * xsi3 * xsi1;
     phi[4] = 4.0 * xsi1 * xsi2;
     phi[5] = 4.0 * xsi2 * xsi3;

     return;
}


template<>
void QuadShapeFunction<2>::evaluateIso(double *xsi, double *phi, double *wpcs, int *inc, std::vector<IParameters_ *> &iparameters ,int Npatch) {
    
    std::vector<IParameters_ *> *ipar_; // isogeometric parameters    
    ipar_ = &iparameters;

    double Xsi[2];

    Xsi[0] = xsi[0];
    Xsi[1] = xsi[1];
       
    int degm = (*ipar_)[Npatch] -> getDegree(0);
    int degn = (*ipar_)[Npatch] -> getDegree(1);
    int npcm = (*ipar_)[Npatch] -> getNcp(0);
    int npcn = (*ipar_)[Npatch] -> getNcp(1);   

    int numLocalBF = (degm+1)*(degn+1);

    double *uknot;
    double *vknot;
    double  uLeft [3];
    double  uRight[3];
    double  vLeft [3];
    double  vRight[3];
    double  phiL1 [3];
    double  phiL2 [3];
    double  phit  [9];
    double  uBF[3][3] = {}; // stores de base functions in each direction
    double  vBF[3][3] = {}; // stores de base functions in each direction
    double saved,temp;

    uknot = (*ipar_)[Npatch]->getuKnot();
    vknot = (*ipar_)[Npatch]->getvKnot();  

    // indexes from index space where a function related with first local point of element start 
    int uind = inc[0];
    int vind = inc[1]; 

    // parametric coordinates that defines the element
    double u1 = uknot[uind];
    double u2 = uknot[uind+1];
    double v1 = vknot[vind];
    double v2 = vknot[vind +1];

    // relates integration space with the parametric one
    const double xsi1 =((u2 - u1) * Xsi[0] + (u2 + u1)) * 0.5;
    const double xsi2 =((v2 - v1) * Xsi[1] + (v2 + v1)) * 0.5;

    //Functions in u direction
    uBF[0][0] = 1.0;

    for (int j = 1; j < (degm+1); j++) {
        uLeft[j] = xsi1 - uknot[uind + 1 - j];
        uRight[j] = uknot[uind + j] - xsi1;
        saved = 0.;
        for (int r = 0; r<j; r++) {
            //lower triangle (Piegl and Tiller pg. 68-71)
            uBF[j][r] = uRight[r+1] + uLeft[j-r];
            temp = uBF[r][j-1]/uBF[j][r];
            //upper triangle
            uBF[r][j] = saved + uRight[r+1]*temp;
            saved = uLeft[j-r]*temp;
        };
        uBF[j][j] = saved;
    };

    for (int i = 0; i < degm +1 ; i++ ) {
        phiL1[i] = uBF[i][degm];
    };

    //Functions in v direction
    vBF[0][0] = 1.0;

    for (int j = 1; j < (degn +1); j++) {
        vLeft[j] = xsi2 - vknot[vind + 1 - j];
        vRight[j] = vknot[vind + j] - xsi2;
        saved = 0.;
        for (int r = 0; r < j; r++) {
            //lower triangle (Piegl and Tiller pg. 68-71)
            vBF[j][r] = vRight[r+1] + vLeft[j-r];
            temp = vBF[r][j-1]/vBF[j][r];
            //upper triangle
            vBF[r][j] = saved + vRight[r+1]*temp;
            saved = vLeft[j-r]*temp;
        };
        vBF[j][j] = saved;
    };

    for (int i = 0; i < degn +1 ; i++ ) {
        phiL2[i] = vBF[i][degn];
    };
    

    int index = 0;
    for (int j = 0; j < degn + 1 ; j ++) {
        for (int k = 0; k < degm + 1; k ++) {
            phit[index] = phiL1[k] * phiL2 [j];
            index ++;
        };
    };

    double sum = 0.0;
    for (int i = 0; i < numLocalBF; i++) {
        sum += phit[i] * wpcs[i];
    };

    for (int i = 0; i < numLocalBF; i++) {
        phi[i] = phit[i] * wpcs[i]/sum;
    };

    return;
};


template<>
void QuadShapeFunction<3>::evaluateFem(double *xsi, double *phi) const {
    
     const double xsi1 = xsi[0];
     const double xsi2 = xsi[1];
     const double xsi3 = xsi[2];

     phi[0] = 2.0 * (xsi1 - 0.50) * xsi1; 
     phi[1] = 2.0 * (xsi2 - 0.50) * xsi2; 
     phi[2] = 2.0 * (xsi3 - 0.50) * xsi3; 
     phi[3] = (2.0 - 2.0 * xsi3 - 2.0 * xsi2 - 2.0 * xsi1 - 1.0)  \
              * (1.0 - xsi1 - xsi2 - xsi3); 
     phi[4] = 4.0 * xsi1 * xsi2; 
     phi[5] = 4.0 * xsi2 * xsi3;
     phi[6] = 4.0 * xsi1 * xsi3;
     phi[7] = 4.0 * xsi1 * (1.0 - xsi1 - xsi2 - xsi3);
     phi[8] = 4.0 * xsi3 * (1.0 - xsi1 - xsi2 - xsi3);
     phi[9] = 4.0 * xsi2 * (1.0 - xsi1 - xsi2 - xsi3);
     
     return;
}

template<>
void QuadShapeFunction<3>::evaluateIso(double *xsi, double *phi, double *wpcs, int *inc, std::vector<IParameters_ *> &iparameters ,int Npatch) {
    
    double Xsi[3];
    Xsi[0] = xsi[0];
    Xsi[1] = xsi[1];
    Xsi[2] = xsi[2];

    std::vector<IParameters_ *> *ipar_; // isogeometric parameters    
    ipar_ = &iparameters;
       
    int degm = (*ipar_)[Npatch] -> getDegree(0);
    int degn = (*ipar_)[Npatch] -> getDegree(1);
    int degq = (*ipar_)[Npatch] -> getDegree(2);
    int npcm = (*ipar_)[Npatch] -> getNcp(0);
    int npcn = (*ipar_)[Npatch] -> getNcp(1);   
    int npcq = (*ipar_)[Npatch] -> getNcp(2); 


    int numLocalBF = (degm+1)*(degn+1)*(degq+1);

    double *uknot;
    double *vknot;
    double *tknot;
    double uLeft [3];
    double uRight[3];
    double vLeft [3];
    double vRight[3];
    double tLeft [3];
    double tRight[3];
    double phiL1 [3];
    double phiL2 [3];
    double phiL3 [3];
    double phit  [27];
    double uBF[3][3] = {}; // stores de base functions in each direction
    double vBF[3][3] = {};  // stores de base functions in each direction
    double tBF[3][3] = {};  // stores de base functions in each direction
    double saved,temp;

    uknot = (*ipar_)[Npatch]->getuKnot();
    vknot = (*ipar_)[Npatch]->getvKnot(); 
    tknot = (*ipar_)[Npatch]->gettKnot(); 

    // indexes from index space where a function related with first local point of element start 
    int uind = inc[0];
    int vind = inc[1];
    int tind = inc[2]; 

    // parametric coordinates that defines the element
    double u1 = uknot[uind];
    double u2 = uknot[uind+1];
    double v1 = vknot[vind];
    double v2 = vknot[vind +1];
    double t1 = tknot[tind];
    double t2 = tknot[tind +1];

    // relates integration space with the parametric one
    const double xsi1 =((u2 - u1) * Xsi[0] + (u2 + u1)) * 0.5;
    const double xsi2 =((v2 - v1) * Xsi[1] + (v2 + v1)) * 0.5;
    const double xsi3 =((t2 - t1) * Xsi[2] + (t2 + t1)) * 0.5;

    //Functions in parametric direction xsi 
    uBF[0][0] = 1.0;

    for (int j = 1; j < (degm+1); j++) {
        uLeft[j] = xsi1 - uknot[uind + 1 - j];
        uRight[j] = uknot[uind + j] - xsi1;
        saved = 0.;
        for (int r = 0; r<j; r++) {
            //lower triangle (Piegl and Tiller pg. 68-71)
            uBF[j][r] = uRight[r+1] + uLeft[j-r];
            temp = uBF[r][j-1]/uBF[j][r];
            //upper triangle
            uBF[r][j] = saved + uRight[r+1]*temp;
            saved = uLeft[j-r]*temp;
        };
        uBF[j][j] = saved;
    };

    for (int i = 0; i < degm +1 ; i++ ) {
        phiL1[i] = uBF[i][degm];
    };

    //Functions in parametric direction eta 
    vBF[0][0] = 1.0;

    for (int j = 1; j < (degn +1); j++) {
        vLeft[j] = xsi2 - vknot[vind + 1 - j];
        vRight[j] = vknot[vind + j] - xsi2;
        saved = 0.;
        for (int r = 0; r < j; r++) {
            //lower triangle (Piegl and Tiller pg. 68-71)
            vBF[j][r] = vRight[r+1] + vLeft[j-r];
            temp = vBF[r][j-1]/vBF[j][r];
            //upper triangle
            vBF[r][j] = saved + vRight[r+1]*temp;
            saved = vLeft[j-r]*temp;
        };
        vBF[j][j] = saved;
    };

    for (int i = 0; i < degn +1 ; i++ ) {
        phiL2[i] = vBF[i][degn];
    };

    //Functions in parametric direction zeta 
    tBF[0][0] = 1.0;

    for (int j = 1; j < (degq +1); j++) {
        tLeft[j] = xsi3 - tknot[tind + 1 - j];
        tRight[j] = tknot[tind + j] - xsi3;
        saved = 0.;
        for (int r = 0; r < j; r++) {
            //lower triangle (Piegl and Tiller pg. 68-71)
            tBF[j][r] = tRight[r+1] + tLeft[j-r];
            temp = tBF[r][j-1]/tBF[j][r];
            //upper triangle
            tBF[r][j] = saved + tRight[r+1]*temp;
            saved = tLeft[j-r]*temp;
        };
        tBF[j][j] = saved;
    };

    for (int i = 0; i < degq +1 ; i++ ) {
        phiL3[i] = tBF[i][degq];
    };
    

    int index = 0;
    for (int i = 0; i < degq + 1; i++){
        for (int j = 0; j < degn + 1 ; j ++) {
            for (int k = 0; k < degm + 1; k ++) {
                phit[index] = phiL1[k] * phiL2 [j] * phiL3[i];
                index ++;
            };
        };
    };
        
    double sum = 0.0;
    for (int i = 0; i < numLocalBF; i++) {
        sum += phit[i] * wpcs[i];
    };

    for (int i = 0; i < numLocalBF; i++) {
        phi[i] = phit[i] * wpcs[i]/sum;
    };

    return;
};


//------------------------------------------------------------------------------
//-------------------COMPUTE SHAPE FUNCTION DERIVATIVE VALUE--------------------
//------------------------------------------------------------------------------
template<>
void QuadShapeFunction<2>::evaluateGradientFem(double *xsi, double **dphi) const {

     const double xsi1 = xsi[0];
     const double xsi2 = xsi[1];
     const double xsi3 = 1. - xsi1 - xsi2;

     dphi[0][0] = -4. * xsi3 + 1.;
     dphi[1][0] = -4. * xsi3 + 1.;

     dphi[0][1] = 4. * xsi1 - 1.;
     dphi[1][1] = 0.;

     dphi[0][2] = 0.;
     dphi[1][2] = 4. * xsi2 - 1.;
 
     dphi[0][3] = 4. * (xsi3 - xsi1);
     dphi[1][3] = -4. * xsi1;

     dphi[0][4] = 4. * xsi2;
     dphi[1][4] = 4. * xsi1;

     dphi[0][5] = -4. * xsi2;
     dphi[1][5] = 4. * (xsi3 - xsi2);

     
    return;
}

template<>
void QuadShapeFunction<2>::evaluateGradientIso(double *xsi, double **dphi,double *wpcs, int *inc, std::vector<IParameters_ *> &iparameters ,int Npatch) {

    std::vector<IParameters_ *> *ipar_;  // isogeometric parameters

    ipar_ = &iparameters;

    double Xsi[2];
    Xsi[0] = xsi[0];
    Xsi[1] = xsi[1];
       
    int degm = (*ipar_)[Npatch] -> getDegree(0);
    int degn = (*ipar_)[Npatch] -> getDegree(1);
    int npcm = (*ipar_)[Npatch] -> getNcp(0);
    int npcn = (*ipar_)[Npatch] -> getNcp(1);  
 
    int numLocalBF = (degm+1)*(degn+1);

    double *uknot;
    double *vknot;
    double uLeft [3];
    double uRight[3];
    double vLeft [3];
    double vRight[3];
    double phiL1 [3];
    double phiL2 [3];
    double phit  [9];
    double dphiL1[3];;
    double dphiL2[3];;
    double dphit [2][9];
    double aux1[2][2]; // stores the two lines more recently computed
    double aux2[2][2]; // stores the two lines more recently computed
    double uBF[3][3] = {}; // stores de base functions in each direction
    double vBF[3][3] = {}; // stores de base functions in each direction
    double saved,temp;

    uknot = (*ipar_)[Npatch] ->  getuKnot();
    vknot = (*ipar_)[Npatch] ->  getvKnot();  

    // indexes from index space where a function related with first local point of element start 
    int uind = inc [0];
    int vind = inc [1]; 

    // parametric coordinates that defines the element
    double u1 = uknot[uind];
    double u2 = uknot[uind+1];
    double v1 = vknot[vind];
    double v2 = vknot[vind +1];

    // relates integration space with the parametric one
    const double xsi1 =((u2 - u1) * Xsi[0] + (u2 + u1)) * 0.5;
    const double xsi2 =((v2 - v1) * Xsi[1] + (v2 + v1)) * 0.5;

    // Base function u direction
    uBF[0][0] = 1.0;
    for (int j = 1; j < (degm + 1); j++) {
        uLeft[j] = xsi1 - uknot[uind + 1 - j];
        uRight[j] = uknot[uind + j] - xsi1;
        saved = 0.;
        for (int r = 0; r < j; r++) {
            //lower triangle (Piegl and Tiller pg. 68-71)
            uBF[j][r] = uRight[r + 1] + uLeft[j - r];
            temp = uBF[r][j - 1]/uBF[j][r];
            //upper triangle
            uBF[r][j] = saved + uRight[r + 1]*temp;
            saved = uLeft[j - r] * temp;
        };
        uBF[j][j] = saved;
    };

    for (int i = 0; i < (degm +1); i++ ) {
        phiL1[i] = uBF[i][degm];
    };

    // Base function v direction

    vBF[0][0] = 1.0;
    for (int j = 1; j< (degn + 1); j++) {
        vLeft[j] = xsi2 - vknot[vind + 1 - j];
        vRight[j] = vknot[vind + j] - xsi2;
        saved = 0.;
        for (int r = 0; r < j; r++) {
            //lower triangle (Piegl and Tiller pg. 68-71)
            vBF[j][r] = vRight[r + 1] + vLeft[j - r];
            temp = vBF[r][j - 1]/vBF[j][r];
            //upper triangle
            vBF[r][j] = saved + vRight[r + 1]*temp;
            saved = vLeft[j - r]*temp;
        };
        vBF[j][j] = saved;
    };

    for (int i = 0; i < (degn + 1) ; i++ ) {
        phiL2[i] = vBF[i][degn];
    };

    int index = 0;
    for (int j = 0; j< (degn + 1) ; j ++) {
        for (int k = 0; k < (degm + 1); k ++) {
            phit[index] = phiL1[k] * phiL2 [j];
            index ++;
        };
    };

    // derivatives u direction;
    int s1,s2,rk,pk,j1,j2,cor,k;
    double d;

    for (int r = 0; r < (degm + 1); r++) {
        s1 = 0;
        s2 = 1;
        aux1 [0][0] = 1.0;
        k = 1;
        d = 0.0;
        rk = r - k;
        pk = degm - k;
        if (r >= k) {
            aux1[s2][0] = aux1[s1][0] / uBF[pk + 1][rk];
            d = aux1[s2][0] * uBF[rk][pk];
        };
        if (rk >= -1){
            j1 = 1;
        } else {
            j1 = -rk;
        };
        if ((r-1) <= pk) {
            j2 = k - 1 ;
        } else {
            j2 = degm - r;
        };
        for (int j = j1; j <= j2; j++) {
            aux1[s2][j] = (aux1[s1][j] - aux1[s1][j - 1]) / uBF[pk + 1][rk + j];
            d = d + aux1[s2][j] * uBF[rk + j][pk];
        };
        if (r <= pk){
            aux1[s2][k] = -aux1[s1][k - 1] / uBF[pk + 1][r];
            d = d + aux1[s2][k] * uBF[r][pk];
        };
        dphiL1[r] = d;
        int j = s1;
        s1 = s2;
        s2 = j;
    };

    for (int i = 0; i < (degm + 1); i++) {
        dphiL1[i] = dphiL1[i] * degm;
    };

    // derivatives v direction;
    for (int r = 0; r < (degn + 1); r++) {
        s1 = 0;
        s2 = 1;
        aux2 [0][0] = 1.0;
        k = 1;
        d = 0.0;
        rk = r - k;
        pk = degn - k;
        if (r >= k) {
            aux2[s2][0] = aux2[s1][0] / vBF[pk + 1][rk];
            d = aux2[s2][0] * vBF[rk][pk];
        };
        if (rk >= -1){
            j1 = 1;
        } else{
            j1 = -rk;
        };
        if ((r-1) <= pk) {
            j2 = k - 1 ;
        } else {
            j2 = degn - r;
        };
        for (int j = j1; j <= j2; j++) {
            aux2[s2][j] = (aux2[s1][j] - aux2[s1][j - 1]) / vBF[pk + 1][rk + j];
            d = d + aux2[s2][j] * vBF[rk + j][pk];
        };
        if (r <= pk){
            aux2[s2][k] = -aux2[s1][k - 1] / vBF[pk + 1][r];
            d = d + aux2[s2][k] * vBF[r][pk];
        };
        dphiL2 [r] = d;
        int j = s1;
        s1 = s2;
        s2 = j;
    };

    for (int i = 0; i < (degn + 1); i++) {
        dphiL2[i] = dphiL2[i] * degn;
    };


    index = 0;
    for (int j = 0; j< (degn + 1) ; j ++) {
        for (int k = 0; k < (degm + 1); k ++) {
            dphit[0][index] = dphiL1[k] * phiL2 [j];
            dphit[1][index] = phiL1[k] * dphiL2 [j];
            index ++;
        };
    };

    double sum = 0.0;
    double sum1 = 0.0;
    double sum2 = 0.0;
    
    for (int i = 0; i < numLocalBF; i++) {
        sum += phit[i] * wpcs[i];
        sum1 += dphit[0][i] * wpcs[i];
        sum2 += dphit[1][i] * wpcs[i];
    };

    for (int i = 0; i < numLocalBF; i++) {
        dphi[0][i] = ((dphit[0][i] * wpcs[i] * sum - phit[i] * wpcs[i] * sum1) / (sum * sum)); 
        dphi[1][i] = ((dphit[1][i] * wpcs[i] * sum - phit[i] * wpcs[i] * sum2) / (sum * sum)); 
    };

    return;
};
 

template<>
void QuadShapeFunction<3>::evaluateGradientFem(double *xsi, double **dphi) const {

     const double xsi1 = xsi[0];
     const double xsi2 = xsi[1];
     const double xsi3 = xsi[2];

     dphi[0][0] = 4. * xsi1 - 1.;
     dphi[1][0] = 0.;
     dphi[2][0] = 0.;

     dphi[0][1] = 0.;
     dphi[1][1] = 4. * xsi2 - 1.;
     dphi[2][1] = 0.;
     
     dphi[0][2] = 0.;
     dphi[1][2] = 0.;
     dphi[2][2] = 4. * xsi3 - 1.;

     dphi[0][3] = 4. * (xsi1 + xsi2 + xsi3) - 3.;
     dphi[1][3] = 4. * (xsi1 + xsi2 + xsi3) - 3.;
     dphi[2][3] = 4. * (xsi1 + xsi2 + xsi3) - 3.;

     dphi[0][4] = 4. * xsi2;
     dphi[1][4] = 4. * xsi1;
     dphi[2][4] = 0.;

     dphi[0][5] = 0.;
     dphi[1][5] = 4. * xsi3;
     dphi[2][5] = 4. * xsi2;

     dphi[0][6] = 4. * xsi3;
     dphi[1][6] = 0.;
     dphi[2][6] = 4. * xsi1;

     dphi[0][7] = 4. * (1. - 2. * xsi1 - xsi2 - xsi3);
     dphi[1][7] = -4. * xsi1;
     dphi[2][7] = -4. * xsi1;

     dphi[1][8] = -4. * xsi3;
     dphi[0][8] = -4. * xsi3;
     dphi[2][8] = 4. * (1. - 2. * xsi3 - xsi2 - xsi1);

     dphi[0][9] = -4. * xsi2;
     dphi[1][9] = 4. * (1. - 2. * xsi2 - xsi1 - xsi3);
     dphi[2][9] = -4. * xsi2;

     
    return;
}

template<>
void QuadShapeFunction<3>::evaluateGradientIso(double *xsi, double **dphi, double *wpcs, int *inc, std::vector<IParameters_ *> &iparameters ,int Npatch) {

    double Xsi[3];
    Xsi[0] = xsi[0];
    Xsi[1] = xsi[1];
    Xsi[2] = xsi[2];

    std::vector<IParameters_ *> *ipar_;  // isogeometric parameters
    ipar_ = &iparameters;
       
    int degm = (*ipar_)[Npatch] -> getDegree(0);
    int degn = (*ipar_)[Npatch] -> getDegree(1);
    int degq = (*ipar_)[Npatch] -> getDegree(2);
    int npcm = (*ipar_)[Npatch] -> getNcp(0);
    int npcn = (*ipar_)[Npatch] -> getNcp(1);  
    int npcq = (*ipar_)[Npatch] -> getNcp(2);
 
    int numLocalBF = (degm+1)*(degn+1)*(degq+1);

    double *uknot;
    double *vknot;
    double *tknot;
    double  uLeft [3];
    double uRight[3];
    double vLeft [3];
    double vRight[3];
    double tLeft [3];
    double tRight[3];
    double phiL1 [3];
    double phiL2 [3];
    double phiL3 [3];
    double phit  [27];
    double dphiL1[3];
    double dphiL2[3];
    double dphiL3[3];
    double dphit [3][27];
    double aux1[2][2]; // stores the two lines more recently computed
    double aux2[2][2]; // stores the two lines more recently computed
    double aux3[2][2]; // stores the two lines more recently computed
    double uBF[3][3] = {};// stores de base functions in each direction
    double vBF[3][3] = {}; // stores de base functions in each direction
    double tBF[3][3] = {}; // stores de base functions in each direction
    double saved,temp;

    uknot = (*ipar_)[Npatch] ->  getuKnot();
    vknot = (*ipar_)[Npatch] ->  getvKnot();  
    tknot = (*ipar_)[Npatch] ->  gettKnot();

    // indexes from index space where a function related with first local point of element start 
    int uind = inc [0];
    int vind = inc [1]; 
    int tind = inc [2];

    // parametric coordinates that defines the element
    double u1 = uknot[uind];
    double u2 = uknot[uind+1];
    double v1 = vknot[vind];
    double v2 = vknot[vind +1];
    double t1 = tknot[tind];
    double t2 = tknot[tind +1];

    // relates integration space with the parametric one
    const double xsi1 =((u2 - u1) * Xsi[0] + (u2 + u1)) * 0.5;
    const double xsi2 =((v2 - v1) * Xsi[1] + (v2 + v1)) * 0.5;
    const double xsi3 =((t2 - t1) * Xsi[2] + (t2 + t1)) * 0.5;
 
    // Base function u direction
    uBF[0][0] = 1.0;
    for (int j = 1; j < (degm + 1); j++) {
        uLeft[j] = xsi1 - uknot[uind + 1 - j];
        uRight[j] = uknot[uind + j] - xsi1;
        saved = 0.;
        for (int r = 0; r < j; r++) {
            //lower triangle (Piegl and Tiller pg. 68-71)
            uBF[j][r] = uRight[r + 1] + uLeft[j - r];
            temp = uBF[r][j - 1]/uBF[j][r];
            //upper triangle
            uBF[r][j] = saved + uRight[r + 1]*temp;
            saved = uLeft[j - r] * temp;
        };
        uBF[j][j] = saved;
    };

    for (int i = 0; i < (degm +1); i++ ) {
        phiL1[i] = uBF[i][degm];
    };

    // Base function v direction
    vBF[0][0] = 1.0;
    for (int j = 1; j< (degn + 1); j++) {
        vLeft[j] = xsi2 - vknot[vind + 1 - j];
        vRight[j] = vknot[vind + j] - xsi2;
        saved = 0.;
        for (int r = 0; r < j; r++) {
            //lower triangle (Piegl and Tiller pg. 68-71)
            vBF[j][r] = vRight[r + 1] + vLeft[j - r];
            temp = vBF[r][j - 1]/vBF[j][r];
            //upper triangle
            vBF[r][j] = saved + vRight[r + 1]*temp;
            saved = vLeft[j - r]*temp;
        }
        vBF[j][j] = saved;
    };

    for (int i = 0; i < (degn + 1) ; i++ ) {
        phiL2[i] = vBF[i][degn];
    };

    // Base function t direction
    tBF[0][0] = 1.0;
    for (int j = 1; j< (degq + 1); j++) {
        tLeft[j] = xsi3 - tknot[tind + 1 - j];
        tRight[j] = tknot[tind + j] - xsi3;
        saved = 0.;
        for (int r = 0; r < j; r++) {
            //lower triangle (Piegl and Tiller pg. 68-71)
            tBF[j][r] = tRight[r + 1] + tLeft[j - r];
            temp = tBF[r][j - 1]/tBF[j][r];
            //upper triangle
            tBF[r][j] = saved + tRight[r + 1]*temp;
            saved = tLeft[j - r]*temp;
        };
        tBF[j][j] = saved;
    };

    for (int i = 0; i < (degq + 1) ; i++ ) {
        phiL3[i] = tBF[i][degq];
    };

    int index = 0;

    for (int i = 0; i< (degq + 1) ; i ++) {
        for (int j = 0; j< (degn + 1) ; j ++) {
            for (int k = 0; k < (degm + 1); k ++) {
                phit[index] = phiL1[k] * phiL2 [j] * phiL3 [i];
                index ++;
            };
        };
    };
        

    // derivatives u direction;
    int s1,s2,rk,pk,j1,j2,cor,k;
    double d;

    for (int r = 0; r < (degm + 1); r++) {
        s1 = 0;
        s2 = 1;
        aux1 [0][0] = 1.0;
        k = 1;
        d = 0.0;
        rk = r - k;
        pk = degm - k;
        if (r >= k) {
            aux1[s2][0] = aux1[s1][0] / uBF[pk + 1][rk];
            d = aux1[s2][0] * uBF[rk][pk];
        };
        if (rk >= -1){
            j1 = 1;
        } else {
            j1 = -rk;
        };
        if ((r-1) <= pk) {
            j2 = k - 1 ;
        } else {
            j2 = degm - r;
        };
        for (int j = j1; j <= j2; j++) {
            aux1[s2][j] = (aux1[s1][j] - aux1[s1][j - 1]) / uBF[pk + 1][rk + j];
            d = d + aux1[s2][j] * uBF[rk + j][pk];
        };
        if (r <= pk){
            aux1[s2][k] = -aux1[s1][k - 1] / uBF[pk + 1][r];
            d = d + aux1[s2][k] * uBF[r][pk];
        };
        dphiL1[r] = d;
        int j = s1;
        s1 = s2;
        s2 = j;
    };

    for (int i = 0; i < (degm + 1); i++) {
        dphiL1[i] = dphiL1[i] * degm;
    };

    // derivatives v direction;
    for (int r = 0; r < (degn + 1); r++) {
        s1 = 0;
        s2 = 1;
        aux2 [0][0] = 1.0;
        k = 1;
        d = 0.0;
        rk = r - k;
        pk = degn - k;
        if (r >= k) {
            aux2[s2][0] = aux2[s1][0] / vBF[pk + 1][rk];
            d = aux2[s2][0] * vBF[rk][pk];
        };
        if (rk >= -1){
            j1 = 1;
        } else{
            j1 = -rk;
        };
        if ((r-1) <= pk) {
            j2 = k - 1 ;
        } else {
            j2 = degn - r;
        };
        for (int j = j1; j <= j2; j++) {
            aux2[s2][j] = (aux2[s1][j] - aux2[s1][j - 1]) / vBF[pk + 1][rk + j];
            d = d + aux2[s2][j] * vBF[rk + j][pk];
        };
        if (r <= pk){
            aux2[s2][k] = -aux2[s1][k - 1] / vBF[pk + 1][r];
            d = d + aux2[s2][k] * vBF[r][pk];
        };
        dphiL2 [r] = d;
        int j = s1;
        s1 = s2;
        s2 = j;
    };

    for (int i = 0; i < (degn + 1); i++) {
        dphiL2[i] = dphiL2[i] * degn;
    };

    // derivatives t direction;
    for (int r = 0; r < (degq + 1); r++) {
        s1 = 0;
        s2 = 1;
        aux3 [0][0] = 1.0;
        k = 1;
        d = 0.0;
        rk = r - k;
        pk = degq - k;
        if (r >= k) {
            aux3[s2][0] = aux3[s1][0] / tBF[pk + 1][rk];
            d = aux3[s2][0] * tBF[rk][pk];
        };
        if (rk >= -1){
            j1 = 1;
        } else{
            j1 = -rk;
        };
        if ((r-1) <= pk) {
            j2 = k - 1 ;
        } else {
            j2 = degq - r;
        };
        for (int j = j1; j <= j2; j++) {
            aux3[s2][j] = (aux3[s1][j] - aux3[s1][j - 1]) / tBF[pk + 1][rk + j];
            d = d + aux3[s2][j] * tBF[rk + j][pk];
        };
        if (r <= pk){
            aux3[s2][k] = -aux3[s1][k - 1] / tBF[pk + 1][r];
            d = d + aux3[s2][k] * tBF[r][pk];
        };
        dphiL3 [r] = d;
        int j = s1;
        s1 = s2;
        s2 = j;
    };

    for (int i = 0; i < (degq + 1); i++) {
        dphiL3[i] = dphiL3[i] * degq;
    };


    index = 0;
    for (int i = 0; i< (degq + 1) ; i ++) {
        for (int j = 0; j< (degn + 1) ; j ++) {
            for (int k = 0; k < (degm + 1); k ++) {
                dphit[0][index] = dphiL1[k] * phiL2 [j] * phiL3 [i];
                dphit[1][index] = phiL1[k] * dphiL2 [j] * phiL3 [i];
                dphit[2][index] = phiL1[k] * phiL2 [j] * dphiL3 [i];
                index ++;
            };
        };
    };
    

    double sum = 0.0;
    double sum1 = 0.0;
    double sum2 = 0.0;
    double sum3 = 0.0;
    
    for (int i = 0; i < numLocalBF; i++) {
        sum += phit[i] * wpcs[i];
        sum1 += dphit[0][i] * wpcs[i];
        sum2 += dphit[1][i] * wpcs[i];
        sum3 += dphit[2][i] * wpcs[i];
    };

    for (int i = 0; i < numLocalBF; i++) {
        dphi[0][i] = ((dphit[0][i] * wpcs[i] * sum - phit[i] * wpcs[i] * sum1) / (sum * sum)); 
        dphi[1][i] = ((dphit[1][i] * wpcs[i] * sum - phit[i] * wpcs[i] * sum2) / (sum * sum)); 
        dphi[2][i] = ((dphit[2][i] * wpcs[i] * sum - phit[i] * wpcs[i] * sum3) / (sum * sum)); 
    };

    return;
};

template<>
void QuadShapeFunction<2>::evaluateHessianFem(double *xsi, double ***ddphi){

     ddphi[0][0][0] = 4.;
     ddphi[0][1][0] = 4.;
     ddphi[1][0][0] = 4.;
     ddphi[1][1][0] = 4.;

     ddphi[0][0][1] = 4.;
     ddphi[0][1][1] = 0.;
     ddphi[1][0][1] = 0.;
     ddphi[1][1][1] = 0.;

     ddphi[0][0][2] = 0.;
     ddphi[0][1][2] = 0.;
     ddphi[1][0][2] = 0.;
     ddphi[1][1][2] = 4.;

     ddphi[0][0][3] = -8.;
     ddphi[0][1][3] = -4.;
     ddphi[1][0][3] = -4.;
     ddphi[1][1][3] = 0.;

     ddphi[0][0][4] = 0.;
     ddphi[0][1][4] = 4.;
     ddphi[1][0][4] = 4.;
     ddphi[1][1][4] = 0.;

     ddphi[0][0][5] = 0.;
     ddphi[0][1][5] = -4.;
     ddphi[1][0][5] = -4.;
     ddphi[1][1][5] = -8.;

};

template<>
void QuadShapeFunction<2>::evaluateHessianIso(double *xsi, double*** ddphi, double* wpcs, int *inc, std::vector<IParameters_ *> &iparameters ,int Npatch){


    std::vector<IParameters_ *> *ipar_;  // isogeometric parameters

    ipar_ = &iparameters;


    double Xsi[2];
    Xsi[0] = xsi[0];
    Xsi[1] = xsi[1];
       
    int degm = (*ipar_)[Npatch] -> getDegree(0);
    int degn = (*ipar_)[Npatch] -> getDegree(1);
    int npcm = (*ipar_)[Npatch] -> getNcp(0);
    int npcn = (*ipar_)[Npatch] -> getNcp(1); 


    int numLocalBF = (degm+1)*(degn+1);

    double uLeft [3];
    double uRight[3];
    double vLeft [3];
    double vRight[3];
    double phiL1 [3];
    double phiL2 [3];
    double dersL1[3][2];
    double dersL2[3][2];
    double phit  [9];
    double dphit [2][9];
    double ddphit [2][2][9];
    double aux1[2][3] = {}; // stores the two lines more recently computed
    double aux2[2][3] = {}; // stores the two lines more recently computed
    double uBF[3][3] = {}; // stores de base functions in each direction
    double vBF[3][3] = {}; // stores de base functions in each direction
    double saved,temp;

    double *uknot = (*ipar_)[Npatch] ->  getuKnot();
    double *vknot = (*ipar_)[Npatch] ->  getvKnot();  

    // indexes from index space where a function related with first local point of element start 
    int uind = inc [0];
    int vind = inc [1]; 

    // parametric coordinates that defines the element
    double u1 = uknot[uind];
    double u2 = uknot[uind+1];
    double v1 = vknot[vind];
    double v2 = vknot[vind +1];

    // relates integration space with the parametric one
    const double xsi1 =((u2 - u1) * Xsi[0] + (u2 + u1)) * 0.5;
    const double xsi2 =((v2 - v1) * Xsi[1] + (v2 + v1)) * 0.5;

    // Base function u direction
    uBF[0][0] = 1.0;
    for (int j = 1; j < (degm + 1); j++) {
        uLeft[j] = xsi1 - uknot[uind + 1 - j];
        uRight[j] = uknot[uind + j] - xsi1;
        saved = 0.;
        for (int r = 0; r < j; r++) {
            //lower triangle (Piegl and Tiller pg. 68-71)
            uBF[j][r] = uRight[r + 1] + uLeft[j - r];
            temp = uBF[r][j - 1]/uBF[j][r];
            //upper triangle
            uBF[r][j] = saved + uRight[r + 1]*temp;
            saved = uLeft[j - r] * temp;
        };
        uBF[j][j] = saved;
    };

    for (int i = 0; i <= degm ; i++ ) {
        phiL1[i] = uBF[i][degm];
    };

   // Base function v direction
    vBF[0][0] = 1.0;
    for (int j = 1; j< (degn + 1); j++) {
        vLeft[j] = xsi2 - vknot[vind + 1 - j];
        vRight[j] = vknot[vind + j] - xsi2;
        saved = 0.;
        for (int r = 0; r < j; r++) {
            //lower triangle (Piegl and Tiller pg. 68-71)
            vBF[j][r] = vRight[r + 1] + vLeft[j - r];
            temp = vBF[r][j - 1]/vBF[j][r];
            //upper triangle
            vBF[r][j] = saved + vRight[r + 1]*temp;
            saved = vLeft[j - r]*temp;
        };
        vBF[j][j] = saved;
    };

    for (int i = 0; i <= degn ; i++ ) {
        phiL2[i] = vBF[i][degn];
    };



    // derivatives u direction;
    int s1,s2,rk,pk,j1,j2,cor,k;
    double d;

    for (int r = 0; r <= degm; r++) { //number of functions
        s1 = 0;
        s2 = 1;
        aux1 [0][0] = 1.0;
        
        for (int k = 1; k <= degm; k++){ //ithDerivative = 2
            
            d = 0.0;
            rk = r - k;
            pk = degm - k;
            if (r >= k) {
                aux1[s2][0] = aux1[s1][0] / uBF[pk + 1][rk];
                d = aux1[s2][0] * uBF[rk][pk];
            };
            if (rk >= -1){
                j1 = 1;
            } else {
                j1 = -rk;
            };
            if ((r-1) <= pk) {
                j2 = k - 1 ;
            } else {
                j2 = degm - r;
            };
            for (int j = j1; j <= j2; j++) {
                aux1[s2][j] = (aux1[s1][j] - aux1[s1][j - 1]) / uBF[pk + 1][rk + j];
                d = d + aux1[s2][j] * uBF[rk + j][pk];
            };
            if (r <= pk){
                aux1[s2][k] = -aux1[s1][k - 1] / uBF[pk + 1][r];
                d = d + aux1[s2][k] * uBF[r][pk];
            }
            dersL1[r][k-1] = d;
            int j = s1;
            s1 = s2;
            s2 = j;
        };
    };

    int r1 = degm;
    for (int k = 1; k <= degm; k++) { //ithDerivative 
       
        for (int i = 0; i <= degm; i++) dersL1[i][k-1] *= r1;

        r1 *= (degm - k);    
    };


    // derivatives v direction;
    for (int r = 0; r <= degn; r++) {
        s1 = 0;
        s2 = 1;
        aux2 [0][0] = 1.0;
        for (int k = 1; k <= degn; k++){ //ith derivative
            d = 0.0;
            rk = r - k;
            pk = degn - k;
            if (r >= k) {
                aux2[s2][0] = aux2[s1][0] / vBF[pk + 1][rk];
                d = aux2[s2][0] * vBF[rk][pk];
            };
            if (rk >= -1){
                j1 = 1;
            } else{
                j1 = -rk;
            };
            if ((r-1) <= pk) {
                j2 = k - 1 ;
            } else {
                j2 = degn - r;
            };
            for (int j = j1; j <= j2; j++) {
                aux2[s2][j] = (aux2[s1][j] - aux2[s1][j - 1]) / vBF[pk + 1][rk + j];
                d = d + aux2[s2][j] * vBF[rk + j][pk];
            };
            if (r <= pk){
                aux2[s2][k] = -aux2[s1][k - 1] / vBF[pk + 1][r];
                d = d + aux2[s2][k] * vBF[r][pk];
            };
            dersL2 [r][k-1] = d;
            int j = s1;
            s1 = s2;
            s2 = j;
        };
    };

    r1 = degn;
    for (int k = 1; k <= degn; k++) {
        for (int i = 0; i <= degn; i++) dersL2[i][k-1] *= r1;
        r1 *= (degn - k);    
    };



    int index = 0;
    for (int j = 0; j<= degn; j ++) {
        for (int k = 0; k <= degm; k ++) {
           
            phit[index] = phiL1[k] * phiL2 [j];

            dphit[0][index] = dersL1[k][0] * phiL2 [j];
            dphit[1][index] = phiL1[k] * dersL2 [j][0];
            
            ddphit[0][0][index] = dersL1[k][1] * phiL2 [j]; 
            ddphit[0][1][index] = dersL1[k][0] * dersL2[j][0];
            ddphit[1][0][index] = dersL1[k][0] * dersL2[j][0];
            ddphit[1][1][index] = phiL1[k] * dersL2[j][1];

            index ++;
        };
    };


    double sumf = 0.0;
    double sumdqsi = 0.0;
    double sumdeta = 0.0;
    double sumdqsi_dqsi = 0.0;
    double sumdqsi_deta = 0.0;
    double sumdeta_dqsi = 0.0;
    double sumdeta_deta = 0.0; 



    for (int i = 0; i < numLocalBF; i++) {
        sumf += phit[i] * wpcs[i];
        sumdqsi += dphit[0][i] * wpcs[i];
        sumdeta += dphit[1][i] * wpcs[i];
        sumdqsi_dqsi += ddphit[0][0][i] * wpcs[i];
        sumdqsi_deta += ddphit[0][1][i] * wpcs[i];
        sumdeta_dqsi += ddphit[1][0][i] * wpcs[i];
        sumdeta_deta += ddphit[1][1][i] * wpcs[i]; 
    };



    for (int i = 0; i < numLocalBF; i++) {

        double u_ = dphit[0][i] * wpcs[i] * sumf - phit[i] * wpcs[i] * sumdqsi;
        double ul_ = ddphit[0][0][i] * wpcs[i] * sumf - phit[i] * wpcs[i] * sumdqsi_dqsi;
        double v_ = sumf * sumf;
        double vl_ = 2 * sumf * sumdqsi;
        ddphi[0][0][i] = ((ul_ * v_ - u_ * vl_)/ (v_ * v_));

        
        u_ = dphit[0][i] * wpcs[i] * sumf - phit[i] * wpcs[i] * sumdqsi;;
        ul_ = ddphit[0][1][i] * wpcs[i] * sumf + dphit[0][i] * wpcs[i] * sumdeta - 
              dphit[1][i] * wpcs[i] * sumdqsi - phit[i] * wpcs[i] * sumdqsi_deta;
        v_ = sumf * sumf;
        vl_ = 2 * sumf * sumdeta;
        ddphi[0][1][i] = ((ul_ * v_ - u_ * vl_)/ (v_ * v_));


        u_ = dphit[1][i] * wpcs[i] * sumf - phit[i] * wpcs[i] * sumdeta;;
        ul_ = ddphit[1][0][i] * wpcs[i] * sumf + dphit[1][i] * wpcs[i] * sumdqsi - 
              dphit[0][i] * wpcs[i] * sumdeta - phit[i] * wpcs[i] * sumdeta_dqsi;
        v_ = sumf * sumf;
        vl_ = 2 * sumf * sumdqsi;
        ddphi[1][0][i] = ((ul_ * v_ - u_ * vl_)/ (v_ * v_));


        u_ = dphit[1][i] * wpcs[i] * sumf - phit[i] * wpcs[i] * sumdeta;
        ul_ = ddphit[1][1][i] * wpcs[i] * sumf - phit[i] * wpcs[i] * sumdeta_deta;
        v_ = sumf * sumf;
        vl_ = 2 * sumf * sumdeta;
        ddphi[1][1][i] = ((ul_ * v_ - u_ * vl_)/ (v_ * v_));
    };

    return;

};

