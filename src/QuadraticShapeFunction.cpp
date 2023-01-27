#include "QuadraticShapeFunction.h"

//------------------------------------------------------------------------------
//-------------------------COMPUTE SHAPE FUNCTION VALUE-------------------------
//------------------------------------------------------------------------------
// Defines quadratic shape functions and its derivatives 
// for triangles and tetrahedrons


template<int DIM>
void QuadShapeFunction<DIM>::basisFunctionsBS(double &xsi, int &deg, double *knot, int &inc, double *phiL){

    double  uLeft [deg+1];
    double  uRight[deg+1];
    double  uBF[deg+1][deg+1] = {}; // stores de base functions in each direction

    double saved,temp;

    // parametric coordinates that defines the element
    double u1 = knot[inc];
    double u2 = knot[inc+1];

    // relates integration space with the parametric one
    double xsi1 =((u2 - u1) * xsi + (u2 + u1)) * 0.5;

    //Functions in u direction
    uBF[0][0] = 1.0;

    for (int j = 1; j < (deg+1); j++) {
        uLeft[j] = xsi1 - knot[inc+1-j];
        uRight[j] = knot[inc+j] - xsi1;
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

    for (int i = 0; i < deg +1 ; i++ ) {
        phiL[i] = uBF[i][deg];
    };

    return;
};

template<int DIM>
void QuadShapeFunction<DIM>::derBasisFunctionsBS(int &degree, double &xsi, int &deg, double *knot, int &inc, double *dphiL){

    double  uLeft [deg+1];
    double  uRight[deg+1];
    double  uBF[deg+1][deg+1] = {}; // stores de base functions in each direction
    double  aux[deg][deg+1] = {}; // stores the two lines more recently computed
    double  saved,temp,d;
    int s1,s2,rk,pk,j1,j2,cor,k;

    //parametric coordinates that defines the element
    double u1 = knot[inc];
    double u2 = knot[inc+1];

    //relates integration space with the parametric one
    double xsi1 =((u2 - u1) * xsi + (u2 + u1)) * 0.5;

    //Basis Functions
    uBF[0][0] = 1.0;
    for (int j = 1; j < (deg+1); j++) {
        uLeft[j] = xsi1 - knot[inc+1-j];
        uRight[j] = knot[inc+j] - xsi1;
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

    //Derivatives
    double dersL[deg+1][degree];
    for (int r = 0; r <= deg; r++) { //number of functions
        s1 = 0;
        s2 = 1;
        aux [0][0] = 1.0;
        
        for (int k = 1; k <= degree; k++){ //ithDerivative = 2
            
            d = 0.0;
            rk = r - k;
            pk = deg - k;
            if (r >= k) {
                aux[s2][0] = aux[s1][0] / uBF[pk + 1][rk];
                d = aux[s2][0] * uBF[rk][pk];
            };
            if (rk >= -1){
                j1 = 1;
            } else {
                j1 = -rk;
            };
            if ((r-1) <= pk) {
                j2 = k - 1 ;
            } else {
                j2 = deg - r;
            };
            for (int j = j1; j <= j2; j++) {
                aux[s2][j] = (aux[s1][j] - aux[s1][j - 1]) / uBF[pk + 1][rk + j];
                d = d + aux[s2][j] * uBF[rk + j][pk];
            };
            if (r <= pk){
                aux[s2][k] = -aux[s1][k - 1] / uBF[pk + 1][r];
                d = d + aux[s2][k] * uBF[r][pk];
            }
            dersL[r][k-1] = d;
            int j = s1;
            s1 = s2;
            s2 = j;
        };
    };

    int r1 = deg;
    for (int k = 1; k <= degree; k++) { //ithDerivative 
       
        for (int i = 0; i <= deg; i++) dersL[i][k-1] *= r1;

        r1 *= (deg - k);    
    };

    for (int i = 0; i < deg+1; i++){
        dphiL[i] = dersL[i][degree-1];
    }

    
    return;
};



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
};

template<>
void QuadShapeFunction<3>::evaluateFem(double *xsi, double *phi) const {
    
     const double xsi1 = xsi[0];
     const double xsi2 = xsi[1];
     const double xsi3 = xsi[2];

     phi[0] = 2.0 * (xsi1 - 0.50) * xsi1; 
     phi[1] = 2.0 * (xsi2 - 0.50) * xsi2; 
     phi[2] = 2.0 * (xsi3 - 0.50) * xsi3; 
     phi[3] = (2.0 - 2.0 * xsi3 - 2.0 * xsi2 - 2.0 * xsi1 - 1.0) * (1.0 - xsi1 - xsi2 - xsi3); 
     phi[4] = 4.0 * xsi1 * xsi2; 
     phi[5] = 4.0 * xsi2 * xsi3;
     phi[6] = 4.0 * xsi1 * xsi3;
     phi[7] = 4.0 * xsi1 * (1.0 - xsi1 - xsi2 - xsi3);
     phi[8] = 4.0 * xsi3 * (1.0 - xsi1 - xsi2 - xsi3);
     phi[9] = 4.0 * xsi2 * (1.0 - xsi1 - xsi2 - xsi3);

     return;
};


template<int DIM>
void QuadShapeFunction<DIM>::evaluateIso(double *xsi, double *phi, double *wpcs, int *inc, std::vector<IParameters_ *> &iparameters ,int Npatch) {
    
    std::vector<IParameters_ *> *ipar_ = &iparameters;

    int deg[DIM],npc[DIM];

    for (int i = 0; i < DIM; i++){
        deg[i] = (*ipar_)[Npatch] -> getDegree(i);
        npc[i] = (*ipar_)[Npatch] -> getNcp(i);
    };

    int numLocalBF = 1.;

    for (int i = 0; i < DIM; i++){
        numLocalBF *= (deg[i]+1);
    }

    //assuming same functions degree in all directions
    double phiL [DIM][deg[0]+1];
    double phit  [numLocalBF];

    for (int i = 0; i < DIM; i++){
        int size = deg[i]+npc[i]+1;
        double knot[size];
        (*ipar_)[Npatch] -> getKnot(i,size,knot);
        double phil[deg[i]+1];
        basisFunctionsBS(xsi[i],deg[i],knot,inc[i],phil);
        for (int j = 0; j < deg[i]+1; j++){
            phiL[i][j] = phil[j];
        };
    };
    
    int index = 0;
    if(DIM == 2){
        for (int j = 0; j< deg[1]+1; j ++) {
            for (int k = 0; k < deg[0]+1; k ++) {
                phit[index] = phiL[0][k] * phiL[1][j];           
                index ++;
            };
        };  
    } else { //DIM = 3
        for (int i = 0; i < deg[2]+1 ; i++){
            for (int j = 0; j< deg[1]+1; j ++) {
                for (int k = 0; k < deg[0]+1; k ++) {
                    phit[index] = phiL[0][k] * phiL[1][j] * phiL[2][i];                
                    index ++;
                };
            }; 
        };
    };

    double sumF = 0.0;
    for (int i = 0; i < numLocalBF; i++) {
        sumF += phit[i] * wpcs[i];
    };

    for (int i = 0; i < numLocalBF; i++) {
        for (int j = 0; j < DIM; j++) phi[i] = phit[i] * wpcs[i]/ (sumF * sumF); 
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

};


template<>
void QuadShapeFunction<2>::evaluateGradientIso(double *xsi, double **dphi, double *wpcs, int *inc, std::vector<IParameters_ *> &iparameters ,int Npatch) {

    const int DIM = 2;
    std::vector<IParameters_ *> *ipar_ = &iparameters;

    int deg[DIM],npc[DIM];

    for (int i = 0; i < DIM; i++){
        deg[i] = (*ipar_)[Npatch] -> getDegree(i);
        npc[i] = (*ipar_)[Npatch] -> getNcp(i);
    };

    int numLocalBF = 1.;

    for (int i = 0; i < DIM; i++){
        numLocalBF *= (deg[i]+1);
    }

    //assuming same functions degree in all directions
    double phiL [DIM][deg[0]+1];
    double dphiL[DIM][deg[0]+1];
    double phit  [numLocalBF];
    double dphit [DIM][numLocalBF];


    for (int i = 0; i < DIM; i++){
        int size = deg[i]+npc[i]+1;
        double knot[size];
        (*ipar_)[Npatch] -> getKnot(i,size,knot);
        double phil[deg[i]+1], dphil[deg[i]+1];
        basisFunctionsBS(xsi[i],deg[i],knot,inc[i],phil);
        int degree = 1;
        derBasisFunctionsBS(degree,xsi[i],deg[i],knot,inc[i],dphil);
        for (int j = 0; j < deg[i]+1; j++){
            phiL[i][j] = phil[j];
            dphiL[i][j] = dphil[j];
        };
    };

    int index = 0;
    for (int j = 0; j< deg[1]+1; j ++) {
        for (int k = 0; k < deg[0]+1; k ++) {
        
            phit[index] = phiL[0][k] * phiL[1][j];

            dphit[0][index] = dphiL[0][k] * phiL[1][j];
            dphit[1][index] = phiL[0][k] * dphiL[1][j];
            
            index ++;
        };
    }; 
     
    double sumF = 0.0;
    double sumDer[DIM] = {};

    for (int i = 0; i < numLocalBF; i++) {
        sumF += phit[i] * wpcs[i];
        for (int j = 0; j < DIM; j++){
            sumDer[j] += dphit[j][i] * wpcs[i];
        };
    };

    for (int i = 0; i < numLocalBF; i++) {
        for (int j = 0; j < DIM; j++) dphi[j][i] = ((dphit[j][i] * wpcs[i] * sumF - phit[i] * wpcs[i] * sumDer[j]) / (sumF * sumF)); 
    };

    return;
};


template<>
void QuadShapeFunction<3>::evaluateGradientIso(double *xsi, double **dphi, double *wpcs, int *inc, std::vector<IParameters_ *> &iparameters ,int Npatch) {


    const int DIM = 3;

    std::vector<IParameters_ *> *ipar_ = &iparameters;

    int deg[DIM],npc[DIM];

    for (int i = 0; i < DIM; i++){
        deg[i] = (*ipar_)[Npatch] -> getDegree(i);
        npc[i] = (*ipar_)[Npatch] -> getNcp(i);
    };

    int numLocalBF = 1.;

    for (int i = 0; i < DIM; i++){
        numLocalBF *= (deg[i]+1);
    }

    //assuming same functions degree in all directions
    double phiL [DIM][deg[0]+1];
    double dphiL[DIM][deg[0]+1];
    double phit  [numLocalBF];
    double dphit [DIM][numLocalBF];


    for (int i = 0; i < DIM; i++){
        int size = deg[i]+npc[i]+1;
        double knot[size];
        (*ipar_)[Npatch] -> getKnot(i,size,knot);
        double phil[deg[i]+1], dphil[deg[i]+1];
        basisFunctionsBS(xsi[i],deg[i],knot,inc[i],phil);
        int degree = 1;
        derBasisFunctionsBS(degree,xsi[i],deg[i],knot,inc[i],dphil);
        for (int j = 0; j < deg[i]+1; j++){
            phiL[i][j] = phil[j];
            dphiL[i][j] = dphil[j];
        };
    };

    int index = 0;

    for (int i = 0; i < deg[2]+1 ; i++){
        for (int j = 0; j< deg[1]+1; j ++) {
            for (int k = 0; k < deg[0]+1; k ++) {
            
                phit[index] = phiL[0][k] * phiL[1][j] * phiL[2][i];

                dphit[0][index] = dphiL[0][k] * phiL[1][j] * phiL[2][i];
                dphit[1][index] = phiL[0][k] * dphiL[1][j] * phiL[2][i];
                dphit[2][index] = phiL[0][k] * phiL[1][j] * dphiL[2][i];
                
                index ++;
            };
        }; 
    };
     
    double sumF = 0.0;
    double sumDer[DIM] = {};

    for (int i = 0; i < numLocalBF; i++) {
        sumF += phit[i] * wpcs[i];
        for (int j = 0; j < DIM; j++){
            sumDer[j] += dphit[j][i] * wpcs[i];
        };
    };

    for (int i = 0; i < numLocalBF; i++) {
        for (int j = 0; j < DIM; j++) dphi[j][i] = ((dphit[j][i] * wpcs[i] * sumF - phit[i] * wpcs[i] * sumDer[j]) / (sumF * sumF)); 
    };

    

    return;
};


template<>
void QuadShapeFunction<2>::evaluateHessianFem(double ***ddphi){

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
void QuadShapeFunction<3>::evaluateHessianFem(double ***ddphi){

    ddphi[0][0][0] = 4.;
    ddphi[0][1][0] = 0.;
    ddphi[0][2][0] = 0.;

    ddphi[1][0][0] = 0.;
    ddphi[1][1][0] = 0.;
    ddphi[1][2][0] = 0.;

    ddphi[2][0][0] = 0.;
    ddphi[2][1][0] = 0.;
    ddphi[2][2][0] = 0.;


    ddphi[0][0][1] = 0.;
    ddphi[0][1][1] = 0.;
    ddphi[0][2][1] = 0.;

    ddphi[1][0][1] = 0.;
    ddphi[1][1][1] = 4.;
    ddphi[1][2][1] = 0.;

    ddphi[2][0][1] = 0.;
    ddphi[2][1][1] = 0.;
    ddphi[2][2][1] = 0.;

    
    ddphi[0][0][2] = 0.;
    ddphi[0][1][2] = 0.;
    ddphi[0][2][2] = 0.;

    ddphi[1][0][2] = 0.;
    ddphi[1][1][2] = 0.;
    ddphi[1][2][2] = 0.;

    ddphi[2][0][2] = 0.;
    ddphi[2][1][2] = 0.;
    ddphi[2][2][2] = 4.;


    ddphi[0][0][3] = 4.;
    ddphi[0][1][3] = 4.;
    ddphi[0][2][3] = 4.;

    ddphi[1][0][3] = 4.;
    ddphi[1][1][3] = 4.;
    ddphi[1][2][3] = 4.;

    ddphi[2][0][3] = 4.;
    ddphi[2][1][3] = 4.;
    ddphi[2][2][3] = 4.;


    ddphi[0][0][4] = 0.;
    ddphi[0][1][4] = 4.;
    ddphi[0][2][4] = 0.;

    ddphi[1][0][4] = 4.;
    ddphi[1][1][4] = 0.;
    ddphi[1][2][4] = 0.;

    ddphi[2][0][4] = 0.;
    ddphi[2][1][4] = 0.;
    ddphi[2][2][4] = 0.;

    
    ddphi[0][0][5] = 0.;
    ddphi[0][1][5] = 0.;
    ddphi[0][2][5] = 0.;

    ddphi[1][0][5] = 0.;
    ddphi[1][1][5] = 0.;
    ddphi[1][2][5] = 4.;

    ddphi[2][0][5] = 0.;
    ddphi[2][1][5] = 4.;
    ddphi[2][2][5] = 0.;


    ddphi[0][0][6] = 0.;
    ddphi[0][1][6] = 0.;
    ddphi[0][2][6] = 4.;

    ddphi[1][0][6] = 0.;
    ddphi[1][1][6] = 0.;
    ddphi[1][2][6] = 0.;

    ddphi[2][0][6] = 4.;
    ddphi[2][1][6] = 0.;
    ddphi[2][2][6] = 0.;


    ddphi[0][0][7] = -8.;
    ddphi[0][1][7] = -4.;
    ddphi[0][2][7] = -4.;

    ddphi[1][0][7] = -4.;
    ddphi[1][1][7] = 0.;
    ddphi[1][2][7] = 0.;

    ddphi[2][0][7] = -4.;
    ddphi[2][2][7] = 0.;


    ddphi[0][0][8] =  0.;
    ddphi[0][1][8] =  0.;
    ddphi[0][2][8] = -4.;

    ddphi[1][0][8] =  0.;
    ddphi[1][1][8] =  0.;
    ddphi[1][2][8] = -4.;

    ddphi[2][0][8] = -4.;
    ddphi[2][1][8] = -4.;
    ddphi[2][2][8] = -8.;


    ddphi[0][0][9] = 0.;
    ddphi[0][1][9] = -4.;
    ddphi[0][2][9] = 0.;
    
    ddphi[1][0][9] =-4.;
    ddphi[1][0][9] = -8.;
    ddphi[1][0][9] = -4.;
    
    ddphi[2][0][9] = -0.;
    ddphi[2][0][9] = -4.;
    ddphi[2][0][9] = 0.;

};


template<int DIM>
void QuadShapeFunction<DIM>::evaluateHessianIso(double *xsi, double*** ddphi, double* wpcs, int *inc, std::vector<IParameters_ *> &iparameters ,int Npatch){
   
    std::vector<IParameters_ *> *ipar_ = &iparameters;

    int deg[DIM],npc[DIM];

    for (int i = 0; i < DIM; i++){
        deg[i] = (*ipar_)[Npatch] -> getDegree(i);
        npc[i] = (*ipar_)[Npatch] -> getNcp(i);
    };

    int numLocalBF = 1.;

    for (int i = 0; i < DIM; i++){
        numLocalBF *= (deg[i]+1);
    }

    //assuming same functions degree in all directions
    double phiL [DIM][deg[0]+1];
    double dphiL[DIM][deg[0]+1];
    double ddphiL [DIM][deg[0]+1];
    double phit  [numLocalBF];
    double dphit [DIM][numLocalBF];
    double ddphit [DIM][DIM][numLocalBF];


    for (int i = 0; i < DIM; i++){
        int size = deg[i]+npc[i]+1;
        double knot[size];
        (*ipar_)[Npatch] -> getKnot(i,size,knot);
        double phil[deg[i]+1], dphil[deg[i]+1], ddphil[deg[i]+1];
        basisFunctionsBS(xsi[i],deg[i],knot,inc[i],phil);
        int degree = 1;
        derBasisFunctionsBS(degree,xsi[i],deg[i],knot,inc[i],dphil);
        degree = 2;
        derBasisFunctionsBS(degree,xsi[i],deg[i],knot,inc[i],ddphil);
        for (int j = 0; j < deg[i]+1; j++){
            phiL[i][j] = phil[j];
            dphiL[i][j] = dphil[j];
            ddphiL[i][j] = ddphil[j];
        };
    };

    int index = 0;
    if (DIM == 2){
        for (int j = 0; j< deg[1]+1; j ++) {
            for (int k = 0; k < deg[0]+1; k ++) {
            
                phit[index] = phiL[0][k] * phiL[1][j];

                dphit[0][index] = dphiL[0][k] * phiL[1][j];
                dphit[1][index] = phiL[0][k] * dphiL[1][j];
                
                ddphit[0][0][index] = ddphiL[0][k] * phiL[1][j]; 
                ddphit[0][1][index] = dphiL[0][k] * dphiL[1][j];
                
                ddphit[1][0][index] = dphiL[0][k] * dphiL[1][j];
                ddphit[1][1][index] = phiL[0][k] * ddphiL[1][j];

                index ++;
            };
        };
    } else { //DIM = 3
        for (int i = 0; i < deg[2]+1; i++){
            for (int j = 0; j< deg[1]+1; j ++) {
                for (int k = 0; k < deg[0]+1; k ++) {
                
                    phit[index] = phiL[0][k] * phiL[1][j] * phiL[2][i];

                    dphit[0][index] = dphiL[0][k] * phiL[1][j] * phiL[2][i];
                    dphit[1][index] = phiL[0][k] * dphiL[1][j] * phiL[2][i];
                    dphit[2][index] = phiL[0][k] * phiL[1][j] * dphiL[2][i];
                    
                    ddphit[0][0][index] = ddphiL[0][k] * phiL[1][j] * phiL[2][i]; 
                    ddphit[0][1][index] = dphiL[0][k] * dphiL[1][j] * phiL[2][i];
                    ddphit[0][2][index] = dphiL[0][k] * phiL[1][j] * dphiL[2][i];
                
                    ddphit[1][0][index] = dphiL[0][k] * dphiL[1][j] * phiL[2][i];
                    ddphit[1][1][index] = phiL[0][k] * ddphiL[1][j] * phiL[2][i];
                    ddphit[1][2][index] = phiL[0][k] * dphiL[1][j] * dphiL[2][i];

                    ddphit[2][0][index] = dphiL[0][k] * phiL[1][j] * dphiL[2][i];
                    ddphit[2][1][index] = phiL[0][k] * dphiL[1][j] * dphiL[2][i];
                    ddphit[2][2][index] = phiL[0][k] * phiL[1][j] * ddphiL[2][i];

                    index ++;
                };
            };
        };
    };
  
    double sumF = 0.0;
    double sumDer[DIM] = {};
    double sumSecDer[DIM][DIM] = {};

    for (int i = 0; i < numLocalBF; i++) {
        sumF += phit[i] * wpcs[i];
        for (int j = 0; j < DIM; j++){
            sumDer[j] += dphit[j][i] * wpcs[i];
            for (int k = 0; k < DIM; k++){
                sumSecDer[j][k] += ddphit[j][k][i] * wpcs[i];
            };
        };
    };

    double v_ = sumF * sumF;
    for (int i = 0; i < numLocalBF; i++) {
        for (int j = 0; j < DIM; j++){
            double u_ = dphit[j][i] * wpcs[i] * sumF - phit[i] * wpcs[i] * sumDer[j];
            for (int k = 0; k < DIM; k++){
                double ul_ = ddphit[j][k][i] * wpcs[i] * sumF + dphit[j][i] * wpcs[i] * sumDer[k] - 
                             dphit[k][i] * wpcs[i] * sumDer[j] - phit[i] * wpcs[i] * sumSecDer[j][k];;
                double vl_ = 2 * sumF * sumDer[k];
                ddphi[j][k][i] = ((ul_ * v_ - u_ * vl_)/ (v_ * v_));
            };
        };
    };

    return;

};


template<>
double QuadShapeFunction<2>::ParametricCoordBezier(int i, int j){

    paramCoord[0][0] = -1.; paramCoord[0][1] = -1.;
    paramCoord[1][0] = 0.; paramCoord[1][1] = -1.;
    paramCoord[2][0] = 1.; paramCoord[2][1] = -1.;
    paramCoord[3][0] = -1.; paramCoord[3][1] = 0.;
    paramCoord[4][0] = 0.; paramCoord[4][1] = 0.;
    paramCoord[5][0] = 1.; paramCoord[5][1] = 0.;
    paramCoord[6][0] = -1.; paramCoord[6][1] = 1.;
    paramCoord[7][0] = 0.; paramCoord[7][1] = 1.;
    paramCoord[8][0] = 1.; paramCoord[8][1] = 1.;  

    return paramCoord[i][j];
};

template<>
double QuadShapeFunction<3>::ParametricCoordBezier(int i, int j){

    paramCoord[0][0] = -1.; paramCoord[0][1] = -1.; paramCoord[0][2] = -1.;
    paramCoord[1][0] = 0.; paramCoord[1][1] = -1.; paramCoord[1][2] = -1.;
    paramCoord[2][0] = 1.; paramCoord[2][1] = -1.; paramCoord[2][2] = -1.;
    paramCoord[3][0] = -1.; paramCoord[3][1] = 0.; paramCoord[3][2] = -1.;
    paramCoord[4][0] = 0.; paramCoord[4][1] = 0.; paramCoord[4][2] = -1.;
    paramCoord[5][0] = 1.; paramCoord[5][1] = 0.; paramCoord[5][2] = -1.;
    paramCoord[6][0] = -1.; paramCoord[6][1] = 1.; paramCoord[6][2] = -1.;
    paramCoord[7][0] = 0.; paramCoord[7][1] = 1.; paramCoord[7][2] = -1.;
    paramCoord[8][0] = 1.; paramCoord[8][1] = 1.; paramCoord[8][2] = -1.;  

    paramCoord[9][0] = -1.; paramCoord[9][1] = -1.; paramCoord[9][2] = 0.;
    paramCoord[10][0] = 0.; paramCoord[10][1] = -1.; paramCoord[10][2] = 0.;
    paramCoord[11][0] = 1.; paramCoord[11][1] = -1.; paramCoord[11][2] = 0.;
    paramCoord[12][0] = -1.; paramCoord[12][1] = 0.; paramCoord[12][2] = 0.;
    paramCoord[13][0] = 0.; paramCoord[13][1] = 0.; paramCoord[13][2] = 0.;
    paramCoord[14][0] = 1.; paramCoord[14][1] = 0.; paramCoord[14][2] = 0.;
    paramCoord[15][0] = -1.; paramCoord[15][1] = 1.; paramCoord[15][2] = 0.;
    paramCoord[16][0] = 0.; paramCoord[16][1] = 1.; paramCoord[16][2] = 0.;
    paramCoord[17][0] = 1.; paramCoord[17][1] = 1.; paramCoord[17][2] = 0.;  

    paramCoord[18][0] = -1.; paramCoord[18][1] = -1.; paramCoord[18][2] = 1.;
    paramCoord[19][0] = 0.; paramCoord[19][1] = -1.; paramCoord[19][2] = 1.;
    paramCoord[20][0] = 1.; paramCoord[20][1] = -1.; paramCoord[20][2] = 1.;
    paramCoord[21][0] = -1.; paramCoord[21][1] = 0.; paramCoord[21][2] = 1.;
    paramCoord[22][0] = 0.; paramCoord[22][1] = 0.; paramCoord[22][2] = 1.;
    paramCoord[23][0] = 1.; paramCoord[23][1] = 0.; paramCoord[23][2] = 1.;
    paramCoord[24][0] = -1.; paramCoord[24][1] = 1.; paramCoord[24][2] = 1.;
    paramCoord[25][0] = 0.; paramCoord[25][1] = 1.; paramCoord[25][2] = 1.;
    paramCoord[26][0] = 1.; paramCoord[26][1] = 1.; paramCoord[26][2] = 1.;  

    return paramCoord[i][j];
}

template class QuadShapeFunction<2>;
template class QuadShapeFunction<3>;

