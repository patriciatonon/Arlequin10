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
}


template<>
void QuadShapeFunction<2>::evaluateIso(double *xsi, double *phi, double *wpcs, int *inc, std::vector<IParameters_ *> &iparameters ,int Npatch) {
    

    std::vector<IParameters_ *> *ipar_ = &iparameters;
       
    int degm = (*ipar_)[Npatch] -> getDegree(0);
    int degn = (*ipar_)[Npatch] -> getDegree(1);
    int npcm = (*ipar_)[Npatch] -> getNcp(0);
    int npcn = (*ipar_)[Npatch] -> getNcp(1);   
    int numLocalBF = (degm+1)*(degn+1);

    double *uknot = (*ipar_)[Npatch]->getuKnot();
    double *vknot = (*ipar_)[Npatch]->getvKnot(); 

    double  phiL1 [degm+1];
    double  phiL2 [degn+1];
    double  phit  [numLocalBF];


    basisFunctionsBS(xsi[0],degm,uknot,inc[0],phiL1);
    basisFunctionsBS(xsi[1],degn,vknot,inc[1],phiL2);
    

    int index = 0;
    for (int j = 0; j < degn + 1 ; j ++) {
        for (int k = 0; k < degm + 1; k ++) {
            phit[index] = phiL1[k] * phiL2 [j];
            index ++;
        }
    }

    double sum = 0.0;
    for (int i = 0; i < numLocalBF; i++) {
        sum += phit[i] * wpcs[i];
    }

    for (int i = 0; i < numLocalBF; i++) {
        phi[i] = phit[i] * wpcs[i]/sum;
    }

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
    
    std::vector<IParameters_ *> *ipar_ = &iparameters;
       
    int degm = (*ipar_)[Npatch] -> getDegree(0);
    int degn = (*ipar_)[Npatch] -> getDegree(1);
    int degq = (*ipar_)[Npatch] -> getDegree(2);
    int npcm = (*ipar_)[Npatch] -> getNcp(0);
    int npcn = (*ipar_)[Npatch] -> getNcp(1);   
    int npcq = (*ipar_)[Npatch] -> getNcp(2); 

    int numLocalBF = (degm+1)*(degn+1)*(degq+1);

    double phiL1 [degm+1];
    double phiL2 [degn+1];
    double phiL3 [degq+1];
    double phit  [numLocalBF];

    double *uknot = (*ipar_)[Npatch]->getuKnot();
    double *vknot = (*ipar_)[Npatch]->getvKnot(); 
    double *tknot = (*ipar_)[Npatch]->gettKnot(); 

    basisFunctionsBS(xsi[0],degm,uknot,inc[0],phiL1);
    basisFunctionsBS(xsi[1],degn,vknot,inc[1],phiL2);
    basisFunctionsBS(xsi[2],degq,tknot,inc[2],phiL3);


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
void QuadShapeFunction<2>::evaluateGradientIso(double *xsi, double **dphi, double *wpcs, int *inc, std::vector<IParameters_ *> &iparameters ,int Npatch) {

    std::vector<IParameters_ *> *ipar_ = &iparameters;
      
    int degm = (*ipar_)[Npatch] -> getDegree(0);
    int degn = (*ipar_)[Npatch] -> getDegree(1);
    int npcm = (*ipar_)[Npatch] -> getNcp(0);
    int npcn = (*ipar_)[Npatch] -> getNcp(1);  
 
    int numLocalBF = (degm+1)*(degn+1);

    double phiL1 [degm+1];
    double phiL2 [degn+1];
    double phit  [numLocalBF];
    double dphiL1[degm+1];
    double dphiL2[degn+1];
    double dphit [2][numLocalBF];

    double *uknot = (*ipar_)[Npatch] -> getuKnot();
    double *vknot = (*ipar_)[Npatch] -> getvKnot();  

    basisFunctionsBS(xsi[0],degm,uknot,inc[0],phiL1);
    basisFunctionsBS(xsi[1],degn,vknot,inc[1],phiL2);

    int index = 0;
    for (int j = 0; j< (degn + 1) ; j ++) {
        for (int k = 0; k < (degm + 1); k ++) {
            phit[index] = phiL1[k] * phiL2 [j];
            index ++;
        };
    };

    int degree = 1; //firs derivative

    derBasisFunctionsBS(degree,xsi[0],degm,uknot,inc[0],dphiL1);
    derBasisFunctionsBS(degree,xsi[1],degn,vknot,inc[1],dphiL2);

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


    std::vector<IParameters_ *> *ipar_;  // isogeometric parameters
    ipar_ = &iparameters;
       
    int degm = (*ipar_)[Npatch] -> getDegree(0);
    int degn = (*ipar_)[Npatch] -> getDegree(1);
    int degq = (*ipar_)[Npatch] -> getDegree(2);
    int npcm = (*ipar_)[Npatch] -> getNcp(0);
    int npcn = (*ipar_)[Npatch] -> getNcp(1);  
    int npcq = (*ipar_)[Npatch] -> getNcp(2);
 
    int numLocalBF = (degm+1)*(degn+1)*(degq+1);

    double phiL1 [degm+1];
    double phiL2 [degn+1];
    double phiL3 [degq+1];
    double phit  [numLocalBF];
    double dphiL1[degm+1];
    double dphiL2[degn+1];
    double dphiL3[degq+1];
    double dphit [3][numLocalBF];

    double *uknot = (*ipar_)[Npatch] ->  getuKnot();
    double *vknot = (*ipar_)[Npatch] ->  getvKnot();  
    double *tknot = (*ipar_)[Npatch] ->  gettKnot();

    basisFunctionsBS(xsi[0],degm,uknot,inc[0],phiL1);
    basisFunctionsBS(xsi[1],degn,vknot,inc[1],phiL2);
    basisFunctionsBS(xsi[2],degq,tknot,inc[2],phiL3);

    int index = 0;
    for (int i = 0; i< (degq + 1) ; i ++) {
        for (int j = 0; j< (degn + 1) ; j ++) {
            for (int k = 0; k < (degm + 1); k ++) {
                phit[index] = phiL1[k] * phiL2 [j] * phiL3 [i];
                index ++;
            };
        };
    };
        
    int degree = 1; //first derivative
    derBasisFunctionsBS(degree,xsi[0],degm,uknot,inc[0],dphiL1);
    derBasisFunctionsBS(degree,xsi[1],degn,vknot,inc[1],dphiL2);
    derBasisFunctionsBS(degree,xsi[2],degq,tknot,inc[2],dphiL3);

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
void QuadShapeFunction<2>::evaluateHessianIso(double *xsi, double*** ddphi, double* wpcs, int *inc, std::vector<IParameters_ *> &iparameters ,int Npatch){

    std::vector<IParameters_ *> *ipar_ = &iparameters;
       
    int degm = (*ipar_)[Npatch] -> getDegree(0);
    int degn = (*ipar_)[Npatch] -> getDegree(1);
    int npcm = (*ipar_)[Npatch] -> getNcp(0);
    int npcn = (*ipar_)[Npatch] -> getNcp(1); 

    int numLocalBF = (degm+1)*(degn+1);

    double phiL1 [degm+1];
    double phiL2 [degn+1];
    double dphiL1 [degm+1];
    double dphiL2 [degn+1];
    double ddphiL1 [degm+1];
    double ddphiL2 [degn+1];
    double phit  [numLocalBF];
    double dphit [2][numLocalBF];
    double ddphit [2][2][numLocalBF];

    double *uknot = (*ipar_)[Npatch] ->  getuKnot();
    double *vknot = (*ipar_)[Npatch] ->  getvKnot();  

    basisFunctionsBS(xsi[0],degm,uknot,inc[0],phiL1);
    basisFunctionsBS(xsi[1],degn,vknot,inc[1],phiL2);

    int degree = 1; //first derivative
    derBasisFunctionsBS(degree,xsi[0],degm,uknot,inc[0],dphiL1);
    derBasisFunctionsBS(degree,xsi[1],degn,vknot,inc[1],dphiL2);

    degree = 2; //second derivatives
    derBasisFunctionsBS(degree,xsi[0],degm,uknot,inc[0],ddphiL1);
    derBasisFunctionsBS(degree,xsi[1],degn,vknot,inc[1],ddphiL2);

    int index = 0;
    for (int j = 0; j<= degn; j ++) {
        for (int k = 0; k <= degm; k ++) {
           
            phit[index] = phiL1[k] * phiL2 [j];

            dphit[0][index] = dphiL1[k] * phiL2 [j];
            dphit[1][index] = phiL1[k] * dphiL2[j];
            
            ddphit[0][0][index] = ddphiL1[k] * phiL2 [j]; 
            ddphit[0][1][index] = dphiL1[k] * dphiL2[j];
            ddphit[1][0][index] = dphiL1[k] * dphiL2[j];
            ddphit[1][1][index] = phiL1[k] * ddphiL2[j];

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


