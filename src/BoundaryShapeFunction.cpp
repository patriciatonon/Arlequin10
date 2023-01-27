
#include "BoundaryShapeFunction.h"

//------------------------------------------------------------------------------
//-------------------------COMPUTE SHAPE FUNCTION VALUE-------------------------
//------------------------------------------------------------------------------

// Evaluates the shape function value in a line boundary
template<>
void BoundShapeFunction<2>::evaluateBoundaryFem(double *Xsi, double *phiB_) {

    double xsi;
    xsi = Xsi[0];

    double aux;
    double Nnos = 3;
    double AdimCoord[3];
    AdimCoord[0] = -1.;
    AdimCoord[1] =  0.;
    AdimCoord[2] =  1.;

    for (int j = 0; j < Nnos; j++) {
        phiB_[j] = 1.0;
        for (int i = 0; i < Nnos; i++) {
            aux = 1.0;
            if(i != j){                
                phiB_[j]= phiB_[j] * (xsi - AdimCoord[i]) / 
                    (AdimCoord[j] - AdimCoord[i]);                
            };
        };
    };


    return;
}; 



/// Returns the shape function value in a surface boundary
template<>
void BoundShapeFunction<3>::evaluateBoundaryFem(double *Xsi, double *phiB_) {
    
    const double xsi1 = Xsi[0];
    const double xsi2 = Xsi[1];
    const double xsi3 = 1. - xsi1 - xsi2;

    phiB_[0] = xsi3 * (2.0 * xsi3 - 1.0);
    phiB_[1] = xsi1 * (2.0 * xsi1 - 1.0);
    phiB_[2] = xsi2 * (2.0 * xsi2 - 1.0);
    phiB_[3] = 4.0 * xsi3 * xsi1;
    phiB_[4] = 4.0 * xsi1 * xsi2;
    phiB_[5] = 4.0 * xsi2 * xsi3;

};

template<>
void BoundShapeFunction<2>::evaluateBoundaryIso(double *Xsi, double *phiB_, double *wpcs, int side, int *inc, std::vector<IParameters_ *> &iparameters, int Npatch){

    int           s = side;    
    //exemplo of type and side
            // side 2
            //7  8  9
    //side 3//4  5  6 //side 1
            //1  2  3
            //side 0

    double xsi;
    xsi = Xsi[0];

    std::vector<IParameters_ *> *ipar_;
    ipar_ = &iparameters;

    if ((s == 0) || (s == 2)) {

        // isogemetric data
        int degm = (*ipar_)[Npatch] -> getDegree(0);
        int npcm = (*ipar_)[Npatch]-> getNcp(0);
        double *uknot;
        double uLeft [3];
        double uRight[3];
        double phiL1 [3];
        double uBF [3][3] = {}; // stores de base functions in each direction
        double saved,temp;

        uknot = (*ipar_)[Npatch] -> getuKnot();

        // indexes from index space where a function related with first local point of element start 
        int uind = inc[0];

        // parametric coordinates that defines the element
        double u1 = uknot[uind];
        double u2 = uknot[uind+1];

        // relates integration space with the parametric one
        const double xsi1 =((u2 - u1) * xsi + (u2 + u1)) * 0.5;
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

        for (int i = 0; i < degm + 1 ; i++ ) {
            phiL1[i] = uBF[i][degm];
        };

        double sum = 0.;
        for (int i = 0; i < degm + 1 ; i++ ) {
            sum += phiL1[i] * wpcs[i];
        };

        for (int i = 0; i < degm + 1 ; i++ ) {
             phiB_[i] = phiL1[i] * wpcs[i] / sum;
        };
    };

    if ((s == 1) || (s == 3)) {

        // isogemetric data
        int degn = (*ipar_)[Npatch] -> getDegree(1);
        int npcn = (*ipar_)[Npatch] -> getNcp(1);   
        int dimv = degn + npcn +1;
        double *vknot;
        double vLeft [3];
        double vRight[3];
        double phiL2 [3];
        double vBF[3][3]= {}; // stores de base functions in each direction
        double saved,temp;

        vknot = (*ipar_)[Npatch] -> getvKnot();  

        // indexes from index space where a function related with first local point of element start 
        int vind = inc[1]; 

        // parametric coordinates that defines the element
        double v1 = vknot[vind];
        double v2 = vknot[vind +1];

        // relates integration space with the parametric one
        const double xsi2 =((v2 - v1) * xsi + (v2 + v1)) * 0.5; 

        vBF[0][0] = 1.0;
        for (int j = 1; j < (degn + 1); j++) {
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

        for (int i = 0; i < degn + 1 ; i++ ) {
            phiL2[i] = vBF[i][degn];
        };

        double sum = 0.0;

        for (int i = 0; i < degn + 1 ; i++ ) {
            sum += phiL2[i] * wpcs[i];
        };

        for (int i = 0; i < degn + 1 ; i++ ) {
             phiB_[i] = phiL2[i] * wpcs[i]/sum;
        };

    };

};

template<>
void BoundShapeFunction<3>::evaluateBoundaryIso(double *Xsi, double *phiB_, double *wpcs, int side, int *inc, std::vector<IParameters_ *> &iparameters, int Npatch){

    int           s = side;    
    std::vector<IParameters_ *> *ipar_;
    ipar_ = &iparameters;

    double xsi[2];

    xsi[0] = Xsi[0];
    xsi[1] = Xsi[1];

    // Front and Back side of the element (plan xy)
    if ((s == 0) || (s == 2)) {
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
        double uBF[3][3]= {}; // stores de base functions in each direction
        double vBF[3][3]= {}; // stores de base functions in each direction
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
        const double xsi1 =((u2 - u1) * xsi[0] + (u2 + u1)) * 0.5;
        const double xsi2 =((v2 - v1) * xsi[1] + (v2 + v1)) * 0.5;

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
            phiB_[i] = phit[i] * wpcs[i]/sum;
        };

    };

    // Left and right (plan vt)
    if ((s == 1) || (s == 3)) {
        int degq = (*ipar_)[Npatch] -> getDegree(2);
        int degn = (*ipar_)[Npatch] -> getDegree(1);
        int npcq = (*ipar_)[Npatch] -> getNcp(2);
        int npcn = (*ipar_)[Npatch] -> getNcp(1);   

        int numLocalBF = (degq+1)*(degn+1);

        double *uknot;
        double *vknot;
        double  uLeft[3];
        double uRight[3];
        double vLeft [3];
        double vRight[3];
        double phiL1 [3];
        double phiL2 [3];
        double phit  [9];
        double uBF[3][3] = {}; // stores de base functions in each direction
        double vBF[3][3] = {}; // stores de base functions in each direction
        double saved,temp;

        uknot = (*ipar_)[Npatch]->gettKnot();
        vknot = (*ipar_)[Npatch]->getvKnot();  

        // indexes from index space where a function related with first local point of element start 
        int uind = inc[2];
        int vind = inc[1]; 

        // parametric coordinates that defines the element
        double u1 = uknot[uind];
        double u2 = uknot[uind+1];
        double v1 = vknot[vind];
        double v2 = vknot[vind +1];

        int degm = degq;
        int npcm = npcq;

        // relates integration space with the parametric one
        const double xsi1 =((u2 - u1) * xsi[0] + (u2 + u1)) * 0.5;
        const double xsi2 =((v2 - v1) * xsi[1] + (v2 + v1)) * 0.5;

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
            phiB_[i] = phit[i] * wpcs[i]/sum;
        };
        
    };

    //Bottom and upper side
    if ((s == 4) || (s == 5)) {
        int degm = (*ipar_)[Npatch] -> getDegree(0);
        int degq = (*ipar_)[Npatch] -> getDegree(2);
        int npcm = (*ipar_)[Npatch] -> getNcp(0);
        int npcq = (*ipar_)[Npatch] -> getNcp(2);   

        int numLocalBF = (degm+1)*(degq+1);

        double *uknot;
        double *vknot;
        double  uLeft[3];
        double uRight[3];
        double vLeft [3];
        double vRight[3];
        double phiL1 [3];
        double phiL2 [3];
        double phit  [9];
        double uBF[3][3] = {}; // stores de base functions in each direction
        double vBF[3][3] = {}; // stores de base functions in each direction
        double saved,temp;

        uknot = (*ipar_)[Npatch]->getuKnot();
        vknot = (*ipar_)[Npatch]->gettKnot();  

        // indexes from index space where a function related with first local point of element start 
        int uind = inc[0];
        int vind = inc[2]; 

        // parametric coordinates that defines the element
        double u1 = uknot[uind];
        double u2 = uknot[uind+1];
        double v1 = vknot[vind];
        double v2 = vknot[vind +1];

        // relates integration space with the parametric one
        const double xsi1 =((u2 - u1) * xsi[0] + (u2 + u1)) * 0.5;
        const double xsi2 =((v2 - v1) * xsi[1] + (v2 + v1)) * 0.5;

        int degn = degq;
        int npcn = npcq;

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
            phiB_[i] = phit[i] * wpcs[i]/sum;
        };
    };
    return;

};


//------------------------------------------------------------------------------
//-------------------COMPUTE SHAPE FUNCTION DERIVATIVE VALUE--------------------
//------------------------------------------------------------------------------
template<>
void BoundShapeFunction<2>::evaluateGradientBoundaryFem(double *Xsi, double **dphiB_) {

    double xsi = Xsi[0];
    double aux;
    double Nnos = 3;
    double AdimCoord[3];
    AdimCoord[0] = -1.;
    AdimCoord[1] =  0.;
    AdimCoord[2] =  1.;


    for (int j = 0; j < Nnos; j++) {
        dphiB_[0][j] = 0.0;
        for (int i = 0; i < Nnos; i++) {
            aux = 1.0;
            if(i != j){                             
                for (int k = 0; k < Nnos; k++) {
                    if ((i != k) && (j != k)) aux = aux*(xsi - AdimCoord[k]);
                };
                dphiB_[0][j] += aux;
            };
        };
    };

    for (int i = 0; i < Nnos; i++) {
        for (int k = 0; k < Nnos; k++) {
            if (i != k) dphiB_[0][i] = dphiB_[0][i] / (AdimCoord[i] - AdimCoord[k]);
        };
    };  


    return;
};

template<>
void BoundShapeFunction<3>::evaluateGradientBoundaryFem(double *Xsi, double **dphiB_) {

    const double xsi1 = Xsi[0];
    const double xsi2 = Xsi[1];
    const double xsi3 = 1. - xsi1 - xsi2;
    
    dphiB_[0][0] = -4. * xsi3 + 1.;
    dphiB_[1][0] = -4. * xsi3 + 1.;

    dphiB_[0][1] = 4. * xsi1 - 1.;
    dphiB_[1][1] = 0.;

    dphiB_[0][2] = 0.;
    dphiB_[1][2] = 4. * xsi2 - 1.;

    dphiB_[0][3] = 4. * (xsi3 - xsi1);
    dphiB_[1][3] = -4. * xsi1;

    dphiB_[0][4] = 4. * xsi2;
    dphiB_[1][4] = 4. * xsi1;

    dphiB_[0][5] = -4. * xsi2;
    dphiB_[1][5] = 4. * (xsi3 - xsi2);


    return;
};

template<>
void BoundShapeFunction<2>::evaluateGradientBoundaryIso(double *Xsi, double **dphiB_, double *wpcs, int side, int *inc, std::vector<IParameters_ *> &iparameters, int Npatch){

    double xsi = Xsi[0];

    int           s = side;    
    //exemplo of type and side
            // side 2
            //7  8  9
    //side 3//4  5  6 //side 1
            //1  2  3
            //side 0

    std::vector<IParameters_ *> *ipar_;  // isogeometric parameters
    ipar_ = &iparameters;

    if ((s == 0) || (s == 2)) {

        // isogemetric data
        int degm = (*ipar_)[Npatch] -> getDegree(0);
        int npcm = (*ipar_)[Npatch] -> getNcp(0);       
        double *uknot;
        double uLeft [3];
        double uRight[3];
        double phiL1 [3];
        double dphiL1[3];
        double aux1[2][2]; // stores the two lines more recently computed
        double uBF[3][3] = {}; // stores de base functions in each direction
        double saved,temp;

        uknot = (*ipar_)[Npatch] -> getuKnot();

        // indexes from index space where a function related with first local point of element start 
        int uind = inc[0];

        // parametric coordinates that defines the element
        double u1 = uknot[uind];
        double u2 = uknot[uind+1];
       
        // relates integration space with the parametric one
        const double xsi1 =((u2 - u1) * xsi + (u2 + u1)) * 0.5;
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

        double sum0 = 0.0;
        for (int i = 0; i < degm + 1 ; i++ ) {
            sum0 += phiL1[i] * wpcs[i];
        };

        double sum = 0.0;
        for (int i = 0; i < degm + 1 ; i++ ) {
            sum += dphiL1[i] * wpcs[i];
        };

        for (int i = 0; i < degm + 1 ; i++ ) {
             dphiB_[0][i] = (sum0 * dphiL1[i] * wpcs[i] - sum * phiL1[i] * wpcs[i])/(sum0*sum0);
        };

    };

    if ((s == 1) || (s == 3)) {

        // isogemetric data
        int degn = (*ipar_)[Npatch] -> getDegree(1);
        int npcn = (*ipar_)[Npatch] -> getNcp(1);   
        double *vknot;
        double vLeft [3];
        double vRight[3];
        double phiL2 [3];   
        double dphiL2[3];
        double aux2[2][2]; // stores the two lines more recently computed
        double vBF[3][3] = {}; // stores de base functions in each direction
        double saved,temp;

        vknot = (*ipar_)[Npatch] -> getvKnot();  

        // indexes from index space where a function related with first local point of element start 
        int vind = inc[1]; 

        // parametric coordinates that defines the element
        double v1 = vknot[vind];
        double v2 = vknot[vind +1];

        // relates integration space with the parametric one
        const double xsi2 =((v2 - v1) * xsi + (v2 + v1)) * 0.5;
        
        // Base function u direction       
        vBF[0][0] = 1.0;

        for (int j = 1; j < (degn + 1); j++) {
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

        for (int i = 0; i < (degn + 1) ; i++ ) {
            phiL2[i] = vBF[i][degn];
        };

        // derivatives u direction;
        int s1,s2,rk,pk,j1,j2,cor,k;
        double d;

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

        double sum0 = 0.0;
        for (int i = 0; i < degn + 1 ; i++ ) {
            sum0 += phiL2[i] * wpcs[i];
        };

        double sum = 0.0;
        for (int i = 0; i < degn + 1 ; i++ ) {
            sum += dphiL2[i] * wpcs[i];
        };

        for (int i = 0; i < degn + 1 ; i++ ) {
             dphiB_[0][i] = (sum0 * dphiL2[i] * wpcs[i] - sum * phiL2[i] * wpcs[i])/(sum0 * sum0);
        };

    };

    return;

};
 
template<>
void BoundShapeFunction<3>::evaluateGradientBoundaryIso(double *Xsi, double **dphiB_, double *wpcs, int side, int *inc, std::vector<IParameters_ *> &iparameters, int Npatch) {

    int           s = side;    
    std::vector<IParameters_ *> *ipar_;
    ipar_ = &iparameters;

    double xsi[2];

    xsi[0] = Xsi[0];
    xsi[1] = Xsi[1];

    // Front and Back side of the element (plan xy)
    if ((s == 0) || (s == 2)) {
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
        double dphiL1[3];
        double phiL2 [3];
        double dphiL2[3];
        double phit  [9];
        double dphit [2][9];
        double aux1[2][3];
        double aux2[2][3];
        double uBF[3][3] = {}; // stores de base functions in each direction
        double vBF[3][3] = {}; // stores de base functions in each direction
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
        const double xsi1 =((u2 - u1) * xsi[0] + (u2 + u1)) * 0.5;
        const double xsi2 =((v2 - v1) * xsi[1] + (v2 + v1)) * 0.5;

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
	        dphiB_[0][i] = ((dphit[0][i] * wpcs[i] * sum - phit[i] * wpcs[i] * sum1) / (sum * sum)); 
	        dphiB_[1][i] = ((dphit[1][i] * wpcs[i] * sum - phit[i] * wpcs[i] * sum2) / (sum * sum)); 
	    };

    };

    // Left and right (plan vt)
    if ((s == 1) || (s == 3)) {
        int degq = (*ipar_)[Npatch] -> getDegree(2);
        int degn = (*ipar_)[Npatch] -> getDegree(1);
        int npcq = (*ipar_)[Npatch] -> getNcp(2);
        int npcn = (*ipar_)[Npatch] -> getNcp(1);   

        int dimt = degq + npcq +1;
        int dimv = degn + npcn +1;
        int numLocalBF = (degq+1)*(degn+1);
        
        double *uknot;
       	double *vknot;
        double uLeft [3];
        double uRight[3];
        double vLeft [3];
        double vRight[3];
        double phiL1 [3];
        double dphiL1[3];
        double phiL2 [3];
        double dphiL2[3];
        double phit  [9];
        double dphit [2][9];
        double aux1[2][3];
        double aux2[2][3];
        double uBF[3][3] = {}; // stores de base functions in each direction
        double vBF[3][3] = {}; // stores de base functions in each direction
        double saved,temp;

        uknot = (*ipar_)[Npatch]->gettKnot();
        vknot = (*ipar_)[Npatch]->getvKnot();  

        // indexes from index space where a function related with first local point of element start 
        int uind = inc[2];
        int vind = inc[1]; 

        // parametric coordinates that defines the element
        double u1 = uknot[uind];
        double u2 = uknot[uind+1];
        double v1 = vknot[vind];
        double v2 = vknot[vind +1];

        int degm = degq;
        int npcm = npcq;

        // relates integration space with the parametric one
        const double xsi1 =((u2 - u1) * xsi[0] + (u2 + u1)) * 0.5;
        const double xsi2 =((v2 - v1) * xsi[1] + (v2 + v1)) * 0.5;

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
	        }
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
	        dphiB_[0][i] = ((dphit[0][i] * wpcs[i] * sum - phit[i] * wpcs[i] * sum1) / (sum * sum)); 
	        dphiB_[1][i] = ((dphit[1][i] * wpcs[i] * sum - phit[i] * wpcs[i] * sum2) / (sum * sum)); 
	    };

    };

    //Bottom and upper side
    if ((s == 4) || (s == 5)) {
        int degm = (*ipar_)[Npatch] -> getDegree(0);
        int degq = (*ipar_)[Npatch] -> getDegree(2);
        int npcm = (*ipar_)[Npatch] -> getNcp(0);
        int npcq = (*ipar_)[Npatch] -> getNcp(2);   

        int numLocalBF = (degm+1)*(degq+1);

        double *uknot;
       	double *vknot;
        double uLeft [3];
        double uRight[3];
        double vLeft [3];
        double vRight[3];
        double phiL1 [3];
        double dphiL1[3];
        double phiL2 [3];
        double dphiL2[3];
        double phit  [9];
        double dphit [2][9];
        double aux1[2][3];
        double aux2[2][3];
        double uBF[3][3] = {}; // stores de base functions in each direction
        double vBF[3][3] = {}; // stores de base functions in each direction
        double saved,temp;

        uknot = (*ipar_)[Npatch]->getuKnot();
        vknot = (*ipar_)[Npatch]->gettKnot();  

        // indexes from index space where a function related with first local point of element start 
        int uind = inc[0];
        int vind = inc[2]; 

        // parametric coordinates that defines the element
        double u1 = uknot[uind];
        double u2 = uknot[uind+1];
        double v1 = vknot[vind];
        double v2 = vknot[vind +1];

        int degn = degq;
        int npcn = npcq;

        // relates integration space with the parametric one
        const double xsi1 =((u2 - u1) * xsi[0] + (u2 + u1)) * 0.5;
        const double xsi2 =((v2 - v1) * xsi[1] + (v2 + v1)) * 0.5;

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
	        dphiB_[0][i] = ((dphit[0][i] * wpcs[i] * sum - phit[i] * wpcs[i] * sum1) / (sum * sum)); 
	        dphiB_[1][i] = ((dphit[1][i] * wpcs[i] * sum - phit[i] * wpcs[i] * sum2) / (sum * sum)); 
	    };
    };

};


template<>
void BoundShapeFunction<3>::evaluateHessianBoundaryFem(double ***ddphiB_) {

     ddphiB_[0][0][0] = 4.;
     ddphiB_[0][1][0] = 4.;
     ddphiB_[1][0][0] = 4.;
     ddphiB_[1][1][0] = 4.;

     ddphiB_[0][0][1] = 4.;
     ddphiB_[0][1][1] = 0.;
     ddphiB_[1][0][1] = 0.;
     ddphiB_[1][1][1] = 0.;

     ddphiB_[0][0][2] = 0.;
     ddphiB_[0][1][2] = 0.;
     ddphiB_[1][0][2] = 0.;
     ddphiB_[1][1][2] = 4.;

     ddphiB_[0][0][3] = -8.;
     ddphiB_[0][1][3] = -4.;
     ddphiB_[1][0][3] = -4.;
     ddphiB_[1][1][3] = 0.;

     ddphiB_[0][0][4] = 0.;
     ddphiB_[0][1][4] = 4.;
     ddphiB_[1][0][4] = 4.;
     ddphiB_[1][1][4] = 0.;

     ddphiB_[0][0][5] = 0.;
     ddphiB_[0][1][5] = -4.;
     ddphiB_[1][0][5] = -4.;
     ddphiB_[1][1][5] = -8.;


    return;
};