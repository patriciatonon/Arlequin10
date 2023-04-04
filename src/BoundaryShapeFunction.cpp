
#include "BoundaryShapeFunction.h"

//------------------------------------------------------------------------------
//-------------------------COMPUTE SHAPE FUNCTION VALUE-------------------------
//------------------------------------------------------------------------------
template<int DIM>
void BoundShapeFunction<DIM>::basisFunctionsBS(double &xsi, int &deg, double *knot, int &inc, double *phiL){

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
void BoundShapeFunction<DIM>::derBasisFunctionsBS(int &degree, double &xsi, int &deg, double *knot, int &inc, double *dphiL){

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

    double xsi = Xsi[0];

    std::vector<IParameters_ *> *ipar_ = &iparameters;

    if ((s == 0) || (s == 2)) {

        // isogemetric data
        int degm = (*ipar_)[Npatch] -> getDegree(0);
        double *uknot = (*ipar_)[Npatch] -> getuKnot();
        double phiL1 [3];

        // indexes from index space where a function related with first local point of element start 
        int uind = inc[0];

        basisFunctionsBS(xsi, degm, uknot, uind, phiL1);

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
        double *vknot = (*ipar_)[Npatch] -> getvKnot(); 
        double phiL2 [3];

        // indexes from index space where a function related with first local point of element start 
        int vind = inc[1]; 

        basisFunctionsBS(xsi, degn, vknot, vind, phiL2);

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

        double *uknot = (*ipar_)[Npatch]->getuKnot();
        double *vknot = (*ipar_)[Npatch]->getvKnot(); 
        double phiL1 [3];
        double phiL2 [3];
        double phit  [9];

        // indexes from index space where a function related with first local point of element start 
        int uind = inc[0];
        int vind = inc[1]; 

        basisFunctionsBS(xsi[0], degm, uknot, uind, phiL1);
        basisFunctionsBS(xsi[1], degn, vknot, vind, phiL2);

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
        int degm = (*ipar_)[Npatch] -> getDegree(2);
        int degn = (*ipar_)[Npatch] -> getDegree(1);
        int npcm = (*ipar_)[Npatch] -> getNcp(2);
        int npcn = (*ipar_)[Npatch] -> getNcp(1);   

        int numLocalBF = (degm+1)*(degn+1);

        double phiL1 [3];
        double phiL2 [3];
        double phit  [9];


         double *uknot = (*ipar_)[Npatch]->gettKnot();
         double *vknot = (*ipar_)[Npatch]->getvKnot();  

        // indexes from index space where a function related with first local point of element start 
        int uind = inc[2];
        int vind = inc[1]; 

        basisFunctionsBS(xsi[0], degm, uknot, uind, phiL1);
        basisFunctionsBS(xsi[1], degn, vknot, vind, phiL2);
        
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
        int degn = (*ipar_)[Npatch] -> getDegree(2);
        int npcm = (*ipar_)[Npatch] -> getNcp(0);
        int npcn = (*ipar_)[Npatch] -> getNcp(2);   

        int numLocalBF = (degm+1)*(degn+1);

        double phiL1 [3];
        double phiL2 [3];
        double phit  [9];

        double *uknot = (*ipar_)[Npatch]->getuKnot();
        double *vknot = (*ipar_)[Npatch]->gettKnot();  

        // indexes from index space where a function related with first local point of element start 
        int uind = inc[0];
        int vind = inc[2]; 

        basisFunctionsBS(xsi[0], degm, uknot, uind, phiL1);
        basisFunctionsBS(xsi[1], degn, vknot, vind, phiL2);

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
        double phiL1 [3];
        double dphiL1[3];
        double *uknot = (*ipar_)[Npatch] -> getuKnot();
        // indexes from index space where a function related with first local point of element start 
        int uind = inc[0];

        basisFunctionsBS(xsi,degm,uknot,uind,phiL1);
        int degree = 1;
        derBasisFunctionsBS(degree,xsi,degm,uknot,uind,dphiL1);


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
        double phiL2 [3];   
        double dphiL2[3];
        double *vknot = (*ipar_)[Npatch] -> getvKnot();  
        // indexes from index space where a function related with first local point of element start 
        int vind = inc[1]; 

        basisFunctionsBS(xsi,degn,vknot,vind,phiL2);
        int degree = 1;
        derBasisFunctionsBS(degree,xsi,degn,vknot,vind,dphiL2);

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
        int numLocalBF = (degm+1)*(degn+1);
        double phiL1 [3];
        double dphiL1[3];
        double phiL2 [3];
        double dphiL2[3];
        double phit  [9];
        double dphit [2][9];

        double *uknot = (*ipar_)[Npatch]->getuKnot();
        double *vknot = (*ipar_)[Npatch]->getvKnot();  

        int uind = inc[0];
        int vind = inc[1]; 
        
        basisFunctionsBS(xsi[0],degm,uknot,uind,phiL1);
        int degree = 1;
        derBasisFunctionsBS(degree,xsi[0],degm,uknot,uind,dphiL1);

        basisFunctionsBS(xsi[1],degn,vknot,vind,phiL2);
        derBasisFunctionsBS(degree,xsi[1],degn,vknot,vind,dphiL2);
    

        int index = 0;
        for (int j = 0; j < degn + 1 ; j ++) {
            for (int k = 0; k < degm + 1; k ++) {
                phit[index] = phiL1[k] * phiL2 [j];
                index ++;
            };
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
        int degm = (*ipar_)[Npatch] -> getDegree(2);
        int degn = (*ipar_)[Npatch] -> getDegree(1);
        int numLocalBF = (degm+1)*(degn+1);
        double phiL1 [3];
        double dphiL1[3];
        double phiL2 [3];
        double dphiL2[3];
        double phit  [9];
        double dphit [2][9];
        double *uknot = (*ipar_)[Npatch]->gettKnot();
        double *vknot = (*ipar_)[Npatch]->getvKnot();  
        // indexes from index space where a function related with first local point of element start 
        int uind = inc[2];
        int vind = inc[1]; 

        basisFunctionsBS(xsi[0],degm,uknot,uind,phiL1);
        int degree = 1;
        derBasisFunctionsBS(degree,xsi[0],degm,uknot,uind,dphiL1);

        basisFunctionsBS(xsi[1],degn,vknot,vind,phiL2);
        derBasisFunctionsBS(degree,xsi[1],degn,vknot,vind,dphiL2);
        
        int index = 0;
        for (int j = 0; j < degn + 1 ; j ++) {
            for (int k = 0; k < degm + 1; k ++) {
                phit[index] = phiL1[k] * phiL2 [j];
                index ++;
            };
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
        int degn = (*ipar_)[Npatch] -> getDegree(2);
        int numLocalBF = (degm+1)*(degn+1);
        double phiL1 [3];
        double dphiL1[3];
        double phiL2 [3];
        double dphiL2[3];
        double phit  [9];
        double dphit [2][9];
        double *uknot = (*ipar_)[Npatch]->getuKnot();
        double *vknot = (*ipar_)[Npatch]->gettKnot();  
        // indexes from index space where a function related with first local point of element start 
        int uind = inc[0];
        int vind = inc[2]; 

        basisFunctionsBS(xsi[0],degm,uknot,uind,phiL1);
        int degree = 1;
        derBasisFunctionsBS(degree,xsi[0],degm,uknot,uind,dphiL1);

        basisFunctionsBS(xsi[1],degn,vknot,vind,phiL2);
        derBasisFunctionsBS(degree,xsi[1],degn,vknot,vind,dphiL2);
        
        int index = 0;
        for (int j = 0; j < degn + 1 ; j ++) {
            for (int k = 0; k < degm + 1; k ++) {
                phit[index] = phiL1[k] * phiL2 [j];
                index ++;
            };
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

template<int DIM>
void BoundShapeFunction<DIM>::evaluateHessianBoundaryIso(double *Xsi, double ***ddphiB_, double *wpcs, int side, int *inc, std::vector<IParameters_ *> &iparameters, int Npatch) {

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
        int numLocalBF = (degm+1)*(degn+1);
        double phiL1 [3];
        double dphiL1[3];
        double ddphiL1[3];
        double phiL2 [3];
        double dphiL2[3];
        double ddphiL2[3];
        double phit  [9];
        double dphit [2][9];
        double ddphit [2][2][9];

        double *uknot = (*ipar_)[Npatch]->getuKnot();
        double *vknot = (*ipar_)[Npatch]->getvKnot();  

        int uind = inc[0];
        int vind = inc[1]; 
        
        basisFunctionsBS(xsi[0],degm,uknot,uind,phiL1);
        int degree = 1;
        derBasisFunctionsBS(degree,xsi[0],degm,uknot,uind,dphiL1);
        degree = 2;
        derBasisFunctionsBS(degree,xsi[0],degm,uknot,uind,ddphiL1);

        basisFunctionsBS(xsi[1],degn,vknot,vind,phiL2);
        degree = 1;
        derBasisFunctionsBS(degree,xsi[1],degn,vknot,vind,dphiL2);
        degree = 2;
        derBasisFunctionsBS(degree,xsi[1],degn,vknot,vind,ddphiL2);
	    

        int index = 0;
        for (int j = 0; j< degn; j ++) {
            for (int k = 0; k < degm ; k ++) {
            
                phit[index] = phiL1[k] * phiL2[j];

                dphit[0][index] = dphiL1[k] * phiL2[j];
                dphit[1][index] = phiL1[k] * dphiL2[j];
                
                ddphit[0][0][index] = ddphiL1[k] * phiL2[j]; 
                ddphit[0][1][index] = dphiL1[k] * dphiL2[j];
                ddphit[1][0][index] = dphiL1[k] * dphiL2[j];
                ddphit[1][1][index] = phiL1[k] * ddphiL2[j];

                index ++;
            };
        };

    
        double sumF = 0.0;
        double sumDer[2] = {};
        double sumSecDer[2][2] = {};

        for (int i = 0; i < numLocalBF; i++) {
            sumF += phit[i] * wpcs[i];
            for (int j = 0; j < 2; j++){
                sumDer[j] += dphit[j][i] * wpcs[i];
                for (int k = 0; k < 2; k++){
                    sumSecDer[j][k] += ddphit[j][k][i] * wpcs[i];
                };
            };
        };

        double v_ = sumF * sumF;
        for (int i = 0; i < numLocalBF; i++) {
            for (int j = 0; j < 2; j++){
                double u_ = dphit[j][i] * wpcs[i] * sumF - phit[i] * wpcs[i] * sumDer[j];
                for (int k = 0; k < 2; k++){
                    double ul_ = ddphit[j][k][i] * wpcs[i] * sumF + dphit[j][i] * wpcs[i] * sumDer[k] - 
                                dphit[k][i] * wpcs[i] * sumDer[j] - phit[i] * wpcs[i] * sumSecDer[j][k];;
                    double vl_ = 2 * sumF * sumDer[k];
                    ddphiB_[j][k][i] = ((ul_ * v_ - u_ * vl_)/ (v_ * v_));
                };
            };
        };

    };

    // Left and right (plan vt)
    if ((s == 1) || (s == 3)) {
        int degm = (*ipar_)[Npatch] -> getDegree(2);
        int degn = (*ipar_)[Npatch] -> getDegree(1);
        int numLocalBF = (degm+1)*(degn+1);
        double phiL1 [3];
        double dphiL1[3];
        double ddphiL1[3];
        double phiL2 [3];
        double dphiL2[3];
        double ddphiL2[3];
        double phit  [9];
        double dphit [2][9];
        double ddphit [2][2][9];
        double *uknot = (*ipar_)[Npatch]->gettKnot();
        double *vknot = (*ipar_)[Npatch]->getvKnot();  
        // indexes from index space where a function related with first local point of element start 
        int uind = inc[2];
        int vind = inc[1]; 

                basisFunctionsBS(xsi[0],degm,uknot,uind,phiL1);
        int degree = 1;
        derBasisFunctionsBS(degree,xsi[0],degm,uknot,uind,dphiL1);
        degree = 2;
        derBasisFunctionsBS(degree,xsi[0],degm,uknot,uind,ddphiL1);

        basisFunctionsBS(xsi[1],degn,vknot,vind,phiL2);
        degree = 1;
        derBasisFunctionsBS(degree,xsi[1],degn,vknot,vind,dphiL2);
        degree = 2;
        derBasisFunctionsBS(degree,xsi[1],degn,vknot,vind,ddphiL2);
	    

        int index = 0;
        for (int j = 0; j< degn; j ++) {
            for (int k = 0; k < degm ; k ++) {
            
                phit[index] = phiL1[k] * phiL2[j];

                dphit[0][index] = dphiL1[k] * phiL2[j];
                dphit[1][index] = phiL1[k] * dphiL2[j];
                
                ddphit[0][0][index] = ddphiL1[k] * phiL2[j]; 
                ddphit[0][1][index] = dphiL1[k] * dphiL2[j];
                ddphit[1][0][index] = dphiL1[k] * dphiL2[j];
                ddphit[1][1][index] = phiL1[k] * ddphiL2[j];

                index ++;
            };
        };

    
        double sumF = 0.0;
        double sumDer[2] = {};
        double sumSecDer[2][2] = {};

        for (int i = 0; i < numLocalBF; i++) {
            sumF += phit[i] * wpcs[i];
            for (int j = 0; j < 2; j++){
                sumDer[j] += dphit[j][i] * wpcs[i];
                for (int k = 0; k < 2; k++){
                    sumSecDer[j][k] += ddphit[j][k][i] * wpcs[i];
                };
            };
        };

        double v_ = sumF * sumF;
        for (int i = 0; i < numLocalBF; i++) {
            for (int j = 0; j < 2; j++){
                double u_ = dphit[j][i] * wpcs[i] * sumF - phit[i] * wpcs[i] * sumDer[j];
                for (int k = 0; k < 2; k++){
                    double ul_ = ddphit[j][k][i] * wpcs[i] * sumF + dphit[j][i] * wpcs[i] * sumDer[k] - 
                                dphit[k][i] * wpcs[i] * sumDer[j] - phit[i] * wpcs[i] * sumSecDer[j][k];;
                    double vl_ = 2 * sumF * sumDer[k];
                    ddphiB_[j][k][i] = ((ul_ * v_ - u_ * vl_)/ (v_ * v_));
                };
            };
        };
    };

    //Bottom and upper side
    if ((s == 4) || (s == 5)) {
        int degm = (*ipar_)[Npatch] -> getDegree(0);
        int degn = (*ipar_)[Npatch] -> getDegree(2);
        int numLocalBF = (degm+1)*(degn+1);
        double phiL1 [3];
        double dphiL1[3];
        double ddphiL1[3];
        double phiL2 [3];
        double dphiL2[3];
        double ddphiL2[3];
        double phit  [9];
        double dphit [2][9];
        double ddphit [2][2][9];
        double *uknot = (*ipar_)[Npatch]->getuKnot();
        double *vknot = (*ipar_)[Npatch]->gettKnot();  
        // indexes from index space where a function related with first local point of element start 
        int uind = inc[0];
        int vind = inc[2]; 

        basisFunctionsBS(xsi[0],degm,uknot,uind,phiL1);
        int degree = 1;
        derBasisFunctionsBS(degree,xsi[0],degm,uknot,uind,dphiL1);
        degree = 2;
        derBasisFunctionsBS(degree,xsi[0],degm,uknot,uind,ddphiL1);

        basisFunctionsBS(xsi[1],degn,vknot,vind,phiL2);
        degree = 1;
        derBasisFunctionsBS(degree,xsi[1],degn,vknot,vind,dphiL2);
        degree = 2;
        derBasisFunctionsBS(degree,xsi[1],degn,vknot,vind,ddphiL2);
	    

        int index = 0;
        for (int j = 0; j< degn; j ++) {
            for (int k = 0; k < degm ; k ++) {
            
                phit[index] = phiL1[k] * phiL2[j];

                dphit[0][index] = dphiL1[k] * phiL2[j];
                dphit[1][index] = phiL1[k] * dphiL2[j];
                
                ddphit[0][0][index] = ddphiL1[k] * phiL2[j]; 
                ddphit[0][1][index] = dphiL1[k] * dphiL2[j];
                ddphit[1][0][index] = dphiL1[k] * dphiL2[j];
                ddphit[1][1][index] = phiL1[k] * ddphiL2[j];

                index ++;
            };
        };

    
        double sumF = 0.0;
        double sumDer[2] = {};
        double sumSecDer[2][2] = {};

        for (int i = 0; i < numLocalBF; i++) {
            sumF += phit[i] * wpcs[i];
            for (int j = 0; j < 2; j++){
                sumDer[j] += dphit[j][i] * wpcs[i];
                for (int k = 0; k < 2; k++){
                    sumSecDer[j][k] += ddphit[j][k][i] * wpcs[i];
                };
            };
        };

        double v_ = sumF * sumF;
        for (int i = 0; i < numLocalBF; i++) {
            for (int j = 0; j < 2; j++){
                double u_ = dphit[j][i] * wpcs[i] * sumF - phit[i] * wpcs[i] * sumDer[j];
                for (int k = 0; k < 2; k++){
                    double ul_ = ddphit[j][k][i] * wpcs[i] * sumF + dphit[j][i] * wpcs[i] * sumDer[k] - 
                                dphit[k][i] * wpcs[i] * sumDer[j] - phit[i] * wpcs[i] * sumSecDer[j][k];;
                    double vl_ = 2 * sumF * sumDer[k];
                    ddphiB_[j][k][i] = ((ul_ * v_ - u_ * vl_)/ (v_ * v_));
                };
            };
        };
    };

};
template class BoundShapeFunction<2>;
template class BoundShapeFunction<3>;
