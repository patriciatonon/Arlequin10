//------------------------------------------------------------------------------
// 
//                      Patricia Tonon and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2019 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------ISOGEMETRIC CONSTANT PARAMETERS-----------------------
//------------------------------------------------------------------------------

#ifndef ISOGEOMETRIC_PARAMETERS_H
#define ISOGEOMETRIC_PARAMETERS_H

using namespace boost::numeric;

/// Defines the fluid boundary shape functions

template<int DIM>
class IsogeometricParameters{

private:
    double *uknot_;                  // Knot vector u direction 
    double *vknot_;                  // Knot vector v direction 
    double *tknot_;                  // Knot vector t direction 
    int     deg_[DIM];
    int     ncp_[DIM];     
    int     ipatch_;                 // Index of the patch
    int     numElemPatchx;           // Number of elements in x direction per patch
    int     numElemPatchy;           // Number of elements in y direction per patch
    int     numElemPatchz;           // Number of elements in z direction per patch
    
public:
    IsogeometricParameters(int ipatch, int *deg, int *ncp){
        
        ipatch_ = ipatch;
        for (int i=0; i< DIM; i++){
            deg_[i] = deg[i];
            ncp_[i] = ncp[i];
        }

    };

    // Sets the number of elements in each patch - x direction
    void setNumElemPatchx(int numElem) {numElemPatchx = numElem;};

    // Gets the number of elements in each patch - x direction
    int getNumElemPatchx(){return numElemPatchx;};

    // Sets the number of elements in each patch - y direction
    void setNumElemPatchy(int numElem) {numElemPatchy = numElem;};

    // Gets the number of elements in each patch - y direction
    int getNumElemPatchy(){return numElemPatchy;};

    // Sets the number of elements in each patch - z direction
    void setNumElemPatchz(int numElem) {numElemPatchz = numElem;};

    // Gets the number of elements in each patch - z direction
    int getNumElemPatchz(){return numElemPatchz;};

    // Sets the knot vector - u direction
    void setuKnot(double *uknot){uknot_ = uknot;};

    // Sets the knot vector - v direction
    void setvKnot(double *vknot){vknot_ = vknot;};

    // Sets the knot vector - w direction
    void settKnot(double *tknot){tknot_ = tknot;};

    // Gets the knot vector - u direction
    double* getuKnot(){return uknot_;};

    // Gets the knot vector - v direction
    double* getvKnot(){return vknot_;};

    // Gets the knot vector - t direction
    double* gettKnot(){return tknot_;};

    // Gets the degree of the functions in each direction per patch
    int getDegree(int dir) {return deg_[dir];}

    // Gets the basis functions number in each direction per patch
    int getNcp(int dir) {return ncp_[dir];}

};


#endif
