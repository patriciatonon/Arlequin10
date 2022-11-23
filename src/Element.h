//------------------------------------------------------------------------------
// 
//        Jeferson W D Fernandes,Rodolfo A K Sanches and Patrícia Tonon
//                             University of Sao Paulo
//                           (C) 2019 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------------------ELEMENT------------------------------------
//-----------------------------------------------------------------------------

#ifndef ELEMENT_H
#define ELEMENT_H

#include <boost/numeric/ublas/vector.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include "Node.h"
#include "BoundaryIntegrationQuadrature.h"
#include "BoundaryIntegrationQuadratureIso.h"
#include "IntegrationQuadrature.h"
#include "SpecialIntegrationQuadrature.h"
#include "FluidParameters.h"
#include "IsogeometricParameters.h"

using namespace boost::numeric;

// Defines the fluid element object and all the element information
template<int DIM>
class Element{
 
public:
    // Defines the classes locally:
    //Class of nodes
    typedef Node<DIM>                              Nodes;
    //Class of normal quadrature
    typedef IntegQuadrature<DIM>                   NormalQuad;
    //Class of special quadrature
    typedef SpecialIntegQuadrature<DIM>            SpecialQuad;
    //Class of Fluid parameters
    typedef FluidParameters<DIM>                   FParameters;   
    //Class of isogeometric parameters
    typedef IsogeometricParameters<DIM>            IParameters;

private:
    
    //Element basic information
	FParameters                *parameters;         // Fluid parameters
    std::vector<IParameters *> *iparameters;        // Isogeometric parameters for each patch
	std::vector<Nodes *>       *nodes_;    			// Nodes
    int*                       connect_;            // Mesh connectivity
    int*                       BeConnect_;          // Bezier mesh connectivity
    int            			   index_;              // Element index
    int                        ElemType_;           // Element type: 0 - FEM elements or 1 - IGA elements
    int            			   Npatch_;             // Number of the patch
    std::vector<int>           elemSide_;           // Elements sides in boundary          
    std::vector<int>           elemBound_;          // Number of the boundary
    double                     xK[DIM], XK[DIM];    // Element boxes coordinates

    //Coarse mesh information when integration is performed in fine mesh
    std::vector<Nodes *>       *nodesC_;            //  Coarse nodes
    std::vector<IParameters *> *iparametersC;       //  Isogeometric parameters for each patch
    int*                       connectC_;           //  Coarse element connectivity 
    int                        NpatchC_;            //  Coarse patch
    
    //Information for analysis in the gluing zone
    bool                       model;               // Model: true - fine mesh/ false - coarse mesh
    bool                       glueZone;            // Defines if the element is in the gluing zone (true = gluing zone)  
    

    // Integration points and element in the coarse mesh correspondent to the fine mesh in gluing zone 
    // double intPointCorrespXsi_ISO[25][DIM];		       
    // double intPointCorrespElem_ISO[25];  
    double intPointCorrespXsi_ISO[18*DIM-27][DIM];            
    double intPointCorrespElem_ISO[18*DIM-27]; 
    
    double intPointCorrespXsi_FEM[8*DIM-9][DIM];            
    double intPointCorrespElem_FEM[8*DIM-9]; 
    // double intPointCorrespXsi_FEM[12][DIM];            
    // double intPointCorrespElem_FEM[12]; 
    
    // Integration point energy weight 
    double intPointWeightFunction_ISO[18*DIM-27];
    double intPointWeightFunctionPrev_ISO[18*DIM-27];
    double intPointWeightFunction_FEM[8*DIM-9];
    double intPointWeightFunctionPrev_FEM[8*DIM-9];
    

    // Integration point energy weight gluing zone
    double intPointWeightFunctionSpecial_ISO[18*DIM-27];
    double intPointWeightFunctionSpecialPrev_ISO[18*DIM-27];
    // double intPointWeightFunctionSpecial_ISO[25];
    // double intPointWeightFunctionSpecialPrev_ISO[25];


    double intPointWeightFunctionSpecial_FEM[8*DIM-9];
    double intPointWeightFunctionSpecialPrev_FEM[8*DIM-9];
    // double intPointWeightFunctionSpecial_FEM[12];
    // double intPointWeightFunctionSpecialPrev_FEM[12];

    int iTimeStep;
    double integScheme;

    
public:
    // fluid element constructor
    Element(int index, int *connect, std::vector<Nodes *> &nodes, int ElemType, FParameters &param, std::vector<IParameters *> &iparam, int ipatch){

        index_ = index;
        connect_ = connect;
        nodes_ = &nodes;
        ElemType_ = ElemType;
        parameters = &param;
        iparameters = &iparam;
        Npatch_ = ipatch;
        glueZone = false;

    	for (int i = 0; i < 8*DIM-9; i++){
    		intPointWeightFunction_FEM[i] = 1.;
    		intPointWeightFunctionPrev_FEM[i] = 1.;
        };
    	for (int i = 0; i < 18*DIM-27; i++){
    		intPointWeightFunction_ISO[i] = 1.;
    		intPointWeightFunctionPrev_ISO[i] = 1.;
        };


        for (int i = 0; i < 8*DIM-9; i++){
            intPointWeightFunctionSpecial_FEM[i] = 1.;
            intPointWeightFunctionSpecialPrev_FEM[i] = 1.;
            for (int j = 0; j < DIM; j++) intPointCorrespXsi_FEM[i][j] = 0.0;
            intPointCorrespElem_FEM[i] = 0;
        };

        // for (int i = 0; i < 12; i++){
        // //     intPointWeightFunctionSpecial_FEM[i] = 1.;
        // //     intPointWeightFunctionSpecialPrev_FEM[i] = 1.;
        //     for (int j = 0; j < DIM; j++) intPointCorrespXsi_FEM[i][j] = 0.0;
        //     intPointCorrespElem_FEM[i] = 0;
        // };
        
        for (int i = 0; i < 18*DIM-27; i++){
            intPointWeightFunctionSpecial_ISO[i] = 1.;
            intPointWeightFunctionSpecialPrev_ISO[i] = 1.;
            for (int j = 0; j < DIM; j++) intPointCorrespXsi_ISO[i][j] = 0.0;
            intPointCorrespElem_ISO[i] = 0;
        };

        //  for (int i = 0; i < 25; i++){
        //     intPointWeightFunctionSpecial_ISO[i] = 1.;
        //     intPointWeightFunctionSpecialPrev_ISO[i] = 1.;
        //     for (int j = 0; j < DIM; j++) intPointCorrespXsi_ISO[i][j] = 0.0;
        //     intPointCorrespElem_ISO[i] = 0;
        // };

        iTimeStep = 0;
          
    };


    //........................Elements basic information.........................
    // Sets the element connectivity
    void setConnectivity(int *connect){connect_ = connect;};
    // Gets the element connectivity
    int* getConnectivity(){return connect_;};
    // Returns the element index
    int getIndex(){return index_;};
    // Returns the element type
    int getElemType(){return ElemType_;};


    //........................Aditional information for IGA elements.............
    // Returns the element patch
    int getPatch(){return Npatch_;};
    //Sets the Bezier element connectivity
    void setBezierConnectivity(int* Becon_){BeConnect_ = Becon_;};
    //Gets the Bezier element connectivity
    int* getBezierConnectivity(){return BeConnect_;};

    
    //...........................Boundary Information...........................
    // Sets the element side in boundary and the correspondent number boundary
    void setElemSideInBoundary(std::vector<int> &side, std::vector<int> &boundary){
         elemSide_ = side;
         elemBound_ = boundary;
    };
    // Gets the element side/number of boundary
    std::pair <std::vector<int>,std::vector<int> > getElemSideInBoundary(){
         return std::make_pair(elemSide_,elemBound_);
    };


    //.....................Intersection element Information......................
    // Gets the element intersection parameters 
    void setIntersectionParameters(double *x, double *X){
        for (int i=0; i<DIM; ++i) {
            xK[i] = x[i];
            XK[i] = X[i];
        }
    };
    // Gets the coordinates intersection parameters
    std::pair<double*,double*> getXIntersectionParameter() {return std::make_pair(xK,XK);};

    
    //........................Arlequin zone informations......................
    // Sets which model the fluid element belongs (true = fine mesh; false = coarse mesh)
    void setModel(bool m){model = m;};
    // Sets if the element is in the gluing zone (true = gluing zone)
    void setGlueZone(){glueZone = true;};
    // Gets if the element is in the gluing zone 
    bool getGlueZone(){return glueZone;};


    //......................Integration Points Information......................
    
    // Returns the number of integration points of the quadrature 
    int getNumberOfIntegrationPoints_FEM(){NormalQuad nQuad = NormalQuad(); return (nQuad.endFem() - nQuad.beginFem());}; //outside of the gluing zone for FEM elements
    int getNumberOfIntegrationPoints_ISO(){NormalQuad nQuad = NormalQuad(); return (nQuad.endIso() - nQuad.beginIso());}; //outside of the gluing zone for IGA elements
    int getNumberOfIntegrationPointsSpecial_ISO(){SpecialQuad sQuad = SpecialQuad(); return (sQuad.endIso() - sQuad.beginIso());}; //inside of gluing zone
   	int getNumberOfIntegrationPointsSpecial_FEM(){SpecialQuad sQuad = SpecialQuad(); return (sQuad.endFem() - sQuad.beginFem());}; //inside of gluing zone

    // Set the integration point and element correspondence in the coarse mesh
    void setIntegrationPointCorrespondence_ISO(int ipoint, double *x, int elem){
        intPointCorrespElem_ISO[ipoint] = elem;
        intPointCorrespXsi_ISO[ipoint][0] = x[0];
        intPointCorrespXsi_ISO[ipoint][1] = x[1]; 
    };
    
    // Returns correspondent element in the coarse mesh
    int getIntegPointCorrespondenceElement_ISO(int index){
    	return intPointCorrespElem_ISO[index];
    };

    // Returns the correspondent integration point coordinates in the coarse mesh
    double getIntegPointCoordinatesValue_ISO(int index, int dir){
        double x = intPointCorrespXsi_ISO[index][dir];
        return x;
    };

        // Set the integration point and element correspondence in the coarse mesh
    void setIntegrationPointCorrespondence_FEM(int ipoint, double *x, int elem){
        intPointCorrespElem_FEM[ipoint] = elem;
        intPointCorrespXsi_FEM[ipoint][0] = x[0];
        intPointCorrespXsi_FEM[ipoint][1] = x[1]; 
    };
    
    // Returns correspondent element in the coarse mesh
    int getIntegPointCorrespondenceElement_FEM(int index){
        return intPointCorrespElem_FEM[index];
    };

    // Returns the correspondent integration point coordinates in the coarse mesh
    double getIntegPointCoordinatesValue_FEM(int index, int dir){
        double x = intPointCorrespXsi_FEM[index][dir];
        return x;
    };

    // Compute the weighted integration points
    void setIntegPointWeightFunction_ISO();
    void setIntegPointWeightFunction_FEM();

    //......................Jacobian Matrixes and Derivatives....................
    // Compute and store the spatial jacobian matrix
    void getJacobianMatrix_FEM(double &djac_, double *xsi, double **Jac, double **ainv_);
    void getJacobianMatrix_COARSE_FEM(double *xsi, double **Jac, double **ainv_);
    void getJacobianMatrix_ISO(double &djac_, double *xsi, double **quadJacMat, double **ainv_);
    void getJacobianMatrix_COARSE_ISO(double *xsi, double **quadJacMat, double **ainv_);
    // Compute and store the quadrature jacobian matrix
    void getQuadJacobianMatrix_ISO(double *xsi, double **quadJacMatInv);


    void getJacobianMatrixValues_FEM(double *xsi, double **ainv_){

        double **Jac;
        Jac= new double*[DIM];
        for (int i = 0; i < DIM; ++i) Jac[i] = new double[DIM];
        
        double djac_;
        getJacobianMatrix_FEM(djac_, xsi,Jac,ainv_);

        for (int i = 0; i < DIM; ++i) delete [] Jac[i];
        delete [] Jac;
    };

    void getJacobianMatrixValues_ISO(double *xsi, double **ainv_ ){
        
        double **quadJacMat;
        quadJacMat = new double*[DIM];
        for (int i = 0; i < DIM; ++i) quadJacMat[i] = new double[DIM];
        
        double djac_;
        getJacobianMatrix_ISO(djac_, xsi,quadJacMat,ainv_);

    	for (int i = 0; i < DIM; ++i) delete [] quadJacMat[i];
    	delete [] quadJacMat;
	}


    double getDetJac_ISO(double *xsi){
        
        double **ainv_;
        ainv_ = new double*[DIM];
        for (int i = 0; i < DIM; ++i) ainv_[i] = new double[DIM];

        double **quadJacMat;
        quadJacMat = new double*[DIM];
        for (int i = 0; i < DIM; ++i) quadJacMat[i] = new double[DIM];

        double djac_;
        getJacobianMatrix_ISO(djac_, xsi,quadJacMat,ainv_);

    	for (int i = 0; i < DIM; ++i) delete [] ainv_[i];
    	delete [] ainv_;

    	for (int i = 0; i < DIM; ++i) delete [] quadJacMat[i];
    	delete [] quadJacMat;

        return djac_;
    }



    // Compute and store the shape function spatial derivatives
    void getSpatialDerivatives_FEM(double *xsi, double **ainv_, double **dphi_dx);
    void getSpatialDerivatives_ISO(double *xsi, double **ainv_, double **dphi_dx);
    void getSpatialDerivatives_COARSE_ISO(double *xsi, double **ainv_, double **dphi_dx);

    void getSecondSpatialDerivatives_FEM(double **ainv_, double ***ddphi_dx);
    void getSecondSpatialDerivatives_ISO(double *xsi, double **ainv_, double ***ddphi_dx);
    void getSecondSpatialDerivatives_COARSE_ISO(double *xsi, double **ainv_, double ***ddphi_dx);


    // Compute and stores the interpolated variables
    void getInterpCoord(int &na, double *phi_, double *x_, double *xPrevious_);
    void getInterpCoord_ISO(int &na, double *phi_, double *x_, double *xPrev_);
    void getInterpVel(int &na, double *phi_, double *u_, double *uPrev_);
    void getInterpVelDer(int &na, double **dphi_dx, double **du_dx, double **duPrev_dx);
    void getInterpVelCoarse(int &na, double *phi_, double *u_, double *uPrev_);
    void getInterpVelDerCoarse(int &na, double **dphi_dx, double **du_dx, double **duPrev_dx); 
    void getInterpSecondVelDer(int &na, double ***dphi_dx, double ***ddu_dxdx, double ***dduPrev_dxdx);
    void getInterpSecondVelDerCoarse(int &na, double ***dphi_dx, double ***ddu_dxdx, double ***dduPrev_dxdx);
    void getInterpMeshVel(int &na, double *phi_, double *uMesh_, double *uMeshPrev_);
    void getInterpMeshVelDer(int &na, double **dphi_dx, double **duMesh_dx, double **duMeshPrev_dx);
    void getInterpMeshVelCoarse(int &na, double *phi_, double *uMesh_, double *uMeshPrev_);
    void getInterpMeshVelDerCoarse(int &na, double **dphi_dx, double **duMesh_dx, double **duMeshPrev_dx);
    void getInterpAccel(int &na, double *phi_, double *accel_, double *accelPrev_); 
    void getInterpAccelDer(int &na, double **dphi_dx, double **daccel_dx, double **daccelPrev_dx);
    void getInterpAccelDerCoarse(int &na, double **dphi_dx, double **daccel_dx, double **daccelPrev_dx);
    void getInterpPress(int &na, double *phi_, double &press_);
    void getInterpPressDer(int &na, double **dphi_dx, double *dpress_dx);
    void getInterpSecondPressDer(int &na, double ***ddphi_dx, double **ddpress_dxdx);
    void getInterpSecondPressDerCoarse(int &na, double ***ddphi_dx, double **ddpress_dxdx);
    void getInterpLambda(int &na, double *phi_, double *lambda_);
    void getInterpLambdaDer(int &na, double **dphi_dx, double **dlambda_dx);

    //......................Stabilization Parameters....................
    // Compute and store the SUPG, PSPG and LSIC stabilization parameters
    void getNewParameterSUPG_FEM(double &tSUPG_, double &tPSPG_, double &tLSIC_, double **Jac , double *phi, double **dphi_dx);
    void getNewParameterSUPG_ISO(double &tSUPG_, double &tPSPG_, double &tLSIC_, double **quadJacMat, double *phi_, double **dphi_dx);
    void getParameterArlequin_ISO(double &tARLQ_, double *phi_, double **dphi_dx);
    void getParameterArlequin_FEM(double &tARLQ_, double *phi_, double **dphi_dx);
    void getParameterArlequin_COARSE_ISO(double &tARLQ_, double *phi_, double **dphi_dx);
    void getParameterArlequin2_FEM(double &djac_, double &weight_, double &tARLQ_, 
                                   double *phi_, double *phiC_,double **dphi_dx, double ***ddphi_dx);
    void getParameterArlequin2_COARSE_ISO(double &djac_, double &weight_, double &tARLQ_, 
                                          double *phi_, double *phiC_, double **dphi_dx, double **dphiC_dx, 
                                          double ***ddphi_dx, double ***ddphiC_dx);


 	//......................Drag and Lift Parameters....................
    // Compute and store the drag and lift forces at the element boundary
    void computeDragAndLiftForces_ISO(int &side, double &pDForce, double &pLForce, double &fDForce, double &fLForce, 
                                      double &dForce, double &lForce, double &aux_Mom, double &aux_Per);
    void computeDragAndLiftForces_FEM(int &side, double &pDForce, double &pLForce, double &fDForce, double &fLForce, 
                                      double &dForce, double &lForce, double &aux_Mom, double &aux_Per);

   
    //.......................Element vectors and matrices.......................
    // Compute and store the element matrix for the incompressible flow problem
    void getElemMatrix(int &na, double &djac_, double &weight_, double &tSUPG_, double &tPSPG_, double &tLSIC_,
                       int &index, double *phi_, double **dphi_dx, double **jacobianNRMatrix);  
    
    //Compute and store the residual vector for the incompressible flow problem
    void getResidualVector(int &na, double &djac_, double &weight_, double &tSUPG_, double &tPSPG_, double &tLSIC_,
                          int &index, double *phi_, double **dphi_dx, double *rhsVector);

    
    //Arlequin Matrixes and Vectores
    void getMatrixAndVectorsSameMesh(int &na, double &djac_, double &weight_,double *phi_, double **dphi_dx, 
                                     double **lagrMultMatrix, double *rhsVector1, double *rhsVector2);
    void getMatrixAndVectorsSameMesh_tSUPG_tPSPG(int &na, double &djac_, double &weight_,double &tSUPG_, double &tPSPG_,double *phi_, 
                                                 double **dphi_dx,double **jacobianNRMatrix,double *rhsVector);
    void getMatrixAndVectorsSameMeshArlqStab(int &na, double &djac_, double &weight_,double &tARLQ_, int &index,
                                             double *phi_, double **dphi_dx, double ***ddphi_dx,
                                             double **arlequinStabD, double *arlequinStabVectorD,
                                             double **arlequinStab1, double *arlequinStabVector1);
   
    void getMatrixAndVectorsDifferentMesh(double &djac_, double &weight_,int &na, double *phi_, 
                                          double **dphi_dx, int &naC, double *phiC_, double **dphiC_dx,
                                          double **lagrMultMatrix, double *rhsVector1, double *rhsVector2);   
    void getMatrixAndVectorsDifferentMeshTemp(double &djac_, double &weight_,double &tARLQ_, int &na, double *phi_, 
                                              double **dphi_dx, int &naC, double *phiC_, double **dphiC_dx,
                                              double **lagrMultMatrix, double *rhsVector1, double *rhsVector2,
                                              double **arlequinStab, double *arlequinStabVector);  
    void getMatrixAndVectorsDifferentMesh_tSUPG_tPSPG(double &djac_, double &weight_,double &tSUPG_, double &tPSPG_,
                                                      int &na, double *phi_, int &naC, double *phiC_, double **dphiC_dx,
                                                      double **jacobianNRMatrix,double *rhsVector);
    void getMatrixAndVectorsDifferentMeshArlqStab(double &djac_, double &weight_,double &tARLQ_,int &index, 
                                                  int &na, double **dphi_dx, int &naC, double *phiC_, double **dphiC_dx,double ***ddphiC_dx,
                                                  double **arlequinStabD, double *arlequinStabVectorD,
                                                  double **arlequinStab0, double *arlequinStabVector0);


    //...............................Problem type...............................
    // Compute the Transient Navier-Stokes problem matrices and vectors
    void getTransientNavierStokes_FEM(double **jacobianNRMatrix, double *rhsVector);
    void getTransientNavierStokes_ISO(double **jacobianNRMatrix, double *rhsVector);

    // Compute and store the Lagrange multiplier operator when integrating the same mesh portion
    void getLagrangeMultipliersSameMesh_FEM(int &index, double **lagrMultMatrix, 
                                            double *lagrMultVector, double *rhsVector);
    void getLagrangeMultipliersSameMeshArlqStab_FEM(int &index, int &ipatchC, std::vector<Nodes *> &nodesCoarse_,int *connecC,
                                                    std::vector<IParameters *> &iparamC, double **arlequinStabD, 
                                                    double *arlequinStabVectorD, double **arlequinStab1, double *arlequinStabVector1);
    void getLagrangeMultipliersSameMesh_tSUPG_tPSPG_FEM(int &index, double **jacobianNRMatrix, double *rhsVector);
    
    
    void getLagrangeMultipliersSameMesh_ISO(double **lagrMultMatrix, double *lagrMultVector, double *rhsVector);
    void getLagrangeMultipliersSameMesh_tSUPG_tPSPG_ISO(double **jacobianNRMatrix, double *rhsVector);
    void getLagrangeMultipliersSameMeshArlqStab_ISO(double **arlequinStabD, double *arlequinStabVectorD,
                                                    double **arlequinStab1, double *arlequinStabVector1);


    //Compute and store the Lagrange multiplier operator when integrationg the different mesh portion
    void getLagrangeMultipliersDifferentMesh_ISO_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_,int *connecC,
                                                    std::vector<IParameters *> &iparamC, int &ielem, 
                                                    double **lagrMultMatrix,double *rhsVector1, double *rhsVector2,
                                                    double **arlequinStab, double *arlequinStabVector);
    void getLagrangeMultipliersDifferentMesh_FEM_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_,int *connecC,
                                                    std::vector<IParameters *> &iparamC, int &ielem, 
                                                    double **lagrMultMatrix,double *rhsVector1, double *rhsVector2);
    void getLagrangeMultipliersDifferentMesh_tSUPG_tPSPG_FEM_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_,int *connecC,
                                                    			std::vector<IParameters *> &iparamC, int &ielem, 
                                                    		    double **jacobianNRMatrix, double *rhsVector);
    void getLagrangeMultipliersDifferentMeshArlqStab_FEM_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_,int *connecC,
                                                             std::vector<IParameters *> &iparamC, int &ielem, 
                                                             double **arlequinStabD, double *arlequinStabVectorD,
                                                             double **arlequinStab0, double *arlequinStabVector0);
    
    void getLagrangeMultipliersDifferentMesh_ISO_FEM(std::vector<Nodes *> &nodesCoarse_,int *connecC,int &ielem, 
                                                    double **lagrMultMatrix,double *rhsVector1, double *rhsVector2);
    void getLagrangeMultipliersDifferentMesh_tSUPG_tPSPG_ISO_FEM(std::vector<Nodes *> &nodesCoarse_,int *connecC,int &ielem, 
                                                                double **jacobianNRMatrix, double *rhsVector);
    void getLagrangeMultipliersDifferentMeshArlqStab_ISO_FEM(std::vector<Nodes *> &nodesCoarse_,int *connecC,int &ielem, 
                                                             double **arlequinStabD, double *arlequinStabVectorD,
                                                             double **arlequinStab0, double *arlequinStabVector0);

    void getLagrangeMultipliersDifferentMesh_FEM_FEM(std::vector<Nodes *> &nodesCoarse_,int *connecC,int &ielem, 
                                                     double **lagrMultMatrix,double *rhsVector1, double *rhsVector2,
                                                     double **arlequinStab, double *arlequinStabVector);


    void setBoundaryConditions_ISO(double **jacobianNRMatrix, double *rhsVector);

    //.......................Bezier Element transformation.......................
    // Computes Bézier extract operator 
    void getMatrixC(double **MatrixC);
    //Computes inverse Bézier extractor operator    
    void getInvMatrixC(double **MatrixCuInv, double **MatrixCvInv);

    FParameters* getFluidParameters(){
        return parameters;
    }

    void setTimeStep(int timeStep) {
    	iTimeStep = timeStep;
    }

};


#endif

