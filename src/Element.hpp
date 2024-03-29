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

#include <boost/algorithm/minmax_element.hpp>
#include "Node.hpp"
#include "BoundaryIntegrationQuadrature.hpp"
#include "BoundaryIntegrationQuadratureIso.hpp"
#include "IntegrationQuadrature.hpp"
#include "SpecialIntegrationQuadrature.hpp"
#include "FluidParameters.hpp"
#include "IsogeometricParameters.hpp"

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
	FParameters                *parameters;          // Fluid parameters
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
    
    //Analysis variables
    double         tSUPG_,tPSPG_,tLSIC_;			// Stabilization parameters
    double         tARLQ_;
    double         weight_;							// Weight for numerical integration
    double         djac_;                           // Jacobian determinant
    double         x_, y_, z_;	           			// Interpolated coordinates
    double         u_, v_, w_, p_;         			// Interpolated velocity and pressure
    double         uprev_, vprev_, wprev_; 			// Interpolated previous velocity 
    double         ax_, ay_, az_;          			// Interpolated acceleration
    double         axprev_, ayprev_, azprev_;       // Interpolated previous acceleration
    double         dax_dx,dax_dy,day_dx,day_dy;     // Interpolated acceleration derivative
    double         daxprev_dx, daxprev_dy,
                   dayprev_dx, dayprev_dy;          // Interpolated previous acceleration derivative
    double         du_dx, du_dy, du_dz, 
    			   dv_dx, dv_dy, dv_dz, 
                   dw_dx, dw_dy, dw_dz;    			// Interpolated fluid velocity spatial derivatives             
    double         ddu_dxdx, ddu_dxdy,
                   ddu_dydx, ddu_dydy,
                   ddv_dxdx, ddv_dxdy,
                   ddv_dydx, ddv_dydy;              // Interpolated fluid spatial velocity second derivatives
    double         duprev_dx, duprev_dy, duprev_dz,
                   dvprev_dx, dvprev_dy, dvprev_dz,
                   dwprev_dx, dwprev_dy, dwprev_dz; // Interpolated previous fluid velocity spatial derivatives
    double         dduprev_dxdx, dduprev_dxdy,
                   dduprev_dydx, dduprev_dydy, 
                   ddvprev_dxdx, ddvprev_dxdy,
                   ddvprev_dydx, ddvprev_dydy;      // Interpolated previous fluid spatial velocity second derivatives
    double         dp_dx, dp_dy, dp_dz;				// Interpolated pressure spatial derivative 
    double         ddp_dxdx, ddp_dxdy, ddp_dydx, 
                   ddp_dydy;                        // Interpolated pressure spatial second derivatives   
    double         lamx_, lamy_;					// Lagrange multipliers
    double         lamx_dx, lamx_dy, 			    // Lagrange multipliers derivatives
    			   lamy_dx, lamy_dy;      
    double         umesh_, vmesh_, wmesh_; 			// Interpolated mesh velocity
    double         umeshprev_, vmeshprev_,
				   wmeshprev_; 			   			// Interpolated previous mesh velocity
	double         dumesh_dx, dumesh_dy, dumesh_dz, 
    			   dvmesh_dx, dvmesh_dy, dvmesh_dz, 
                   dwmesh_dx, dwmesh_dy, dwmesh_dz; // Interpolated mesh velocity spatial derivatives  			   
    double         dumeshprev_dx, dumeshprev_dy, 
    			   dumeshprev_dz,
                   dvmeshprev_dx, dvmeshprev_dy,
                   dvmeshprev_dz,
                   dwmeshprev_dx, dwmeshprev_dy, 
                   dwmeshprev_dz; 					// Interpolated previous mesh velocity spatial derivatives

    //Aerodynamic coeficients     
    double         pressureDragForce;               // Drag force from pressure
    double         pressureLiftForce;               // Lift force from pressure
    double         frictionDragForce;               // Drag force from friction
    double         frictionLiftForce;               // Lift force from friction
    double         dragForce;                       // Total drag force
    double         liftForce;                       // Total lifth force
    double         pitchingMoment;                  // Moment 
    double         perimeter;                       // Perimeter where is calculated the Drag and Lift forces
    

    // Integration points and element in the coarse mesh correspondent to the fine mesh in gluing zone 
    // double intPointCorrespXsi_ISO[25][DIM];		       
    // double intPointCorrespElem_ISO[25];
    // double intPointCorrespXsi_ISO[64][DIM];         
    // double intPointCorrespElem_ISO[64];  
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
    // double intPointWeightFunctionSpecial_ISO[18*DIM-27];
    // double intPointWeightFunctionSpecialPrev_ISO[18*DIM-27];
    // double intPointWeightFunctionSpecial_ISO[25];
    // double intPointWeightFunctionSpecialPrev_ISO[25];
    // double intPointWeightFunctionSpecial_ISO[64];
    // double intPointWeightFunctionSpecialPrev_ISO[64];

    // double intPointWeightFunctionSpecial_FEM[8*DIM-9];
    // // double intPointWeightFunctionSpecialPrev_FEM[8*DIM-9];
    // double intPointWeightFunctionSpecial_FEM[12];
    // double intPointWeightFunctionSpecialPrev_FEM[12];

    int iTimeStep;

    
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
        //     intPointWeightFunctionSpecial_FEM[i] = 1.;
        //     intPointWeightFunctionSpecialPrev_FEM[i] = 1.;
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
        //     intPointWeightFunctionSpecial_ISO[i] = 1.;
        //     intPointWeightFunctionSpecialPrev_ISO[i] = 1.;
            for (int j = 0; j < DIM; j++) intPointCorrespXsi_ISO[i][j] = 0.0;
            intPointCorrespElem_ISO[i] = 0;
        };

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
    void getJacobianMatrix_FEM(double *xsi, double **Jac, double **ainv_);
    void getJacobianMatrix_COARSE_FEM(double *xsi, double **Jac, double **ainv_);
    void getJacobianMatrix_ISO(double *xsi, double **quadJacMat, double **ainv_);
    void getJacobianMatrix_COARSE_ISO(double *xsi, double **quadJacMat, double **ainv_);
    // Compute and store the quadrature jacobian matrix
    void getQuadJacobianMatrix_ISO(double *xsi, double **quadJacMatInv);


    void getJacobianMatrixValues_FEM(double *xsi, double **ainv_){

        double **Jac;
        Jac= new double*[DIM];
        for (int i = 0; i < DIM; ++i) Jac[i] = new double[DIM];

        getJacobianMatrix_FEM(xsi,Jac,ainv_);

        for (int i = 0; i < DIM; ++i) delete [] Jac[i];
        delete [] Jac;
    };

    void getJacobianMatrixValues_ISO(double *xsi, double **ainv_ ){
        
        double **quadJacMat;
        quadJacMat = new double*[DIM];
        for (int i = 0; i < DIM; ++i) quadJacMat[i] = new double[DIM];
        
        getJacobianMatrix_ISO(xsi,quadJacMat,ainv_);

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

        getJacobianMatrix_ISO(xsi,quadJacMat,ainv_);

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

    void getSecondSpatialDerivatives_FEM(double *xsi, double **ainv_, double ***ddphi_dx);
    void getSecondSpatialDerivatives_ISO(double *xsi, double **ainv_, double ***ddphi_dx);
    void getSecondSpatialDerivatives_COARSE_ISO(double *xsi, double **ainv_, double ***ddphi_dx);


    // Compute and stores the interpolated variables
    void getVelAndDerivatives_FEM(double *phi_, double **dphi_dx);
    void getVelAndDerivatives_ISO(double *phi_, double **dphi_dx);

    // Compute and stores the interpolated variables when Arlequin Problem
    void getInterpolatedVariablesSameMesh_FEM(double *phi_, double **dphi_dx);
    void getInterpolatedVariablesSameMeshArlqStab_FEM(double *phi_, double **dphi_dx, double ***ddphi_dx);
    void getInterpolatedVariablesSameMesh_ISO(double *phi_, double **dphi_dx);
    void getInterpolatedVariablesSameMeshArlqStab_ISO(double *phi_, double **dphi_dx, double ***ddphi_dx);
    void getInterpolatedVariablesDifferentMesh_FEM_ISO(double *phi_, double **dphi_dx,double *phiC_, double **dphiC_dx);
    void getInterpolatedVariablesDifferentMesh_ISO_FEM(double *phi_, double **dphi_dx,double *phiC_, double **dphiC_dx);
    void getInterpolatedVariablesDifferentMesh_FEM_FEM(double *phi_, double **dphi_dx,double *phiC_, double **dphiC_dx);
    void getInterpolatedVariablesDifferentMesh_ISO_ISO(double *phi_, double **dphi_dx,double *phiC_, double **dphiC_dx);
    void getInterpolatedVariablesDifferentMeshArlqStab_FEM_ISO(double **dphi_dx,double *phiC_, double **dphiC_dx, double ***ddphiC_dx);
    void getInterpolatedVariablesDifferentMeshArlqStab_ISO_FEM(double **dphi_dx,double *phiC_, double **dphiC_dx, double ***ddphiC_dx);
    void getInterpolateLagMultiplierDerivatives(double **dphi_dx);


    //......................Stabilization Parameters....................
    // Compute and store the SUPG, PSPG and LSIC stabilization parameters
    void getNewParameterSUPG_FEM(double **Jac , double *phi, double **dphi_dx);
    void getNewParameterSUPG_ISO(double **quadJacMat, double *phi_, double **dphi_dx);
    void getParameterArlequin_ISO(double *phi_, double **dphi_dx);
    void getParameterArlequin_FEM(double *phi_, double **dphi_dx);
    void getParameterArlequin2_FEM(double *phi_, double *phiC_, double **dphi_dx, double ***ddphi_dx);


 	//......................Drag and Lift Parameters....................
    // Compute and store the drag and lift forces at the element boundary
    void computeDragAndLiftForces_ISO(int &side, double &pDForce, double &pLForce, double &fDForce, double &fLForce, 
                                      double &dForce, double &lForce, double &aux_Mom, double &aux_Per);
    void computeDragAndLiftForces_FEM(int &side, double &pDForce, double &pLForce, double &fDForce, double &fLForce, 
                                      double &dForce, double &lForce, double &aux_Mom, double &aux_Per);

   
    //.......................Element vectors and matrices.......................
    // Compute and store the element matrix for the incompressible flow problem
    void getElemMatrix_FEM(int &index,double *phi_, double **dphi_dx, double **jacobianNRMatrix);
    void getElemMatrix_ISO(int &index,double *phi_, double **dphi_dx, double **jacobianNRMatrix);   
    //Compute and store the residual vector for the incompressible flow problem
    void getResidualVector_FEM(int &index,double *phi_, double **dphi_dx, double *rhsVector);
    void getResidualVector_ISO(int &index,double *phi_, double **dphi_dx, double *rhsVector);
    
    //Arlequin Matrixes and Vectores
    void getMatrixAndVectorsSameMesh_FEM(double *phi_, double **dphi_dx, double **lagrMultMatrix, 
                                         double *rhsVector1, double *rhsVector2);
    void getMatrixAndVectorsSameMesh_tSUPG_tPSPG_FEM(double *phi_, double **dphi_dx,double **jacobianNRMatrix,
    											    double *rhsVector);
    void getMatrixAndVectorsSameMeshArlqStab_FEM(int &index, double *phi_, double **dphi_dx, double ***ddphi_dx,
                                                 double **arlequinStabD, double *arlequinStabVectorD,
                                                 double **arlequinStab1, double *arlequinStabVector1);
    void getMatrixAndVectorsSameMesh_ISO(double *phi_, double **dphi_dx, double **lagrMultMatrix, 
        	                             double *rhsVector1, double *rhsVector2); 
    void getMatrixAndVectorsSameMesh_tSUPG_tPSPG_ISO(double *phi_, double **dphi_dx,double **jacobianNRMatrix,
                                                    double *rhsVector);
    void getMatrixAndVectorsSameMeshArlqStab_ISO(int &index, double *phi_, double **dphi_dx, double ***ddphi_dx,
                                                 double **arlequinStabD, double *arlequinStabVectorD,
                                                 double **arlequinStab1, double *arlequinStabVector1);
    void getMatrixAndVectorsDifferentMesh_FEM_FEM(double *phi_, double **dphi_dx, double *phiC_, double **dphiC_dx,
                                                  double **lagrMultMatrix, double *rhsVector1, double *rhsVector2,
                                                  double **arlequinStab, double *arlequinStabVector);   
    void getMatrixAndVectorsDifferentMesh_FEM_ISO(double *phi_, double **dphi_dx, double *phiC_, double **dphiC_dx,
                                                  double **lagrMultMatrix, double *rhsVector1, double *rhsVector2);
    void getMatrixAndVectorsDifferentMesh_ISO_FEM(double *phi_, double **dphi_dx, double *phiC_, double **dphiC_dx,
                                                  double **lagrMultMatrix, double *rhsVector1, double *rhsVector2);
    void getMatrixAndVectorsDifferentMesh_ISO_ISO(double *phi_, double **dphi_dx, double *phiC_, double **dphiC_dx,
        									      double **lagrMultMatrix, double *rhsVector1, double *rhsVector2,
                                                  double **arlequinStab, double *arlequinStabVector);
    void getMatrixAndVectorsDifferentMesh_tSUPG_tPSPG_FEM_ISO(double *phi_,double **dphiC_dx,
                                                              double **jacobianNRMatrix,double *rhsVector);
    void getMatrixAndVectorsDifferentMesh_tSUPG_tPSPG_ISO_FEM(double *phi_,double **dphiC_dx,
                                                              double **jacobianNRMatrix,double *rhsVector);
    void getMatrixAndVectorsDifferentMeshArlqStab_FEM_ISO(int &index, double **dphi_dx, double *phiC_, double **dphiC_dx,double ***ddphiC_dx,
                                                          double **arlequinStabD, double *arlequinStabVectorD,
                                                          double **arlequinStab0, double *arlequinStabVector0);
    void getMatrixAndVectorsDifferentMeshArlqStab_ISO_FEM(int &index, double **dphi_dx, double *phiC_, double **dphiC_dx,double ***ddphiC_dx,
                                                          double **arlequinStabD, double *arlequinStabVectorD,
                                                          double **arlequinStab0, double *arlequinStabVector0);

    //...............................Problem type...............................
    // Compute the Transient Navier-Stokes problem matrices and vectors
    void getTransientNavierStokes_FEM(double **jacobianNRMatrix, double *rhsVector);
    void getTransientNavierStokes_ISO(double **jacobianNRMatrix, double *rhsVector);

    // Compute and store the Lagrange multiplier operator when integrating the same mesh portion
    void getLagrangeMultipliersSameMesh_FEM(int &index, double **lagrMultMatrix, double *lagrMultVector, double *rhsVector);
    // void getLagrangeMultipliersSameMesh_FEM(double **lagrMultMatrix, double *lagrMultVector, double *rhsVector);
    void getLagrangeMultipliersSameMeshArlqStab_FEM(int &index, int &ipatchC, std::vector<Nodes *> &nodesCoarse_,int *connecC,
                                                    std::vector<IParameters *> &iparamC, double **arlequinStabD, 
                                                    double *arlequinStabVectorD,double **arlequinStab1, double *arlequinStabVector1);
    // void getLagrangeMultipliersSameMeshArlqStab_FEM(double **arlequinStabD, double *arlequinStabVectorD,
    //                                                 double **arlequinStab1, double *arlequinStabVector1);
    void getLagrangeMultipliersSameMesh_tSUPG_tPSPG_FEM(double **jacobianNRMatrix, double *rhsVector);
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

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------WEIGHTED INTEGRATION POINTS-------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::setIntegPointWeightFunction_FEM(){


	int dim = 2;
    QuadShapeFunction<2>    shapeQuad;

    // NORMAL QUADRATURE
    NormalQuad              nQuad = NormalQuad(); 
    int index = 0;
    for(double* it = nQuad.beginFem(); it != nQuad.endFem(); it++){
    	intPointWeightFunctionPrev_FEM[index] = intPointWeightFunction_FEM[index];
    	intPointWeightFunction_FEM[index] = 0.;
    	index++;
    };
   	
    double xsi[dim],phi_[6];
    index = 0;
    for(double* it = nQuad.beginFem(); it != nQuad.endFem(); it++){
        
        for (int i = 0; i < dim; i++) xsi[i] = nQuad.PointListFem(index,i);
        shapeQuad.evaluateFem(xsi,phi_);

        for (int i = 0; i < 6; i++){
        	intPointWeightFunction_FEM[index] += (*nodes_)[connect_[i]] -> getWeightFunction() * phi_[i];
        };
    	
    	index++;
    };


    //SPECIAL QUADRATURE
    // SpecialQuad sQuad = SpecialQuad();
    // index = 0;
    // for(double* it = sQuad.beginFem(); it != sQuad.endFem(); it++){
    // 	intPointWeightFunctionSpecialPrev_FEM[index] = intPointWeightFunctionSpecial_FEM[index];
    // 	intPointWeightFunctionSpecial_FEM[index] = 0.;
    // 	index++;
    // };

    // index = 0;
    // for(double* it = sQuad.beginFem(); it != sQuad.endFem(); it++){
    	
    // 	for (int i = 0; i < dim; i++) xsi[i] = sQuad.PointListFem(index,i);
    //     shapeQuad.evaluateFem(xsi,phi_);

    //     for (int i = 0; i < 6; i++){
    //     	intPointWeightFunctionSpecial_FEM[index] += (*nodes_)[connect_[i]] -> getWeightFunction() * phi_[i];
    //     };
    	
    // 	index++;
    // };
    
};

template<>
void Element<2>::setIntegPointWeightFunction_ISO(){

    
	int dim = 2;
    QuadShapeFunction<2>    shapeQuad;

    // NORMAL QUADRATURE
    NormalQuad              nQuad = NormalQuad(); 
    int index = 0;
    for(double* it = nQuad.beginIso(); it != nQuad.endIso(); it++){
    	intPointWeightFunctionPrev_ISO[index] = intPointWeightFunction_ISO[index];
    	intPointWeightFunction_ISO[index] = 0.;
    	index++;
    };
   	
   	//data for evaluation of NURBS functions
    double wpc[9],phi_[9],xsi[dim];
    for (int i = 0; i < 9; i ++) wpc[i] = (*nodes_)[connect_[i]] -> getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]] -> getINC();  

    index = 0;
    for(double* it = nQuad.beginIso(); it != nQuad.endIso(); it++){
        
        for (int i = 0; i < dim; i++) xsi[i] = nQuad.PointListIso(index,i);
        shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,(*iparameters),Npatch_);

        for (int i = 0; i < 9; i++){
        	intPointWeightFunction_ISO[index] += (*nodes_)[connect_[i]] -> getWeightFunction() * phi_[i];
        };
    	
    	index++;
    };

    //SPECIAL QUADRATURE
    // SpecialQuad sQuad = SpecialQuad();
    // index = 0;
    // for(double* it = sQuad.beginIso(); it != sQuad.endIso(); it++){
    // 	intPointWeightFunctionSpecialPrev_ISO[index] = intPointWeightFunctionSpecial_ISO[index];
    // 	intPointWeightFunctionSpecial_ISO[index] = 0.;
    // 	index++;
    // };

    // index = 0;
    // for(double* it = sQuad.beginIso(); it != sQuad.endIso(); it++){
    	
    // 	for (int i = 0; i < dim; i++) xsi[i] = sQuad.PointListIso(index,i);
    //     shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,(*iparameters),Npatch_);

    //     for (int i = 0; i < 9; i++){
    //     	intPointWeightFunctionSpecial_ISO[index] += (*nodes_)[connect_[i]] -> getWeightFunction() * phi_[i];
    //     };
    	
    // 	index++;
    // };


};




//------------------------------------------------------------------------------
//-------------------------SPATIAL TRANSFORM - JACOBIAN-------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getJacobianMatrix_FEM(double *xsi,  double **Jac,  double **ainv_) {

    //Computes the spatial Jacobian matrix and its inverse
    QuadShapeFunction<2>    shapeQuad;
    
    double **dphi;
    dphi = new double*[2];
    for (int i = 0; i < 2; ++i) dphi[i] = new double[6];
    
    shapeQuad.evaluateGradientFem(xsi,dphi);

    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;

    double &alpha_f = parameters -> getAlphaF();

    for (int i=0; i<6; i++){
		
		double xna_ = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(0) + 
                      (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(0);
        double yna_ = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(1) + 
                      (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(1);    

        dx_dxsi1 += xna_ * dphi[0][i];
        dx_dxsi2 += xna_ * dphi[1][i];
        dy_dxsi1 += yna_ * dphi[0][i];
        dy_dxsi2 += yna_ * dphi[1][i];   
    };

    // Defining the Jacobian matrix
    Jac[0][0] = dx_dxsi1;
    Jac[0][1] = dx_dxsi2;
    Jac[1][0] = dy_dxsi1;
    Jac[1][1] = dy_dxsi2;

    //Computing the jacobian determinant
    djac_ = Jac[0][0] * Jac[1][1] - Jac[0][1] * Jac[1][0];

    //Computing Jacobian inverse (transposed)
    ainv_[0][0] =  dy_dxsi2 / djac_;
    ainv_[0][1] = -dy_dxsi1 / djac_;
    ainv_[1][0] = -dx_dxsi2 / djac_;
    ainv_[1][1] =  dx_dxsi1 / djac_;

    for (int i = 0; i < 2; ++i) delete [] dphi[i];
    delete [] dphi;

    return;
    
};

template<>
void Element<2>::getJacobianMatrix_COARSE_FEM(double *xsi,  double **Jac,  double **ainv_) {

    //Computes the spatial Jacobian matrix and its inverse
    QuadShapeFunction<2>    shapeQuad;
    
    double **dphi;
    dphi = new double*[2];
    for (int i = 0; i < 2; ++i) dphi[i] = new double[6];
    
    shapeQuad.evaluateGradientFem(xsi,dphi);

    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;

    double &alpha_f = parameters->getAlphaF();

    for (int i=0; i<6; i++){
        
        double xna_ = alpha_f * (*nodesC_)[connectC_[i]] -> getCoordinateValue(0) + 
                      (1. - alpha_f) * (*nodesC_)[connectC_[i]] -> getPreviousCoordinateValue(0);
        double yna_ = alpha_f * (*nodesC_)[connectC_[i]] -> getCoordinateValue(1) + 
                      (1. - alpha_f) * (*nodesC_)[connectC_[i]] -> getPreviousCoordinateValue(1);    

        dx_dxsi1 += xna_ * dphi[0][i];
        dx_dxsi2 += xna_ * dphi[1][i];
        dy_dxsi1 += yna_ * dphi[0][i];
        dy_dxsi2 += yna_ * dphi[1][i];   
    };

    // Defining the Jacobian matrix
    Jac[0][0] = dx_dxsi1;
    Jac[0][1] = dx_dxsi2;
    Jac[1][0] = dy_dxsi1;
    Jac[1][1] = dy_dxsi2;

    //Computing the jacobian determinant
    double djacC_ = Jac[0][0] * Jac[1][1] - Jac[0][1] * Jac[1][0];

    //Computing Jacobian inverse (transposed)
    ainv_[0][0] =  dy_dxsi2 / djacC_;
    ainv_[0][1] = -dy_dxsi1 / djacC_;
    ainv_[1][0] = -dx_dxsi2 / djacC_;
    ainv_[1][1] =  dx_dxsi1 / djacC_;

    for (int i = 0; i < 2; ++i) delete [] dphi[i];
    delete [] dphi;

    return;
    
};


template<>
void Element<2>::getJacobianMatrix_ISO(double *xsi, double **quadJacMat, double **ainv_) {

    QuadShapeFunction<2>    shapeQuad;

    double &alpha_f = parameters->getAlphaF();
   
    double **dphi;
    dphi = new double*[2];
    for (int i = 0; i < 2; ++i) dphi[i] = new double[9];

    int *inc_= (*nodes_)[connect_[8]] -> getINC();
    double wpc[9];
    for (int i = 0; i<9; i++) wpc[i] = (*nodes_)[connect_[i]] -> getWeightPC();
    
    shapeQuad.evaluateGradientIso(xsi,dphi,wpc,inc_,(*iparameters),Npatch_);
    
    //Computes the Jacobian matrix - dx/dxsi
    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;
    
    for (int i=0; i<9; i++){
        
        double xna_ = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(0) + 
                      (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(0);
        double yna_ = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(1) + 
                      (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(1);

        dx_dxsi1 += xna_ * dphi[0][i]/wpc[i];
        dx_dxsi2 += xna_ * dphi[1][i]/wpc[i];
        dy_dxsi1 += yna_ * dphi[0][i]/wpc[i];
        dy_dxsi2 += yna_ * dphi[1][i]/wpc[i];   
    };

    //Computing the jacobian determinant 
    double djacpp = dx_dxsi1 * dy_dxsi2 - dx_dxsi2 * dy_dxsi1;

    //Computing jacobian inverse matrix (transposed) - dxsi/dx
    ainv_[0][0] =  dy_dxsi2 / djacpp;
    ainv_[0][1] = -dy_dxsi1 / djacpp;
    ainv_[1][0] = -dx_dxsi2 / djacpp;
    ainv_[1][1] =  dx_dxsi1 / djacpp;

    // Computing the quadrature jacobian matrix - dx/dxqsi
    int degm = (*iparameters)[Npatch_] -> getDegree(0);
    int degn = (*iparameters)[Npatch_] -> getDegree(1);
    int npcm = (*iparameters)[Npatch_] -> getNcp(0);
    int npcn = (*iparameters)[Npatch_] -> getNcp(1);   
    double *uknot_ = (*iparameters)[Npatch_] -> getuKnot();
    double *vknot_ = (*iparameters)[Npatch_] -> getvKnot();  

    int uind = inc_[0];
    int vind = inc_[1]; 

    double u1 = uknot_[uind];
    double u2 = uknot_[uind+1];
    double v1 = vknot_[vind];
    double v2 = vknot_[vind +1];

    double dxsi1_dqxsi1 = (u2-u1)*0.5;
    double dxsi1_dqxsi2 = 0.;
    double dxsi2_dqxsi1 = 0.;
    double dxsi2_dqxsi2 = (v2-v1)*0.5;

    quadJacMat [0][0] = dxsi1_dqxsi1 * dx_dxsi1 + dxsi2_dqxsi1 * dx_dxsi2;
    quadJacMat [0][1] = dxsi1_dqxsi2 * dx_dxsi1 + dxsi2_dqxsi2 * dx_dxsi2;
    quadJacMat [1][0] = dxsi1_dqxsi1 * dy_dxsi1 + dxsi2_dqxsi1 * dy_dxsi2;
    quadJacMat [1][1] = dxsi1_dqxsi2 * dy_dxsi1 + dxsi2_dqxsi2 * dy_dxsi2;

    djac_ = quadJacMat[0][0] * quadJacMat [1][1] - quadJacMat [0][1] * quadJacMat [1][0];

    for (int i = 0; i < 2; ++i) delete [] dphi[i];
    delete [] dphi;

    return;
    
};

template<>
void Element<2>::getJacobianMatrix_COARSE_ISO(double *xsi, double **quadJacMat, double **ainv_) {

    QuadShapeFunction<2>    shapeQuad;


    double &alpha_f = parameters->getAlphaF();
   
    double **dphi;
    dphi = new double*[2];
    for (int i = 0; i < 2; ++i) dphi[i] = new double[9];

    int *inc_= (*nodesC_)[connectC_[8]] -> getINC();
    double wpc[9];
    for (int i = 0; i<9; i++) wpc[i] = (*nodesC_)[connectC_[i]] -> getWeightPC();
    
    shapeQuad.evaluateGradientIso(xsi,dphi,wpc,inc_,(*iparametersC),NpatchC_);
    
    //Computes the Jacobian matrix - dx/dxsi
    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;
    
    for (int i=0; i<9; i++){
        
        double xna_ = alpha_f * (*nodesC_)[connectC_[i]] -> getCoordinateValue(0) + 
                      (1. - alpha_f) * (*nodesC_)[connectC_[i]] -> getPreviousCoordinateValue(0);
        double yna_ = alpha_f * (*nodesC_)[connectC_[i]]-> getCoordinateValue(1) + 
                      (1. - alpha_f) * (*nodesC_)[connectC_[i]] -> getPreviousCoordinateValue(1);

        dx_dxsi1 += xna_ * dphi[0][i]/wpc[i];
        dx_dxsi2 += xna_ * dphi[1][i]/wpc[i];
        dy_dxsi1 += yna_ * dphi[0][i]/wpc[i];
        dy_dxsi2 += yna_ * dphi[1][i]/wpc[i];   
    };

    //Computing the jacobian determinant 
    double djacpp = dx_dxsi1 * dy_dxsi2 - dx_dxsi2 * dy_dxsi1;

    //Computing jacobian inverse matrix (transposed) - dxsi/dx
    ainv_[0][0] =  dy_dxsi2 / djacpp;
    ainv_[0][1] = -dy_dxsi1 / djacpp;
    ainv_[1][0] = -dx_dxsi2 / djacpp;
    ainv_[1][1] =  dx_dxsi1 / djacpp;

    // Computing the quadrature jacobian matrix - dx/dxqsi
    int degm = (*iparametersC)[NpatchC_] -> getDegree(0);
    int degn = (*iparametersC)[NpatchC_] -> getDegree(1);
    int npcm = (*iparametersC)[NpatchC_] -> getNcp(0);
    int npcn = (*iparametersC)[NpatchC_] -> getNcp(1);   
    double *uknot_ = (*iparametersC)[NpatchC_] -> getuKnot();
    double *vknot_ = (*iparametersC)[NpatchC_]-> getvKnot();  

    int uind = inc_[0];
    int vind = inc_[1]; 

    double u1 = uknot_[uind];
    double u2 = uknot_[uind+1];
    double v1 = vknot_[vind];
    double v2 = vknot_[vind +1];

    double dxsi1_dqxsi1 = (u2-u1)*0.5;
    double dxsi1_dqxsi2 = 0.;
    double dxsi2_dqxsi1 = 0.;
    double dxsi2_dqxsi2 = (v2-v1)*0.5;

    quadJacMat [0][0] = dxsi1_dqxsi1 * dx_dxsi1 + dxsi2_dqxsi1 * dx_dxsi2;
    quadJacMat [0][1] = dxsi1_dqxsi2 * dx_dxsi1 + dxsi2_dqxsi2 * dx_dxsi2;
    quadJacMat [1][0] = dxsi1_dqxsi1 * dy_dxsi1 + dxsi2_dqxsi1 * dy_dxsi2;
    quadJacMat [1][1] = dxsi1_dqxsi2 * dy_dxsi1 + dxsi2_dqxsi2 * dy_dxsi2;

    for (int i = 0; i < 2; ++i) delete [] dphi[i];
    delete [] dphi;

    return;
    
};


template<>
void Element<2>::getQuadJacobianMatrix_ISO(double *xsi, double **quadJacMatInv){

    QuadShapeFunction<2>    shapeQuad;

    double &alpha_f = parameters->getAlphaF();
    
    //data from iga element
    int *inc_ = (*nodes_)[connect_[8]] -> getINC();
    double wpc[9] = {};
    for (int i = 0; i<9; i++) wpc[i] = (*nodes_)[connect_[i]] -> getWeightPC();

    double **dphi;
    dphi = new double*[2];
    for (int i = 0; i < 2; ++i) dphi[i] = new double[9];
    
    shapeQuad.evaluateGradientIso(xsi,dphi,wpc,inc_,(*iparameters),Npatch_);

    double dx_dxsi1 = 0.;
    double dx_dxsi2 = 0.;
    double dy_dxsi1 = 0.;
    double dy_dxsi2 = 0.;
    
    //Computing Jacobian matrix (physical to parametric)- dx/dxsi
    for (int i=0; i<9; i++){
        
        double xna_ = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(0) + 
                      (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(0);
        double yna_ = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(1) + 
                      (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(1);

        dx_dxsi1 += xna_ * dphi[0][i]/wpc[i];
        dx_dxsi2 += xna_ * dphi[1][i]/wpc[i];
        dy_dxsi1 += yna_ * dphi[0][i]/wpc[i];
        dy_dxsi2 += yna_ * dphi[1][i]/wpc[i];   
    };


    // Computing the Jacobian matrix (physical to quadrature) - dx/dxqsi
    int degm = (*iparameters)[Npatch_] -> getDegree(0);
    int degn = (*iparameters)[Npatch_] -> getDegree(1);
    int npcm = (*iparameters)[Npatch_] -> getNcp(0);
    int npcn = (*iparameters)[Npatch_] -> getNcp(1);   

    double *uknot_ = (*iparameters)[Npatch_] -> getuKnot();
    double *vknot_ = (*iparameters)[Npatch_] -> getvKnot();  

    double u1 = uknot_[inc_[0]];
    double u2 = uknot_[inc_[0]+1];
    double v1 = vknot_[inc_[1]];
    double v2 = vknot_[inc_[1] +1];

    double dxsi1_dqxsi1 = (u2-u1)*0.5;
    double dxsi1_dqxsi2 = 0.;
    double dxsi2_dqxsi1 = 0.;
    double dxsi2_dqxsi2 = (v2-v1)*0.5;

    double quadJacMat[2][2];

    quadJacMat [0][0] = dxsi1_dqxsi1 * dx_dxsi1 + dxsi2_dqxsi1 * dx_dxsi2;
    quadJacMat [0][1] = dxsi1_dqxsi2 * dx_dxsi1 + dxsi2_dqxsi2 * dx_dxsi2;
    quadJacMat [1][0] = dxsi1_dqxsi1 * dy_dxsi1 + dxsi2_dqxsi1 * dy_dxsi2;
    quadJacMat [1][1] = dxsi1_dqxsi2 * dy_dxsi1 + dxsi2_dqxsi2 * dy_dxsi2;

    djac_ = quadJacMat[0][0] * quadJacMat [1][1] - quadJacMat [0][1] * quadJacMat [1][0];

    //Computing inverse Jacobian (quadrature to physical) - dxqsi/dx
    quadJacMatInv[0][0] =  quadJacMat [1][1]/djac_;
    quadJacMatInv[0][1] = -quadJacMat [0][1]/djac_;
    quadJacMatInv[1][0] = -quadJacMat [1][0]/djac_;
    quadJacMatInv[1][1] =  quadJacMat [0][0]/djac_;

    for (int i = 0; i < 2; ++i) delete [] dphi[i];
    delete [] dphi;
}


//------------------------------------------------------------------------------
//-----------------------------SPATIAL DERIVATIVES------------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getSpatialDerivatives_FEM(double *xsi, double **ainv_,double **dphi_dx) {
    
    QuadShapeFunction<2>    shapeQuad;

    double **dphi;
    dphi = new double*[2];
    for (int i = 0; i < 2; ++i) dphi[i] = new double[6];
    
    shapeQuad.evaluateGradientFem(xsi,dphi);
    
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 6; ++j)
            dphi_dx[i][j] = 0.;

    //Quadratic shape functions spatial first derivatives
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 6; k++)
                dphi_dx[i][k] += ainv_[i][j]*dphi[j][k];

    for (int i = 0; i < 2; ++i) delete [] dphi[i];
    delete [] dphi;

    return;
};


template<>
void Element<2>::getSpatialDerivatives_ISO(double *xsi, double **ainv_, double **dphi_dx) {
    

    QuadShapeFunction<2>    shapeQuad;
    
    double **dphi;
    dphi = new double*[2];
    for (int i = 0; i < 2; ++i) dphi[i] = new double[9];

    double wpc[9] = {};
    for (int i = 0; i<9; i++) wpc[i] = (*nodes_)[connect_[i]] -> getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]] -> getINC();
    
    shapeQuad.evaluateGradientIso(xsi,dphi,wpc,inc_,(*iparameters),Npatch_);

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 9; ++j)
            dphi_dx[i][j] = 0.;

    //Quadratic shape functions spatial first derivatives
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 9; k++)
                dphi_dx[i][k] += ainv_[i][j]*dphi[j][k];

    for (int i = 0; i < 2; ++i) delete [] dphi[i];
    delete [] dphi;
    return;
};


template<>
void Element<2>::getSpatialDerivatives_COARSE_ISO(double *xsi, double **ainv_, double **dphi_dx) {
    

    QuadShapeFunction<2>    shapeQuad;
    
    double **dphi;
    dphi = new double*[2];
    for (int i = 0; i < 2; ++i) dphi[i] = new double[9];

    double wpc[9] = {};
    for (int i = 0; i<9; i++) wpc[i] = (*nodesC_)[connectC_[i]] -> getWeightPC();
    int *inc_ = (*nodesC_)[connectC_[8]] -> getINC();
    
    shapeQuad.evaluateGradientIso(xsi,dphi,wpc,inc_,(*iparametersC),NpatchC_);

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 9; ++j)
            dphi_dx[i][j] = 0.;

    //Quadratic shape functions spatial first derivatives
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 9; k++)
                dphi_dx[i][k] += ainv_[i][j]*dphi[j][k];

    for (int i = 0; i < 2; ++i) delete [] dphi[i];
    delete [] dphi;
    return;
};

template<>
void Element<2>::getSecondSpatialDerivatives_FEM(double *xsi, double **ainv_, double ***ddphi_dx){

   QuadShapeFunction<2>    shapeQuad;

    double inter[2][2][6] = {};

    double ainvT_[2][2] = {};

    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            ainvT_[i][j] = ainv_[j][i];
        };
    };

    double ***ddphi;
    ddphi = new double**[2];
    for (int i = 0; i < 2; ++i) {
        ddphi[i] = new double*[2];
        for (int j = 0; j < 2; j++) ddphi[i][j] = new double[6];
    }
        
    shapeQuad.evaluateHessianFem(xsi,ddphi);
    
    //Quadratic shape functions spatial second derivatives
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                for (int nf = 0; nf < 6; nf++)
                    inter[i][j][nf] += ainvT_[i][k]*ddphi[k][j][nf];

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 6; ++k)
                ddphi_dx[i][j][k] = 0.; 
                          
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                for (int nf = 0; nf < 6; nf++)
                    ddphi_dx[i][j][nf] += inter[i][k][nf] * ainv_[k][j];


    for (int i = 0; i < 2; ++i){
        for (int j = 0; j < 2; j++){
            delete [] ddphi[i][j];
        };
        delete [] ddphi[i];
    };
    delete [] ddphi;    

    return;     
};


template<>
void Element<2>::getSecondSpatialDerivatives_ISO(double *xsi, double **ainv_, double ***ddphi_dx){

    QuadShapeFunction<2>    shapeQuad;

    double ainvT_[2][2] = {};

    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            ainvT_[i][j] = ainv_[j][i];
        };
    };

    double ***ddphi;
    ddphi = new double**[2];
    for (int i = 0; i < 2; ++i) {
        ddphi[i] = new double*[2];
        for (int j = 0; j < 2; j++) ddphi[i][j] = new double[9];
    }

    double wpc[9] = {};
    for (int i = 0; i < 9; i++) wpc[i] = (*nodes_)[connect_[i]] -> getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]] -> getINC();       
    
    shapeQuad.evaluateHessianIso(xsi,ddphi,wpc,inc_,(*iparameters),Npatch_);
    
    //Quadratic shape functions spatial second derivatives
    double inter[2][2][9] = {};
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                for (int nf = 0; nf < 9; nf++)
                    inter[i][j][nf] += ainvT_[i][k]*ddphi[k][j][nf];

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 9; ++k)
                ddphi_dx[i][j][k] = 0.;

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                for (int nf = 0; nf < 9; nf++)
                    ddphi_dx[i][j][nf] += inter[i][k][nf] * ainv_[k][j];


    for (int i = 0; i < 2; ++i){
        for (int j = 0; j < 2; j++){
            delete [] ddphi[i][j];
        };
        delete [] ddphi[i];
    };
    delete [] ddphi;    

    return;

};

template<>
void Element<2>::getSecondSpatialDerivatives_COARSE_ISO(double *xsi, double **ainv_, double ***ddphi_dx){

    QuadShapeFunction<2>    shapeQuad;

    double ainvT_[2][2] = {};

    for (int i = 0; i < 2; i++){
        for (int j = 0; j < 2; j++){
            ainvT_[i][j] = ainv_[j][i];
        };
    };

    double ***ddphi;
    ddphi = new double**[2];
    for (int i = 0; i < 2; ++i) {
        ddphi[i] = new double*[2];
        for (int j = 0; j < 2; j++) ddphi[i][j] = new double[9];
    }

    double wpc[9] = {};
    for (int i = 0; i<9; i++) wpc[i] = (*nodesC_)[connectC_[i]] -> getWeightPC();
    int *inc_ = (*nodesC_)[connectC_[8]] -> getINC();       
    
    shapeQuad.evaluateHessianIso(xsi,ddphi,wpc,inc_,(*iparametersC),NpatchC_);
    
    //Quadratic shape functions spatial second derivatives
    double inter[2][2][9] = {};
    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                for (int nf = 0; nf < 9; nf++)
                    inter[i][j][nf] += ainvT_[i][k]*ddphi[k][j][nf];

    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 9; ++k)
                ddphi_dx[i][j][k] = 0.;

    for (int i = 0; i < 2; i++)
        for (int j = 0; j < 2; j++)
            for (int k = 0; k < 2; k++)
                for (int nf = 0; nf < 9; nf++)
                    ddphi_dx[i][j][nf] += inter[i][k][nf] * ainv_[k][j];


    for (int i = 0; i < 2; ++i){
        for (int j = 0; j < 2; j++){
            delete [] ddphi[i][j];
        };
        delete [] ddphi[i];
    };
    delete [] ddphi;    

    return;

};

//------------------------------------------------------------------------------
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getVelAndDerivatives_FEM(double *phi_, double **dphi_dx) {

	x_ = 0.;          y_ = 0.;
    
    u_ = 0.;          v_ = 0.;
    
    uprev_ = 0.;      vprev_ = 0.;
    
    du_dx = 0.;       du_dy = 0.;    
    dv_dx = 0.;       dv_dy = 0.;
    
    duprev_dx = 0.;   duprev_dy = 0.;
    dvprev_dx = 0.;   dvprev_dy = 0.;
    
    ax_ = 0.;         ay_ = 0.;
    
    axprev_ = 0.;     ayprev_ = 0.;
       
    p_ = 0.;
    
    dp_dx = 0.;       dp_dy = 0.;

    umesh_ = 0.;      vmesh_ = 0.;
    
    umeshprev_ = 0.;  vmeshprev_ = 0.;
    
    dumesh_dx = 0.;   dumesh_dy = 0.;    
    dvmesh_dx = 0.;   dvmesh_dy = 0.;
    
    dumeshprev_dx = 0.;dumeshprev_dy = 0.;
    dvmeshprev_dx = 0.;dvmeshprev_dy = 0.;
      
    double &alpha_f = parameters->getAlphaF();

    //Interpolates the velocity components and its spatial derivatives
    for (int i = 0; i < 6; i++){
        
        //coordinates
        double xna_ = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(0) + 
                      (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(0);
        double yna_ = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(1) + 
                      (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(1);
        
        x_ += xna_ * phi_[i];
        y_ += yna_ * phi_[i];

        //velocity
        u_ += (*nodes_)[connect_[i]] -> getVelocity(0) * phi_[i];
        v_ += (*nodes_)[connect_[i]] -> getVelocity(1) * phi_[i];

        //previous velocity
        uprev_ += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * phi_[i];
        vprev_ += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * phi_[i];

        //velocity derivatives
        du_dx += (*nodes_)[connect_[i]] -> getVelocity(0) * dphi_dx[0][i];
        du_dy += (*nodes_)[connect_[i]] -> getVelocity(0) * dphi_dx[1][i];
        dv_dx += (*nodes_)[connect_[i]] -> getVelocity(1) * dphi_dx[0][i];
        dv_dy += (*nodes_)[connect_[i]] -> getVelocity(1) * dphi_dx[1][i];

        //previous velocity derivatives
        duprev_dx += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * dphi_dx[0][i];
        duprev_dy += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * dphi_dx[1][i];
        dvprev_dx += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * dphi_dx[0][i];
        dvprev_dy += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * dphi_dx[1][i];

        //acceleration
        ax_ += (*nodes_)[connect_[i]] -> getAcceleration(0) * phi_[i];
        ay_ += (*nodes_)[connect_[i]] -> getAcceleration(1) * phi_[i];
        
        //previous acceleration
        axprev_ += (*nodes_)[connect_[i]] -> getPreviousAcceleration(0) * phi_[i];
        ayprev_ += (*nodes_)[connect_[i]] -> getPreviousAcceleration(1) * phi_[i];

        //pressure
        p_ += (*nodes_)[connect_[i]] -> getPressure() * phi_[i];

        //pressure derivatives
        dp_dx += (*nodes_)[connect_[i]] -> getPressure() * dphi_dx[0][i];
        dp_dy += (*nodes_)[connect_[i]] -> getPressure() * dphi_dx[1][i];

        //mesh velocity
        umesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * phi_[i];
        vmesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * phi_[i];
        
        //previous mesh velocity
        umeshprev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * phi_[i];
        vmeshprev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * phi_[i];

        //mesh velocity derivatives
        dumesh_dx += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * dphi_dx[0][i];  
        dumesh_dy += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * dphi_dx[1][i];   
    	dvmesh_dx += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * dphi_dx[0][i];       
    	dvmesh_dy += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * dphi_dx[1][i];
    	
    	//previous mesh velocity derivatives
    	dumeshprev_dx += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * dphi_dx[0][i];   
    	dumeshprev_dy += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * dphi_dx[1][i];
    	dvmeshprev_dx += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * dphi_dx[0][i];   
    	dvmeshprev_dy += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * dphi_dx[1][i];

    };  

    return;
};

template<>
void Element<2>::getVelAndDerivatives_ISO(double *phi_, double **dphi_dx) {

	x_ = 0.;          y_ = 0.;
    
    u_ = 0.;          v_ = 0.;
    
    uprev_ = 0.;      vprev_ = 0.;
    
    du_dx = 0.;       du_dy = 0.;    
    dv_dx = 0.;       dv_dy = 0.;
    
    duprev_dx = 0.;   duprev_dy = 0.;
    dvprev_dx = 0.;   dvprev_dy = 0.;
    
    ax_ = 0.;         ay_ = 0.;
    axprev_ = 0.;     ayprev_ = 0.;
       
    p_ = 0.;
    
    dp_dx = 0.;       dp_dy = 0.;

    umesh_ = 0.;      vmesh_ = 0.;
    
    umeshprev_ = 0.;  vmeshprev_ = 0.;
    
    dumesh_dx = 0.;   dumesh_dy = 0.;    
    dvmesh_dx = 0.;   dvmesh_dy = 0.;
    
    dumeshprev_dx = 0.;dumeshprev_dy = 0.;
    dvmeshprev_dx = 0.;dvmeshprev_dy = 0.;
    
    double &alpha_f = parameters->getAlphaF();

    //Interpolates the velocity components and its spatial derivatives
    for (int i = 0; i < 9; i++){
        
    	//Coordinates
        double xna_ = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(0) + 
                      (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(0);
        double yna_ = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(1) + 
                      (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(1);
        
        double wei = (*nodes_)[connect_[i]] -> getWeightPC();
        
        x_ += xna_ * phi_[i]/wei;
        y_ += yna_ * phi_[i]/wei;

        //velocity
        u_ += (*nodes_)[connect_[i]] -> getVelocity(0) * phi_[i];
        v_ += (*nodes_)[connect_[i]] -> getVelocity(1) * phi_[i];

        //previous velocity
        uprev_ += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * phi_[i];
        vprev_ += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * phi_[i];

        //velocity derivatives
        du_dx += (*nodes_)[connect_[i]] -> getVelocity(0) * dphi_dx[0][i];
        du_dy += (*nodes_)[connect_[i]] -> getVelocity(0) * dphi_dx[1][i];
        dv_dx += (*nodes_)[connect_[i]] -> getVelocity(1) * dphi_dx[0][i];
        dv_dy += (*nodes_)[connect_[i]] -> getVelocity(1) * dphi_dx[1][i];

        //previous velocity derivatives
        duprev_dx += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * dphi_dx[0][i];
        duprev_dy += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * dphi_dx[1][i];
        dvprev_dx += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * dphi_dx[0][i];
        dvprev_dy += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * dphi_dx[1][i];

        //acceleration
        ax_ += (*nodes_)[connect_[i]] -> getAcceleration(0) * phi_[i];
        ay_ += (*nodes_)[connect_[i]] -> getAcceleration(1) * phi_[i];

        //previous acceleration
        axprev_ += (*nodes_)[connect_[i]] -> getPreviousAcceleration(0) * phi_[i];
        ayprev_ += (*nodes_)[connect_[i]] -> getPreviousAcceleration(1) * phi_[i];

        //pressure
        p_ += (*nodes_)[connect_[i]] -> getPressure() * phi_[i];

        //Pressure derivatives
        dp_dx += (*nodes_)[connect_[i]] -> getPressure() * dphi_dx[0][i];
        dp_dy += (*nodes_)[connect_[i]] -> getPressure() * dphi_dx[1][i];

        //mesh velocity
        umesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * phi_[i];
        vmesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * phi_[i];

        //previous mesh velocity
        umeshprev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * phi_[i];
        vmeshprev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * phi_[i];

        //mesh velocity derivatives
        dumesh_dx += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * dphi_dx[0][i];  
        dumesh_dy += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * dphi_dx[1][i];   
    	dvmesh_dx += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * dphi_dx[0][i];       
    	dvmesh_dy += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * dphi_dx[1][i];
    	
    	//previous mesh velocity derivatives
    	dumeshprev_dx += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * dphi_dx[0][i];   
    	dumeshprev_dy += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * dphi_dx[1][i];
    	dvmeshprev_dx += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * dphi_dx[0][i];   
    	dvmeshprev_dy += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * dphi_dx[1][i];

    };  

    return;
};

template<>
void Element<2>::getInterpolatedVariablesSameMesh_FEM(double *phi_, double **dphi_dx) {

    u_ = 0.;          v_ = 0.;
    
    uprev_ = 0.;      vprev_ = 0.;
    
    du_dx = 0.;       du_dy = 0.;    
    dv_dx = 0.;       dv_dy = 0.;
    
    duprev_dx = 0.;   duprev_dy = 0.;
    dvprev_dx = 0.;   dvprev_dy = 0.;
    
    lamx_ = 0.;       lamy_ = 0.;   
    
    lamx_dx = 0.;     lamx_dy = 0.;
    lamy_dx = 0.;     lamy_dy = 0.;     
    
    umesh_ = 0.;      vmesh_ = 0.; 
    
    umeshprev_ = 0.;  vmeshprev_ = 0.;
    
        
    //Interpolates the velocity components and its spatial derivatives
    for (int i = 0; i < 6; i++){
        
        //velocity
        u_ += (*nodes_)[connect_[i]] -> getVelocity(0) * phi_[i];
        v_ += (*nodes_)[connect_[i]] -> getVelocity(1) * phi_[i];

        //previous velocity
        uprev_ += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * phi_[i];
        vprev_ += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * phi_[i];

        //velocity derivatives
        du_dx += (*nodes_)[connect_[i]] -> getVelocity(0) * dphi_dx[0][i];
        du_dy += (*nodes_)[connect_[i]] -> getVelocity(0) * dphi_dx[1][i];
        dv_dx += (*nodes_)[connect_[i]] -> getVelocity(1) * dphi_dx[0][i];
        dv_dy += (*nodes_)[connect_[i]] -> getVelocity(1) * dphi_dx[1][i];

        //previous velocity derivatives
        duprev_dx += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * dphi_dx[0][i];
        duprev_dy += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * dphi_dx[1][i];
        dvprev_dx += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * dphi_dx[0][i];
        dvprev_dy += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * dphi_dx[1][i];

        //lagrange multipliers
        lamx_ += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * phi_[i];       
        lamy_ += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * phi_[i];  

        //lagrange multipliers derivatives
        lamx_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[0][i];     
        lamx_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[1][i];
        lamy_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[0][i];      
        lamy_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[1][i];  

        //mesh velocity
        umesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * phi_[i];
        vmesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * phi_[i];

        //previous mesh velocity
        umeshprev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * phi_[i];
        vmeshprev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * phi_[i];

    };  

    return;
};


template<>
void Element<2>::getInterpolatedVariablesSameMeshArlqStab_FEM(double *phi_, double **dphi_dx, double ***ddphi_dx) {


    u_ = 0.;          v_ = 0.;
    
    uprev_ = 0.;      vprev_ = 0.;
    
    du_dx = 0.;       du_dy = 0.;    
    dv_dx = 0.;       dv_dy = 0.;
    
    duprev_dx = 0.;   duprev_dy = 0.;
    dvprev_dx = 0.;   dvprev_dy = 0.;
    
    dax_dx = 0.;       dax_dy = 0.;
    day_dx = 0.;       day_dy = 0.;
    
    daxprev_dx = 0.;   daxprev_dy = 0.;
    dayprev_dx = 0.;   dayprev_dy = 0.;
    
    ddu_dxdx = 0.;       ddu_dxdy = 0.;    
    ddu_dydx = 0.;       ddu_dydy = 0.; 
    ddv_dxdx = 0.;       ddv_dxdy = 0.;    
    ddv_dydx = 0.;       ddv_dydy = 0.;
    
    dduprev_dxdx = 0.;   dduprev_dxdy = 0.;
    dduprev_dydx = 0.;   dduprev_dydy = 0.;
    ddvprev_dxdx = 0.;   ddvprev_dxdy = 0.;
    ddvprev_dydx = 0.;   ddvprev_dydy = 0.;

    ddp_dxdx = 0.; ddp_dxdy = 0.;
    ddp_dydx = 0.; ddp_dydy = 0.;   

    lamx_dx = 0.;     lamx_dy = 0.;
    lamy_dx = 0.;     lamy_dy = 0.;  
    
    umesh_ = 0.;      vmesh_ = 0.; 
    
    umeshprev_ = 0.;  vmeshprev_ = 0.;
    
    dumesh_dx = 0.;   dumesh_dy = 0.;    
    dvmesh_dx = 0.;   dvmesh_dy = 0.;
    
    dumeshprev_dx = 0.;dumeshprev_dy = 0.;
    dvmeshprev_dx = 0.;dvmeshprev_dy = 0.;
    
    //Interpolates the velocity components and its spatial derivatives
    for (int i = 0; i < 6; i++){
        
        //velocity
        u_ += (*nodes_)[connect_[i]] -> getVelocity(0) * phi_[i];
        v_ += (*nodes_)[connect_[i]] -> getVelocity(1) * phi_[i];

        //previous velocity
        uprev_ += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * phi_[i];
        vprev_ += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * phi_[i];

         //velocity derivatives
        du_dx += (*nodes_)[connect_[i]] -> getVelocity(0) * dphi_dx[0][i];
        du_dy += (*nodes_)[connect_[i]] -> getVelocity(0) * dphi_dx[1][i];
        dv_dx += (*nodes_)[connect_[i]] -> getVelocity(1) * dphi_dx[0][i];
        dv_dy += (*nodes_)[connect_[i]] -> getVelocity(1) * dphi_dx[1][i];

        //previous velocity derivatives
        duprev_dx += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * dphi_dx[0][i];
        duprev_dy += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * dphi_dx[1][i];
        dvprev_dx += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * dphi_dx[0][i];
        dvprev_dy += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * dphi_dx[1][i];

        //acceleration
        dax_dx += (*nodes_)[connect_[i]] -> getAcceleration(0) * dphi_dx[0][i];
        dax_dy += (*nodes_)[connect_[i]] -> getAcceleration(0) * dphi_dx[1][i];
        
        day_dx += (*nodes_)[connect_[i]] -> getAcceleration(1) * dphi_dx[0][i];
        day_dy += (*nodes_)[connect_[i]] -> getAcceleration(1) * dphi_dx[1][i];

        //previous acceleration
        daxprev_dx += (*nodes_)[connect_[i]] -> getPreviousAcceleration(0) * dphi_dx[0][i];
        daxprev_dy += (*nodes_)[connect_[i]] -> getPreviousAcceleration(0) * dphi_dx[1][i];
        
        dayprev_dx += (*nodes_)[connect_[i]] -> getPreviousAcceleration(1) * dphi_dx[0][i];
        dayprev_dy += (*nodes_)[connect_[i]] -> getPreviousAcceleration(1) * dphi_dx[1][i];

        //velocity second derivatives
        ddu_dxdx += (*nodes_)[connect_[i]] -> getVelocity(0) * ddphi_dx[0][0][i];
        ddu_dxdy += (*nodes_)[connect_[i]] -> getVelocity(0) * ddphi_dx[0][1][i];
        ddu_dxdx += (*nodes_)[connect_[i]] -> getVelocity(0) * ddphi_dx[1][0][i];
        ddu_dxdy += (*nodes_)[connect_[i]] -> getVelocity(0) * ddphi_dx[1][1][i];
        ddv_dxdx += (*nodes_)[connect_[i]] -> getVelocity(1) * ddphi_dx[0][0][i];
        ddv_dxdy += (*nodes_)[connect_[i]] -> getVelocity(1) * ddphi_dx[0][1][i];
        ddv_dxdx += (*nodes_)[connect_[i]] -> getVelocity(1) * ddphi_dx[1][0][i];
        ddv_dxdy += (*nodes_)[connect_[i]] -> getVelocity(1) * ddphi_dx[1][1][i];

        //previous second velocity derivatives
        dduprev_dxdx += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * ddphi_dx[0][0][i];
        dduprev_dxdy += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * ddphi_dx[0][1][i];
        dduprev_dydx += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * ddphi_dx[1][0][i];
        dduprev_dydy += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * ddphi_dx[1][1][i];
        ddvprev_dxdx += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * ddphi_dx[0][0][i];
        ddvprev_dxdy += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * ddphi_dx[0][1][i];
        ddvprev_dydx += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * ddphi_dx[1][0][i];
        ddvprev_dydy += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * ddphi_dx[1][1][i];

        //pressure
        ddp_dxdx += (*nodes_)[connect_[i]] -> getPressure()* ddphi_dx[0][0][i]; 
        ddp_dxdy += (*nodes_)[connect_[i]] -> getPressure()* ddphi_dx[0][1][i];
        ddp_dydx += (*nodes_)[connect_[i]] -> getPressure()* ddphi_dx[1][0][i]; 
        ddp_dydy += (*nodes_)[connect_[i]] -> getPressure()* ddphi_dx[1][1][i]; 
        
        //lagrange multipliers derivatives
        lamx_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[0][i];     
        lamx_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[1][i];
        lamy_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[0][i];      
        lamy_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[1][i];  
        

        //mesh velocity
        umesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * phi_[i];
        vmesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * phi_[i];

        //previous mesh velocity
        umeshprev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * phi_[i];
        vmeshprev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * phi_[i];

        //mesh velocity derivatives
        dumesh_dx += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * dphi_dx[0][i];  
        dumesh_dy += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * dphi_dx[1][i];   
        dvmesh_dx += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * dphi_dx[0][i];       
        dvmesh_dy += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * dphi_dx[1][i];
        
        //previous mesh velocity derivatives
        dumeshprev_dx += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * dphi_dx[0][i];   
        dumeshprev_dy += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * dphi_dx[1][i];
        dvmeshprev_dx += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * dphi_dx[0][i];   
        dvmeshprev_dy += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * dphi_dx[1][i];


    };  

    return;
};

template<>
void Element<2>::getInterpolatedVariablesSameMesh_ISO(double *phi_, double **dphi_dx) {

    u_ = 0.;          v_ = 0.;
    
    uprev_ = 0.;      vprev_ = 0.;
    
    du_dx = 0.;       du_dy = 0.;    
    dv_dx = 0.;       dv_dy = 0.;
    
    duprev_dx = 0.;   duprev_dy = 0.;
    dvprev_dx = 0.;   dvprev_dy = 0.;
    
    lamx_ = 0.;       lamy_ = 0.;	
   
    lamx_dx = 0.;     lamx_dy = 0.;
    lamy_dx = 0.;     lamy_dy = 0.;  
    
    umesh_ = 0.;      vmesh_ = 0.; 
    
    umeshprev_ = 0.;  vmeshprev_ = 0.;    
    
        
    //Interpolates the velocity components and its spatial derivatives
    for (int i = 0; i < 9; i++){
        
        //velocity
        u_ += (*nodes_)[connect_[i]] -> getVelocity(0) * phi_[i];
        v_ += (*nodes_)[connect_[i]] -> getVelocity(1) * phi_[i];

        //previous velocity
        uprev_ += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * phi_[i];
        vprev_ += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * phi_[i];

        //velocity derivatives
        du_dx += (*nodes_)[connect_[i]] -> getVelocity(0) * dphi_dx[0][i];
        du_dy += (*nodes_)[connect_[i]] -> getVelocity(0) * dphi_dx[1][i];
        dv_dx += (*nodes_)[connect_[i]] -> getVelocity(1) * dphi_dx[0][i];
        dv_dy += (*nodes_)[connect_[i]] -> getVelocity(1) * dphi_dx[1][i];

        //previous velocity derivatives
        duprev_dx += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * dphi_dx[0][i];
        duprev_dy += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * dphi_dx[1][i];
        dvprev_dx += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * dphi_dx[0][i];
        dvprev_dy += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * dphi_dx[1][i];

        //lagrange multipliers
        lamx_ += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * phi_[i];       
        lamy_ += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * phi_[i];	

    	//lagrange multipliers derivatives
    	lamx_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[0][i];     
    	lamx_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[1][i];
    	lamy_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[0][i];      
    	lamy_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[1][i];   

    	//mesh velocity
        umesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * phi_[i];
        vmesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * phi_[i];

        //previous mesh velocity
        umeshprev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * phi_[i];
        vmeshprev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * phi_[i];

    };

    return;
};

template<>
void Element<2>::getInterpolatedVariablesSameMeshArlqStab_ISO(double *phi_, double **dphi_dx, double ***ddphi_dx) {


    u_ = 0.;          v_ = 0.;
    
    uprev_ = 0.;      vprev_ = 0.;

    du_dx = 0.;       du_dy = 0.;    
    dv_dx = 0.;       dv_dy = 0.;

    duprev_dx = 0.;   duprev_dy = 0.;
    dvprev_dx = 0.;   dvprev_dy = 0.;

    dax_dx = 0.;       dax_dy = 0.;
    day_dx = 0.;       day_dy = 0.;
    
    daxprev_dx = 0.;   daxprev_dy = 0.;
    dayprev_dx = 0.;   dayprev_dy = 0.;
    
    ddu_dxdx = 0.;       ddu_dxdy = 0.;    
    ddu_dydx = 0.;       ddu_dydy = 0.; 
    ddv_dxdx = 0.;       ddv_dxdy = 0.;    
    ddv_dydx = 0.;       ddv_dydy = 0.;
 
    dduprev_dxdx = 0.;   dduprev_dxdy = 0.;
    dduprev_dydx = 0.;   dduprev_dydy = 0.;
    ddvprev_dxdx = 0.;   ddvprev_dxdy = 0.;
    ddvprev_dydx = 0.;   ddvprev_dydy = 0.;

    ddp_dxdx = 0.; ddp_dxdy = 0.;
    ddp_dydx = 0.; ddp_dydy = 0.;     
 
    lamx_dx = 0.;     lamx_dy = 0.;
    lamy_dx = 0.;     lamy_dy = 0.; 

    umesh_ = 0.;      vmesh_ = 0.; 
    
    umeshprev_ = 0.;  vmeshprev_ = 0.;

    dumesh_dx = 0.;   dumesh_dy = 0.;    
    dvmesh_dx = 0.;   dvmesh_dy = 0.;
    
    dumeshprev_dx = 0.;dumeshprev_dy = 0.;
    dvmeshprev_dx = 0.;dvmeshprev_dy = 0.;

    //Interpolates the velocity components and its spatial derivatives
    for (int i = 0; i < 9; i++){
        
        //velocity
        u_ += (*nodes_)[connect_[i]] -> getVelocity(0) * phi_[i];
        v_ += (*nodes_)[connect_[i]] -> getVelocity(1) * phi_[i];

        //previous velocity
        uprev_ += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * phi_[i];
        vprev_ += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * phi_[i];

         //velocity derivatives
        du_dx += (*nodes_)[connect_[i]] -> getVelocity(0) * dphi_dx[0][i];
        du_dy += (*nodes_)[connect_[i]] -> getVelocity(0) * dphi_dx[1][i];
        
        dv_dx += (*nodes_)[connect_[i]] -> getVelocity(1) * dphi_dx[0][i];
        dv_dy += (*nodes_)[connect_[i]] -> getVelocity(1) * dphi_dx[1][i];

        //previous velocity derivatives
        duprev_dx += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * dphi_dx[0][i];
        duprev_dy += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * dphi_dx[1][i];
        
        dvprev_dx += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * dphi_dx[0][i];
        dvprev_dy += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * dphi_dx[1][i];

        //acceleration
        dax_dx += (*nodes_)[connect_[i]] -> getAcceleration(0) * dphi_dx[0][i];
        dax_dy += (*nodes_)[connect_[i]] -> getAcceleration(0) * dphi_dx[1][i];
        
        day_dx += (*nodes_)[connect_[i]] -> getAcceleration(1) * dphi_dx[0][i];
        day_dy += (*nodes_)[connect_[i]] -> getAcceleration(1) * dphi_dx[1][i];

        //previous acceleration
        daxprev_dx += (*nodes_)[connect_[i]] -> getPreviousAcceleration(0) * dphi_dx[0][i];
        daxprev_dy += (*nodes_)[connect_[i]] -> getPreviousAcceleration(0) * dphi_dx[1][i];
        
        dayprev_dx += (*nodes_)[connect_[i]] -> getPreviousAcceleration(1) * dphi_dx[0][i];
        dayprev_dy += (*nodes_)[connect_[i]] -> getPreviousAcceleration(1) * dphi_dx[1][i];

        //velocity second derivatives
        ddu_dxdx += (*nodes_)[connect_[i]] -> getVelocity(0) * ddphi_dx[0][0][i];
        ddu_dxdy += (*nodes_)[connect_[i]] -> getVelocity(0) * ddphi_dx[0][1][i];
        ddu_dxdx += (*nodes_)[connect_[i]] -> getVelocity(0) * ddphi_dx[1][0][i];
        ddu_dxdy += (*nodes_)[connect_[i]] -> getVelocity(0) * ddphi_dx[1][1][i];

        ddv_dxdx += (*nodes_)[connect_[i]] -> getVelocity(1) * ddphi_dx[0][0][i];
        ddv_dxdy += (*nodes_)[connect_[i]] -> getVelocity(1) * ddphi_dx[0][1][i];
        ddv_dxdx += (*nodes_)[connect_[i]] -> getVelocity(1) * ddphi_dx[1][0][i];
        ddv_dxdy += (*nodes_)[connect_[i]] -> getVelocity(1) * ddphi_dx[1][1][i];

        //previous second velocity derivatives
        dduprev_dxdx += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * ddphi_dx[0][0][i];
        dduprev_dxdy += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * ddphi_dx[0][1][i];
        dduprev_dydx += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * ddphi_dx[1][0][i];
        dduprev_dydy += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * ddphi_dx[1][1][i];
        ddvprev_dxdx += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * ddphi_dx[0][0][i];
        ddvprev_dxdy += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * ddphi_dx[0][1][i];
        ddvprev_dydx += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * ddphi_dx[1][0][i];
        ddvprev_dydy += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * ddphi_dx[1][1][i];

                //pressure
        ddp_dxdx += (*nodes_)[connect_[i]] -> getPressure()* ddphi_dx[0][0][i]; 
        ddp_dxdy += (*nodes_)[connect_[i]] -> getPressure()* ddphi_dx[0][1][i];
        ddp_dydx += (*nodes_)[connect_[i]] -> getPressure()* ddphi_dx[1][0][i]; 
        ddp_dydy += (*nodes_)[connect_[i]] -> getPressure()* ddphi_dx[1][1][i]; 
        
        //lagrange multipliers derivatives
        lamx_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[0][i];     
        lamx_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[1][i];
        lamy_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[0][i];      
        lamy_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[1][i];   

        //mesh velocity
        umesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * phi_[i];
        vmesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * phi_[i];

        //previous mesh velocity
        umeshprev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * phi_[i];
        vmeshprev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * phi_[i];

        //mesh velocity derivatives
        dumesh_dx += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * dphi_dx[0][i];  
        dumesh_dy += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * dphi_dx[1][i];   
    	dvmesh_dx += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * dphi_dx[0][i];       
    	dvmesh_dy += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * dphi_dx[1][i];
    	
    	//previous mesh velocity derivatives
    	dumeshprev_dx += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * dphi_dx[0][i];   
    	dumeshprev_dy += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * dphi_dx[1][i];
    	dvmeshprev_dx += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * dphi_dx[0][i];   
    	dvmeshprev_dy += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * dphi_dx[1][i];
       

    };  


    return;
};


template<>
void Element<2>::getInterpolateLagMultiplierDerivatives(double **dphi_dx){

    lamx_dx = 0.;     lamx_dy = 0.;
    lamy_dx = 0.;     lamy_dy = 0.;  

    for (int i = 0; i < 9; i++){
        
        //lagrange multipliers derivatives
        lamx_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[0][i];     
        lamx_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[1][i];
        lamy_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[0][i];      
        lamy_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[1][i];   

    };  

    return;
}

template<>
void Element<2>::getInterpolatedVariablesDifferentMesh_FEM_FEM(double *phi_, double **dphi_dx,double *phiC_, double **dphiC_dx) {

    u_ = 0.;          v_ = 0.;

    uprev_ = 0.;      vprev_ = 0.;

    du_dx = 0.;       du_dy = 0.;    
    dv_dx = 0.;       dv_dy = 0.;
    
    duprev_dx = 0.;   duprev_dy = 0.;
    dvprev_dx = 0.;   dvprev_dy = 0.;
    
    umesh_ = 0.;      vmesh_ = 0.; 
    
    umeshprev_ = 0.;  vmeshprev_ = 0.;    

    lamx_ = 0.;       lamy_ = 0.;   
    
    lamx_dx = 0.;     lamx_dy = 0.;
    lamy_dx = 0.;     lamy_dy = 0.; 
        
    //Interpolates the velocity components and its spatial derivatives
    for (int i = 0; i < 6; i++){

        //velocity
        u_ += (*nodesC_)[connectC_[i]] -> getVelocity(0) * phiC_[i];
        v_ += (*nodesC_)[connectC_[i]] -> getVelocity(1) * phiC_[i];

        //previous velocity
        uprev_ += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * phiC_[i];
        vprev_ += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * phiC_[i];

        //velocity derivatives
        du_dx += (*nodesC_)[connectC_[i]] -> getVelocity(0) * dphiC_dx[0][i];
        du_dy += (*nodesC_)[connectC_[i]] -> getVelocity(0) * dphiC_dx[1][i];
        dv_dx += (*nodesC_)[connectC_[i]] -> getVelocity(1) * dphiC_dx[0][i];
        dv_dy += (*nodesC_)[connectC_[i]] -> getVelocity(1) * dphiC_dx[1][i];

        //previous velocity derivatives
        duprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * dphiC_dx[0][i];
        duprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * dphiC_dx[1][i];
        dvprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * dphiC_dx[0][i];
        dvprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * dphiC_dx[1][i];

        //mesh velocity
        umesh_ += (*nodesC_)[connectC_[i]] -> getMeshVelocity(0) * phiC_[i];
        vmesh_ += (*nodesC_)[connectC_[i]] -> getMeshVelocity(1) * phiC_[i];

        //previous mesh velocity
        umeshprev_ += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(0) * phiC_[i];
        vmeshprev_ += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(1) * phiC_[i];

        //lagrange multipliers
        lamx_ += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * phi_[i];       
        lamy_ += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * phi_[i];  

        //lagrange multipliers derivatives
        lamx_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[0][i];     
        lamx_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[1][i];
        lamy_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[0][i];      
        lamy_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[1][i]; 

    };  


    return;
};

template<>
void Element<2>::getInterpolatedVariablesDifferentMesh_FEM_ISO(double *phi_, double **dphi_dx,double *phiC_, double **dphiC_dx) {

    u_ = 0.;          v_ = 0.;

    uprev_ = 0.;      vprev_ = 0.;

    du_dx = 0.;       du_dy = 0.;    
    dv_dx = 0.;       dv_dy = 0.;
    
    duprev_dx = 0.;   duprev_dy = 0.;
    dvprev_dx = 0.;   dvprev_dy = 0.;
  
    umesh_ = 0.;      vmesh_ = 0.; 
    
    umeshprev_ = 0.;  vmeshprev_ = 0.;    

    lamx_ = 0.;       lamy_ = 0.;   
    
    lamx_dx = 0.;     lamx_dy = 0.;
    lamy_dx = 0.;     lamy_dy = 0.;      
        
    //Interpolates the velocity components and its spatial derivatives
    for (int i = 0; i < 9; i++){

        //velocity
        u_ += (*nodesC_)[connectC_[i]] -> getVelocity(0) * phiC_[i];
        v_ += (*nodesC_)[connectC_[i]] -> getVelocity(1) * phiC_[i];

        //previous velocity
        uprev_ += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * phiC_[i];
        vprev_ += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * phiC_[i];

        //velocity derivatives
        du_dx += (*nodesC_)[connectC_[i]] -> getVelocity(0) * dphiC_dx[0][i];
        du_dy += (*nodesC_)[connectC_[i]] -> getVelocity(0) * dphiC_dx[1][i];
        dv_dx += (*nodesC_)[connectC_[i]] -> getVelocity(1) * dphiC_dx[0][i];
        dv_dy += (*nodesC_)[connectC_[i]] -> getVelocity(1) * dphiC_dx[1][i];

        //previous velocity derivatives
        duprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * dphiC_dx[0][i];
        duprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * dphiC_dx[1][i];
        dvprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * dphiC_dx[0][i];
        dvprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * dphiC_dx[1][i];

        //mesh velocity
        umesh_ += (*nodesC_)[connectC_[i]] -> getMeshVelocity(0) * phiC_[i];
        vmesh_ += (*nodesC_)[connectC_[i]] -> getMeshVelocity(1) * phiC_[i];

        //previous mesh velocity
        umeshprev_ += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(0) * phiC_[i];
        vmeshprev_ += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(1) * phiC_[i];

    };  

    for (int i = 0; i < 6; i++){

        //lagrange multipliers
        lamx_ += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * phi_[i];       
        lamy_ += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * phi_[i];  

        //lagrange multipliers derivatives
        lamx_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[0][i];     
        lamx_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[1][i];
        lamy_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[0][i];      
        lamy_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[1][i]; 
    }


    return;
};

template<>
void Element<2>::getInterpolatedVariablesDifferentMesh_ISO_FEM(double *phi_, double **dphi_dx,double *phiC_, double **dphiC_dx) {

    u_ = 0.;          v_ = 0.;

    uprev_ = 0.;      vprev_ = 0.;

    du_dx = 0.;       du_dy = 0.;    
    dv_dx = 0.;       dv_dy = 0.;
    
    duprev_dx = 0.;   duprev_dy = 0.;
    dvprev_dx = 0.;   dvprev_dy = 0.;
    
    umesh_ = 0.;      vmesh_ = 0.; 
    
    umeshprev_ = 0.;  vmeshprev_ = 0.;    
    
    lamx_ = 0.;       lamy_ = 0.;   
    
    lamx_dx = 0.;     lamx_dy = 0.;
    lamy_dx = 0.;     lamy_dy = 0.;     
        
    //Interpolates the velocity components and its spatial derivatives
    for (int i = 0; i < 6; i++){

        //velocity
        u_ += (*nodesC_)[connectC_[i]] -> getVelocity(0) * phiC_[i];
        v_ += (*nodesC_)[connectC_[i]] -> getVelocity(1) * phiC_[i];

        //previous velocity
        uprev_ += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * phiC_[i];
        vprev_ += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * phiC_[i];

        //velocity derivatives
        du_dx += (*nodesC_)[connectC_[i]] -> getVelocity(0) * dphiC_dx[0][i];
        du_dy += (*nodesC_)[connectC_[i]] -> getVelocity(0) * dphiC_dx[1][i];
        dv_dx += (*nodesC_)[connectC_[i]] -> getVelocity(1) * dphiC_dx[0][i];
        dv_dy += (*nodesC_)[connectC_[i]] -> getVelocity(1) * dphiC_dx[1][i];

        //previous velocity derivatives
        duprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * dphiC_dx[0][i];
        duprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * dphiC_dx[1][i];
        dvprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * dphiC_dx[0][i];
        dvprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * dphiC_dx[1][i];

        //mesh velocity
        umesh_ += (*nodesC_)[connectC_[i]] -> getMeshVelocity(0) * phiC_[i];
        vmesh_ += (*nodesC_)[connectC_[i]] -> getMeshVelocity(1) * phiC_[i];

        //previous mesh velocity
        umeshprev_ += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(0) * phiC_[i];
        vmeshprev_ += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(1) * phiC_[i];

    };  

    for (int i = 0; i < 9; i++){

        //lagrange multipliers
        lamx_ += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * phi_[i];       
        lamy_ += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * phi_[i];  

        //lagrange multipliers derivatives
        lamx_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[0][i];     
        lamx_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[1][i];
        lamy_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[0][i];      
        lamy_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[1][i]; 
    };


    return;



};

template<>
void Element<2>::getInterpolatedVariablesDifferentMeshArlqStab_FEM_ISO(double **dphi_dx,double *phiC_, double **dphiC_dx,
                                                                       double ***ddphiC_dx) {
    u_ = 0.;          v_ = 0.;
    
    uprev_ = 0.;      vprev_ = 0.;

    du_dx = 0.;       du_dy = 0.;    
    dv_dx = 0.;       dv_dy = 0.;

    duprev_dx = 0.;   duprev_dy = 0.;
    dvprev_dx = 0.;   dvprev_dy = 0.;

    dax_dx = 0.;       dax_dy = 0.;
    day_dx = 0.;       day_dy = 0.;
    
    daxprev_dx = 0.;   daxprev_dy = 0.;
    dayprev_dx = 0.;   dayprev_dy = 0.;
    
    ddu_dxdx = 0.;       ddu_dxdy = 0.;    
    ddu_dydx = 0.;       ddu_dydy = 0.; 
    ddv_dxdx = 0.;       ddv_dxdy = 0.;    
    ddv_dydx = 0.;       ddv_dydy = 0.;
 
    dduprev_dxdx = 0.;   dduprev_dxdy = 0.;
    dduprev_dydx = 0.;   dduprev_dydy = 0.;
    ddvprev_dxdx = 0.;   ddvprev_dxdy = 0.;
    ddvprev_dydx = 0.;   ddvprev_dydy = 0.;

    ddp_dxdx = 0.; ddp_dxdy = 0.;
    ddp_dydx = 0.; ddp_dydy = 0.;

    umesh_ = 0.;      vmesh_ = 0.; 
    
    umeshprev_ = 0.;  vmeshprev_ = 0.;

    dumesh_dx = 0.;   dumesh_dy = 0.;    
    dvmesh_dx = 0.;   dvmesh_dy = 0.;
    
    dumeshprev_dx = 0.;dumeshprev_dy = 0.;
    dvmeshprev_dx = 0.;dvmeshprev_dy = 0.;

    lamx_dx = 0.;     lamx_dy = 0.;
    lamy_dx = 0.;     lamy_dy = 0.; 


    //Interpolates the velocity components and its spatial derivatives
    for (int i = 0; i < 9; i++){
        
        //velocity
        u_ += (*nodesC_)[connectC_[i]] -> getVelocity(0) * phiC_[i];
        v_ += (*nodesC_)[connectC_[i]] -> getVelocity(1) * phiC_[i];

        //previous velocity
        uprev_ += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * phiC_[i];
        vprev_ += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * phiC_[i];

         //velocity derivatives
        du_dx += (*nodesC_)[connectC_[i]] -> getVelocity(0) * dphiC_dx[0][i];
        du_dy += (*nodesC_)[connectC_[i]] -> getVelocity(0) * dphiC_dx[1][i];
        dv_dx += (*nodesC_)[connectC_[i]] -> getVelocity(1) * dphiC_dx[0][i];
        dv_dy += (*nodesC_)[connectC_[i]] -> getVelocity(1) * dphiC_dx[1][i];

        //previous velocity derivatives
        duprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * dphiC_dx[0][i];
        duprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * dphiC_dx[1][i];
        dvprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * dphiC_dx[0][i];
        dvprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * dphiC_dx[1][i];

        //acceleration
        dax_dx += (*nodesC_)[connectC_[i]] -> getAcceleration(0) * dphiC_dx[0][i];
        dax_dy += (*nodesC_)[connectC_[i]] -> getAcceleration(0) * dphiC_dx[1][i];
        day_dx += (*nodesC_)[connectC_[i]] -> getAcceleration(1) * dphiC_dx[0][i];
        day_dy += (*nodesC_)[connectC_[i]] -> getAcceleration(1) * dphiC_dx[1][i];

        //previous acceleration
        daxprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousAcceleration(0) * dphiC_dx[0][i];
        daxprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousAcceleration(0) * dphiC_dx[1][i];
        dayprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousAcceleration(1) * dphiC_dx[0][i];
        dayprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousAcceleration(1) * dphiC_dx[1][i];

        //velocity second derivatives
        ddu_dxdx += (*nodesC_)[connectC_[i]] -> getVelocity(0) * ddphiC_dx[0][0][i];
        ddu_dxdy += (*nodesC_)[connectC_[i]] -> getVelocity(0) * ddphiC_dx[0][1][i];
        ddu_dxdx += (*nodesC_)[connectC_[i]] -> getVelocity(0) * ddphiC_dx[1][0][i];
        ddu_dxdy += (*nodesC_)[connectC_[i]] -> getVelocity(0) * ddphiC_dx[1][1][i];

        ddv_dxdx += (*nodesC_)[connectC_[i]] -> getVelocity(1) * ddphiC_dx[0][0][i];
        ddv_dxdy += (*nodesC_)[connectC_[i]] -> getVelocity(1) * ddphiC_dx[0][1][i];
        ddv_dxdx += (*nodesC_)[connectC_[i]] -> getVelocity(1) * ddphiC_dx[1][0][i];
        ddv_dxdy += (*nodesC_)[connectC_[i]] -> getVelocity(1) * ddphiC_dx[1][1][i];

        //previous second velocity derivatives
        dduprev_dxdx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * ddphiC_dx[0][0][i];
        dduprev_dxdy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * ddphiC_dx[0][1][i];
        dduprev_dydx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * ddphiC_dx[1][0][i];
        dduprev_dydy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * ddphiC_dx[1][1][i];
        ddvprev_dxdx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * ddphiC_dx[0][0][i];
        ddvprev_dxdy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * ddphiC_dx[0][1][i];
        ddvprev_dydx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * ddphiC_dx[1][0][i];
        ddvprev_dydy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * ddphiC_dx[1][1][i];
        
        //pressure
        ddp_dxdx += (*nodesC_)[connectC_[i]] -> getPressure()* ddphiC_dx[0][0][i]; 
        ddp_dxdy += (*nodesC_)[connectC_[i]] -> getPressure()* ddphiC_dx[0][1][i];
        ddp_dydx += (*nodesC_)[connectC_[i]] -> getPressure()* ddphiC_dx[1][0][i]; 
        ddp_dydy += (*nodesC_)[connectC_[i]] -> getPressure()* ddphiC_dx[1][1][i]; 


        //mesh velocity
        umesh_ += (*nodesC_)[connectC_[i]] -> getMeshVelocity(0) * phiC_[i];
        vmesh_ += (*nodesC_)[connectC_[i]] -> getMeshVelocity(1) * phiC_[i];

        //previous mesh velocity
        umeshprev_ += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(0) * phiC_[i];
        vmeshprev_ += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(1) * phiC_[i];

        //mesh velocity derivatives
        dumesh_dx += (*nodesC_)[connectC_[i]] -> getMeshVelocity(0) * dphiC_dx[0][i];  
        dumesh_dy += (*nodesC_)[connectC_[i]] -> getMeshVelocity(0) * dphiC_dx[1][i];   
    	dvmesh_dx += (*nodesC_)[connectC_[i]] -> getMeshVelocity(1) * dphiC_dx[0][i];       
    	dvmesh_dy += (*nodesC_)[connectC_[i]] -> getMeshVelocity(1) * dphiC_dx[1][i];
    	
    	//previous mesh velocity derivatives
    	dumeshprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(0) * dphiC_dx[0][i];   
    	dumeshprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(0) * dphiC_dx[1][i];
    	dvmeshprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(1) * dphiC_dx[0][i];   
    	dvmeshprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(1) * dphiC_dx[1][i];
       
    };  

    for (int i = 0; i < 6; i++){

        //lagrange multipliers derivatives
        lamx_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[0][i];     
        lamx_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[1][i];
        lamy_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[0][i];      
        lamy_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[1][i]; 
    };

    return;
};

template<>
void Element<2>::getInterpolatedVariablesDifferentMeshArlqStab_ISO_FEM(double **dphi_dx,double *phiC_, double **dphiC_dx,
                                                                       double ***ddphiC_dx) {
    u_ = 0.;          v_ = 0.;
    
    uprev_ = 0.;      vprev_ = 0.;

    du_dx = 0.;       du_dy = 0.;    
    dv_dx = 0.;       dv_dy = 0.;

    duprev_dx = 0.;   duprev_dy = 0.;
    dvprev_dx = 0.;   dvprev_dy = 0.;

    dax_dx = 0.;       dax_dy = 0.;
    day_dx = 0.;       day_dy = 0.;
    daxprev_dx = 0.;   daxprev_dy = 0.;
    dayprev_dx = 0.;   dayprev_dy = 0.;
    
    ddu_dxdx = 0.;       ddu_dxdy = 0.;    
    ddu_dydx = 0.;       ddu_dydy = 0.; 
    ddv_dxdx = 0.;       ddv_dxdy = 0.;    
    ddv_dydx = 0.;       ddv_dydy = 0.;
 
    dduprev_dxdx = 0.;   dduprev_dxdy = 0.;
    dduprev_dydx = 0.;   dduprev_dydy = 0.;
    ddvprev_dxdx = 0.;   ddvprev_dxdy = 0.;
    ddvprev_dydx = 0.;   ddvprev_dydy = 0.;

    ddp_dxdx = 0.; ddp_dxdy = 0.;
    ddp_dydx = 0.; ddp_dydy = 0.;

    umesh_ = 0.;      vmesh_ = 0.; 

    umeshprev_ = 0.;  vmeshprev_ = 0.;

    dumesh_dx = 0.;   dumesh_dy = 0.;    
    dvmesh_dx = 0.;   dvmesh_dy = 0.;
    
    dumeshprev_dx = 0.;dumeshprev_dy = 0.;
    dvmeshprev_dx = 0.;dvmeshprev_dy = 0.;

    lamx_dx = 0.;     lamx_dy = 0.;
    lamy_dx = 0.;     lamy_dy = 0.; 

    //Interpolates the velocity components and its spatial derivatives
    for (int i = 0; i < 6; i++){
        
        //velocity
        u_ += (*nodesC_)[connectC_[i]] -> getVelocity(0) * phiC_[i];
        v_ += (*nodesC_)[connectC_[i]] -> getVelocity(1) * phiC_[i];

        //previous velocity
        uprev_ += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * phiC_[i];
        vprev_ += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * phiC_[i];

        //mesh velocity
        umesh_ += (*nodesC_)[connectC_[i]] -> getMeshVelocity(0) * phiC_[i];
        vmesh_ += (*nodesC_)[connectC_[i]] -> getMeshVelocity(1) * phiC_[i];

        //previous mesh velocity
        umeshprev_ += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(0) * phiC_[i];
        vmeshprev_ += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(1) * phiC_[i];

         //velocity derivatives
        du_dx += (*nodesC_)[connectC_[i]] -> getVelocity(0) * dphiC_dx[0][i];
        du_dy += (*nodesC_)[connectC_[i]] -> getVelocity(0) * dphiC_dx[1][i];
        dv_dx += (*nodesC_)[connectC_[i]] -> getVelocity(1) * dphiC_dx[0][i];
        dv_dy += (*nodesC_)[connectC_[i]] -> getVelocity(1) * dphiC_dx[1][i];

        //previous velocity derivatives
        duprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * dphiC_dx[0][i];
        duprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * dphiC_dx[1][i];
        dvprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * dphiC_dx[0][i];
        dvprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * dphiC_dx[1][i];

        //acceleration
        dax_dx += (*nodesC_)[connectC_[i]] -> getAcceleration(0) * dphiC_dx[0][i];
        dax_dy += (*nodesC_)[connectC_[i]] -> getAcceleration(0) * dphiC_dx[1][i];
        day_dx += (*nodesC_)[connectC_[i]] -> getAcceleration(1) * dphiC_dx[0][i];
        day_dy += (*nodesC_)[connectC_[i]] -> getAcceleration(1) * dphiC_dx[1][i];

        //previous acceleration
        daxprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousAcceleration(0) * dphiC_dx[0][i];
        daxprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousAcceleration(0) * dphiC_dx[1][i];
        dayprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousAcceleration(1) * dphiC_dx[0][i];
        dayprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousAcceleration(1) * dphiC_dx[1][i];

        //velocity second derivatives
        ddu_dxdx += (*nodesC_)[connectC_[i]] -> getVelocity(0) * ddphiC_dx[0][0][i];
        ddu_dxdy += (*nodesC_)[connectC_[i]] -> getVelocity(0) * ddphiC_dx[0][1][i];
        ddu_dxdx += (*nodesC_)[connectC_[i]] -> getVelocity(0) * ddphiC_dx[1][0][i];
        ddu_dxdy += (*nodesC_)[connectC_[i]] -> getVelocity(0) * ddphiC_dx[1][1][i];

        ddv_dxdx += (*nodesC_)[connectC_[i]] -> getVelocity(1) * ddphiC_dx[0][0][i];
        ddv_dxdy += (*nodesC_)[connectC_[i]] -> getVelocity(1) * ddphiC_dx[0][1][i];
        ddv_dxdx += (*nodesC_)[connectC_[i]] -> getVelocity(1) * ddphiC_dx[1][0][i];
        ddv_dxdy += (*nodesC_)[connectC_[i]] -> getVelocity(1) * ddphiC_dx[1][1][i];

        //previous second velocity derivatives
        dduprev_dxdx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * ddphiC_dx[0][0][i];
        dduprev_dxdy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * ddphiC_dx[0][1][i];
        dduprev_dydx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * ddphiC_dx[1][0][i];
        dduprev_dydy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * ddphiC_dx[1][1][i];
        ddvprev_dxdx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * ddphiC_dx[0][0][i];
        ddvprev_dxdy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * ddphiC_dx[0][1][i];
        ddvprev_dydx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * ddphiC_dx[1][0][i];
        ddvprev_dydy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * ddphiC_dx[1][1][i];
        
        //pressure
        ddp_dxdx += (*nodesC_)[connectC_[i]] -> getPressure()* ddphiC_dx[0][0][i]; 
        ddp_dxdy += (*nodesC_)[connectC_[i]] -> getPressure()* ddphiC_dx[0][1][i];
        ddp_dydx += (*nodesC_)[connectC_[i]] -> getPressure()* ddphiC_dx[1][0][i]; 
        ddp_dydy += (*nodesC_)[connectC_[i]] -> getPressure()* ddphiC_dx[1][1][i]; 

        //mesh velocity
        umesh_ += (*nodesC_)[connectC_[i]] -> getMeshVelocity(0) * phiC_[i];
        vmesh_ += (*nodesC_)[connectC_[i]] -> getMeshVelocity(1) * phiC_[i];

        //previous mesh velocity
        umeshprev_ += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(0) * phiC_[i];
        vmeshprev_ += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(1) * phiC_[i];

        //mesh velocity derivatives
        dumesh_dx += (*nodesC_)[connectC_[i]] -> getMeshVelocity(0) * dphiC_dx[0][i];  
        dumesh_dy += (*nodesC_)[connectC_[i]] -> getMeshVelocity(0) * dphiC_dx[1][i];   
    	dvmesh_dx += (*nodesC_)[connectC_[i]] -> getMeshVelocity(1) * dphiC_dx[0][i];       
    	dvmesh_dy += (*nodesC_)[connectC_[i]] -> getMeshVelocity(1) * dphiC_dx[1][i];
    	
    	//previous mesh velocity derivatives
    	dumeshprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(0) * dphiC_dx[0][i];   
    	dumeshprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(0) * dphiC_dx[1][i];
    	dvmeshprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(1) * dphiC_dx[0][i];   
    	dvmeshprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(1) * dphiC_dx[1][i];

    };  

    for (int i = 0; i < 9; i++){

        //lagrange multipliers derivatives
        lamx_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[0][i];     
        lamx_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[1][i];
        lamy_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[0][i];      
        lamy_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[1][i]; 
    };

    return;
};

template<>
void Element<2>::getInterpolatedVariablesDifferentMesh_ISO_ISO(double *phi_, double **dphi_dx,double *phiC_, double **dphiC_dx) {

	u_ = 0.;          v_ = 0.;

	uprev_ = 0.;      vprev_ = 0.;
    
    du_dx = 0.;       du_dy = 0.;    
    dv_dx = 0.;       dv_dy = 0.;

    duprev_dx = 0.;   duprev_dy = 0.;
    dvprev_dx = 0.;   dvprev_dy = 0.;
     
    umesh_ = 0.;      vmesh_ = 0.; 
    
    umeshprev_ = 0.;  vmeshprev_ = 0.;    

    lamx_ = 0.;       lamy_ = 0.;	
    
    lamx_dx = 0.;     lamx_dy = 0.;
    lamy_dx = 0.;     lamy_dy = 0.;    
        
    //Interpolates the velocity components and its spatial derivatives
    for (int i = 0; i < 9; i++){

    	//velocity
        u_ += (*nodesC_)[connectC_[i]] -> getVelocity(0) * phiC_[i];
        v_ += (*nodesC_)[connectC_[i]] -> getVelocity(1) * phiC_[i];

        //previous velocity
        uprev_ += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * phiC_[i];
        vprev_ += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * phiC_[i];

        //velocity derivatives
        du_dx += (*nodesC_)[connectC_[i]] -> getVelocity(0) * dphiC_dx[0][i];
        du_dy += (*nodesC_)[connectC_[i]] -> getVelocity(0) * dphiC_dx[1][i];
        dv_dx += (*nodesC_)[connectC_[i]] -> getVelocity(1) * dphiC_dx[0][i];
        dv_dy += (*nodesC_)[connectC_[i]] -> getVelocity(1) * dphiC_dx[1][i];

        //previous velocity derivatives
        duprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * dphiC_dx[0][i];
        duprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(0) * dphiC_dx[1][i];
        dvprev_dx += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * dphiC_dx[0][i];
        dvprev_dy += (*nodesC_)[connectC_[i]] -> getPreviousVelocity(1) * dphiC_dx[1][i];

    	//mesh velocity
        umesh_ += (*nodesC_)[connectC_[i]] -> getMeshVelocity(0) * phiC_[i];
        vmesh_ += (*nodesC_)[connectC_[i]] -> getMeshVelocity(1) * phiC_[i];

        //previous mesh velocity
        umeshprev_ += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(0) * phiC_[i];
        vmeshprev_ += (*nodesC_)[connectC_[i]] -> getPreviousMeshVelocity(1) * phiC_[i];

    	//lagrange multipliers
        lamx_ += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * phi_[i];       
        lamy_ += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * phi_[i];	

    	//lagrange multipliers derivatives
    	lamx_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[0][i];     
    	lamx_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(0) * dphi_dx[1][i];
    	lamy_dx += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[0][i];      
    	lamy_dy += (*nodes_)[connect_[i]] -> getLagrangeMultiplier(1) * dphi_dx[1][i];  

    };  


    return;
};



//------------------------------------------------------------------------------
//-----------------------COMPUTE DRAG AND LIFT FORCES --------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::computeDragAndLiftForces_ISO(int &side, double &pDForce, double &pLForce, double &fDForce, double &fLForce, 
                                              double &dForce, double &lForce, double &aux_Mom, double &aux_Per) {

    double    &visc_ = parameters->getViscosity();
    double    localNodesBoundaryIso_[3][2];
    double    pcWeightl[3];
    double    pcWeight[9];
    int       *incL_;
   
    for (int i = 0; i<9; i++) pcWeight[i] = (*nodes_)[connect_[i]] -> getWeightPC();
    int *inc_= (*nodes_)[connect_[8]] -> getINC();     

    if(side == 0){

        pcWeightl[0] = (*nodes_)[connect_[0]] -> getWeightPC();
        pcWeightl[1] = (*nodes_)[connect_[1]] -> getWeightPC();
        pcWeightl[2] = (*nodes_)[connect_[2]] -> getWeightPC();
        
        for (int i=0; i<2; i++){
            localNodesBoundaryIso_[0][i] = ((*nodes_)[connect_[0]] -> getCoordinateValue(i))/pcWeightl[0];
            localNodesBoundaryIso_[1][i] = ((*nodes_)[connect_[1]] -> getCoordinateValue(i))/pcWeightl[1];
            localNodesBoundaryIso_[2][i] = ((*nodes_)[connect_[2]] -> getCoordinateValue(i))/pcWeightl[2];           
        };

        incL_ = (*nodes_)[connect_[2]] -> getINC(); 

    }else{

        if(side == 1){

            pcWeightl[0] = (*nodes_)[connect_[2]] -> getWeightPC();
            pcWeightl[1] = (*nodes_)[connect_[5]] -> getWeightPC();
            pcWeightl[2] = (*nodes_)[connect_[8]] -> getWeightPC();

            for (int i=0; i<2; i++){
                localNodesBoundaryIso_[0][i] = ((*nodes_)[connect_[2]] -> getCoordinateValue(i))/pcWeightl[0];
                localNodesBoundaryIso_[1][i] = ((*nodes_)[connect_[5]] -> getCoordinateValue(i))/pcWeightl[1];
                localNodesBoundaryIso_[2][i] = ((*nodes_)[connect_[8]] -> getCoordinateValue(i))/pcWeightl[2];
            };

            incL_ = (*nodes_)[connect_[8]] -> getINC(); 

        if(side == 2){

            pcWeightl[0] = (*nodes_)[connect_[8]] -> getWeightPC();
            pcWeightl[1] = (*nodes_)[connect_[7]] -> getWeightPC();
            pcWeightl[2] = (*nodes_)[connect_[6]] -> getWeightPC();
                
                for (int i=0; i<2; i++){
                localNodesBoundaryIso_[0][i] = ((*nodes_)[connect_[8]] -> getCoordinateValue(i))/pcWeightl[0];
                localNodesBoundaryIso_[1][i] = ((*nodes_)[connect_[7]] -> getCoordinateValue(i))/pcWeightl[1];
                localNodesBoundaryIso_[2][i] = ((*nodes_)[connect_[6]] -> getCoordinateValue(i))/pcWeightl[2];
            };

            incL_ = (*nodes_)[connect_[8]] -> getINC(); 
        };

        }else{

            pcWeightl[0] = (*nodes_)[connect_[6]] -> getWeightPC();
            pcWeightl[1] = (*nodes_)[connect_[3]] -> getWeightPC();
            pcWeightl[2] = (*nodes_)[connect_[0]] -> getWeightPC();
            
            for (int i=0; i<2; i++){
                localNodesBoundaryIso_[0][i] = ((*nodes_)[connect_[6]] -> getCoordinateValue(i))/pcWeightl[0];
                localNodesBoundaryIso_[1][i] = ((*nodes_)[connect_[3]] -> getCoordinateValue(i))/pcWeightl[1];
                localNodesBoundaryIso_[2][i] = ((*nodes_)[connect_[0]] -> getCoordinateValue(i))/pcWeightl[2];
            };

            incL_ = (*nodes_)[connect_[6]] -> getINC(); 
        };        
    };


    BoundaryIntegQuadratureIso<2>                      bQuad;    
    QuadShapeFunction<2>                           shapeQuad; 
    BoundShapeFunction<2>                          shapeBound;

    double n_vector[2] = {};
    double shearStress[2][2] = {};
    double load_friction[2] = {};
    double load_pressure[2] = {};

    double ident[2][2];
    ident[0][0] = 1.;
    ident[0][1] = 0.;
    ident[1][0] = 0.;
    ident[1][1] = 1.;

    double phi_[9], phiCB_[3], xsi[2], xsiB[1], weightB;

    double **dphiCB_;
    dphiCB_ = new double*[1];
    for (int i = 0; i < 1; ++i) dphiCB_[i] = new double[3];

    double **dphi_dx;
    dphi_dx = new double*[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double*[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **quadJacMat;
    quadJacMat = new double*[2];
    for (int i = 0; i < 2; ++i) quadJacMat[i] = new double[2];

    double moment = 0.;
    double per = 0.;


    int index = 0;
    for(double* it = bQuad.begin(); it != bQuad.end(); it++){
        
        xsiB[0] = bQuad.PointList(index,0);
        weightB = bQuad.WeightList(index);

        if(side == 0){
            xsi[0] = xsiB[0];
            xsi[1] = -1;
        };
        if(side == 1){
            xsi[0] = -1;
            xsi[1] = xsiB[0];
        };
        if(side == 2){
            xsi[0] = xsiB[0];
            xsi[1] = 1;
        };
        if(side == 3){
            xsi[0] = -1;
            xsi[1] = xsiB[0];
        };

        //Computes the shape functions
        shapeQuad.evaluateIso(xsi,phi_,pcWeight,inc_,(*iparameters),Npatch_);
        
        //Computes the jacobian matrix
        getJacobianMatrix_ISO(xsi,quadJacMat,ainv_);

        //Computes spatial derivatives
        getSpatialDerivatives_ISO(xsi,ainv_,dphi_dx);

        //Interpolates velocity and its derivatives values
        getVelAndDerivatives_ISO(phi_,dphi_dx);   

        //Get Functions and Derivatives on the boundary
        shapeBound.evaluateBoundaryIso(xsiB,phiCB_,pcWeightl,side,incL_,(*iparameters),Npatch_);
        shapeBound.evaluateGradientBoundaryIso(xsiB,dphiCB_,pcWeightl,side,incL_,(*iparameters),Npatch_);
       
    
        // (Jacobian parametric space to physical one)
        double Tx=0.; double Ty = 0.;
        for (int i=0; i<3; i++){
            Tx += localNodesBoundaryIso_[i][0] * dphiCB_[0][i];
            Ty += localNodesBoundaryIso_[i][1] * dphiCB_[0][i];
        };

        // (Jacobian parental space to the physical one)
        if ((side == 0) || (side == 2)){
            
            int degm = (*iparameters)[Npatch_] -> getDegree(0);
            int npcm = (*iparameters)[Npatch_] -> getNcp(0);  
            double *uknot_= (*iparameters)[Npatch_] -> getuKnot();

            int uind = incL_[0];

            double u1 = uknot_[uind];
            double u2 = uknot_[uind+1];

            double dxsi1_dqxsi1 = (u2-u1)*0.5;
            
            Tx *= dxsi1_dqxsi1;
            Ty *= dxsi1_dqxsi1;
 
        } else {

            int degn = (*iparameters)[Npatch_] -> getDegree(1);
            int npcn = (*iparameters)[Npatch_] -> getNcp(1);   
            double *vknot_= (*iparameters)[Npatch_] -> getvKnot();  

            int vind = incL_[1]; 

            double v1 = vknot_[vind];
            double v2 = vknot_[vind +1];

            double dxsi2_dqxsi2 = (v2-v1)*0.5;

            Tx *= dxsi2_dqxsi2;
            Ty *= dxsi2_dqxsi2;
        }
        
        double jacb_ = sqrt(Tx*Tx + Ty*Ty);
        
        n_vector[0] =  Ty / jacb_;
        n_vector[1] = -Tx / jacb_;

        shearStress[0][0] = 2. * visc_ * du_dx;
        shearStress[0][1] = visc_ * (du_dy + dv_dx);
        shearStress[1][0] = visc_ * (du_dy + dv_dx);
        shearStress[1][1] = 2. * visc_ * dv_dy;

        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                load_pressure[i] += -p_ * ident[i][j] * n_vector[j] * jacb_ * weightB;
                load_friction[i] += shearStress[i][j] * n_vector[j] * jacb_ * weightB;
            }
        }

        // moment +=  (-(load_pressure[0]+load_friction[0]) * (y_- 4.) 
        //            + (load_pressure[0]+load_friction[0]) * (x_ - 4.)) * jacb_ * weightB[0];
        per += jacb_ * weightB;
        
        index++;
    };

    aux_Per = per;
    aux_Mom = moment;

    pDForce = -load_pressure[0];
    pLForce = -load_pressure[1];
    
    fDForce = -load_friction[0];
    fLForce = -load_friction[1]; 

    dForce = pDForce + fDForce;
    lForce = pLForce + fLForce;


    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = 0; i < 1; ++i) delete [] dphiCB_[i];
    delete [] dphiCB_;   
    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] quadJacMat[i];
    delete [] quadJacMat;


    return;
};


template<>
void Element<2>::computeDragAndLiftForces_FEM(int &side, double &pDForce, double &pLForce, double &fDForce, double &fLForce, 
                                              double &dForce, double &lForce, double &aux_Mom, double &aux_Per) {

    double localNodesBoundary_[3][2];
    double &visc_ = parameters->getViscosity();

    if(side == 0){
        for (int i=0; i<2; i++){
            localNodesBoundary_[0][i] = (*nodes_)[connect_[1]] -> getCoordinateValue(i);
            localNodesBoundary_[1][i] = (*nodes_)[connect_[4]] -> getCoordinateValue(i);
            localNodesBoundary_[2][i] = (*nodes_)[connect_[2]] -> getCoordinateValue(i);
        };
    }else{
        if(side == 1){
            for (int i=0; i<2; i++){
                localNodesBoundary_[0][i] = (*nodes_)[connect_[2]] -> getCoordinateValue(i);
                localNodesBoundary_[1][i] = (*nodes_)[connect_[5]] -> getCoordinateValue(i);
                localNodesBoundary_[2][i] = (*nodes_)[connect_[0]] -> getCoordinateValue(i);
            };
        }else{
            for (int i=0; i<2; i++){
                localNodesBoundary_[0][i] = (*nodes_)[connect_[0]] -> getCoordinateValue(i);
                localNodesBoundary_[1][i] = (*nodes_)[connect_[3]] -> getCoordinateValue(i);
                localNodesBoundary_[2][i] = (*nodes_)[connect_[1]] -> getCoordinateValue(i);
            };
        };        
    };

    BoundaryIntegQuadrature<2>   bQuad;    
    QuadShapeFunction<2>         shapeQuad; 
    BoundShapeFunction<2>        shapeBound;


    double n_vector[2] = {};
    double shearStress[2][2] = {};
    double ident[2][2];
    double load_friction[2] = {};
    double load_pressure[2] = {};

    ident[0][0] = 1.;
    ident[0][1] = 0.;
    ident[1][0] = 0.;
    ident[1][1] = 1.;

    double phi_[6], phiB_[3], xsi[2], xsiB[1], weightB;

    double **dphiB_;
    dphiB_ = new double*[1];
    for (int i = 0; i < 1; ++i) dphiB_[i] = new double[3];

    double **dphi_dx;
    dphi_dx = new double*[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[6];

    double **ainv_;
    ainv_ = new double*[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **Jac;
    Jac = new double*[2];
    for (int i = 0; i < 2; ++i) Jac[i] = new double[2];

    double moment = 0.;
    double per = 0.;
    
    int index = 0;

    for(double* it = bQuad.begin(); it != bQuad.end(); it++){    
        
        xsiB[0] = bQuad.PointList(index,0);
        weightB= bQuad.WeightList(index);

        if(side == 2){
            xsi[0] = (-xsiB[0] + 1.) / 2.;
            xsi[1] = 0.;
        };
        if(side == 1){
            xsi[1] = (xsiB[0] + 1.) / 2.;
            xsi[0] = 0.;
        };
        if(side == 0){
            xsi[0] = (xsiB[0] + 1.) / 2.;
            xsi[1] = 1. - xsi[0];
        };

        //Computes the velocity shape functions
        shapeQuad.evaluateFem(xsi,phi_);
        
        //Computes the jacobian matrix
        getJacobianMatrix_FEM(xsi,Jac,ainv_);

        //Computes spatial derivatives
        getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

        //Interpolates velocity and its derivatives values
        getVelAndDerivatives_FEM(phi_,dphi_dx);    

        shapeBound.evaluateBoundaryFem(xsiB,phiB_);
        shapeBound.evaluateGradientBoundaryFem(xsiB,dphiB_);

        //computes the normal vector in the linear element
        double Tx=0.; double Ty = 0.;
        for (int i = 0; i < 3; i++){
            Tx += localNodesBoundary_[i][0] * dphiB_[0][i];
            Ty += localNodesBoundary_[i][1] * dphiB_[0][i];
        };

        double jacb_ = sqrt(Tx*Tx + Ty*Ty);
        
        n_vector[0] =  Ty/jacb_;
        n_vector[1] = -Tx/jacb_;

        shearStress[0][0] = 2. * visc_ * du_dx;
        shearStress[0][1] = visc_ * (du_dy + dv_dx);
        shearStress[1][0] = visc_ * (du_dy + dv_dx);
        shearStress[1][1] = 2. * visc_ * dv_dy;

        for (int i = 0; i < 2; i++){
            for (int j = 0; j < 2; j++){
                load_pressure[i] += -p_ * ident[i][j] * n_vector[j] * jacb_ * weightB;
                load_friction[i] += shearStress[i][j] * n_vector[j] * jacb_ * weightB;
            }
        }

        // moment += ((-p_ + shearStress[0][0] + shearStress[1][0]) * y_ 
        //          -(-p_ + shearStress[0][1] + shearStress[1][1]) * (x_ - 0.248792267683901))
        //          * jacb_ * weightB;
        per += jacb_ * weightB;
        
        index++;
    };

    aux_Per = per;
    aux_Mom = moment;

    pDForce = -load_pressure[0];
    pLForce = -load_pressure[1];
    
    fDForce = -load_friction[0];
    fLForce = -load_friction[1]; 

    dForce = pDForce + fDForce;
    lForce = pLForce + fLForce;


    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = 0; i < 1; ++i) delete [] dphiB_[i];
    delete [] dphiB_;   
    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] Jac[i];
    delete [] Jac;
  
    return;

};

//------------------------------------------------------------------------------
//--------------------COMPUTE BEZIER TRANSFORMATION MATRIX----------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getMatrixC(double **MatrixC){

    int degm = (*iparameters)[Npatch_] -> getDegree(0); // Degree in u direction
    int degn = (*iparameters)[Npatch_] -> getDegree(1); // Degree in v direction
    int npcm = (*iparameters)[Npatch_] -> getNcp(0); // Number of control points in u direction
    int npcn = (*iparameters)[Npatch_] -> getNcp(1); // Number of control points in v direction  
    double *uknot_ = (*iparameters)[Npatch_] -> getuKnot();
    double *vknot_ = (*iparameters)[Npatch_] -> getvKnot(); // Knots v direction
    int dimu = degm + npcm + 1;
    int dimv = degn + npcn + 1;
    double    MatrixCu[3][3] = {};
    double    MatrixCuNext[3][3] = {};
    double    MatrixCv[3][3] = {};
    double    MatrixCvNext[3][3] = {};
    double    alphasu[3] = {};
    double    alphasv[3] = {};

    int *inc_ = (*nodes_)[connect_[8]] -> getINC(); 

    int uind = inc_[0];
    int vind = inc_[1]; 

    // BEZIER EXTRACTION C U DIRECTION
    int au_ = degm;
    int bu_ = au_ + 1;
    int nb = 0;          //number of Bézier Elements
    int i_,mult,r,s,save;  //auxiliar     
    double numer,alpha;

    for (int i = 0 ; i <= degm; i++){
        MatrixCu[i][i] = 1.;
        MatrixCuNext[i][i] = 1.;
    };
    while (bu_ <= uind + 1){
        for (int i = 0 ; i <= degm; i++){
            for (int j = 0 ; j <= degm; j++){
                MatrixCu[i][j] = MatrixCuNext[i][j];
            };
        };
        for (int i = 0 ; i <= degm; i++){
            MatrixCuNext[i][i] = 1.;
        };
        i_ = bu_;
        while ((bu_ < dimu - 1) && (uknot_[bu_+ 1] == uknot_[bu_])) bu_ ++;
        mult = bu_ - i_ + 1;
        if (mult < degm){
            numer = uknot_[bu_] - uknot_[au_];
            for (int j = degm; j > mult; j--){
                alphasu[j-mult-1] = numer/(uknot_[au_+j] - uknot_[au_]);        
            };
            r = degm - mult; // Insert knot r times
            for (int j = 1; j<= r; j++){
                save = r - j;
                s = mult + j;
                for (int k = degm; k>= s; k--){
                    alpha = alphasu[k-s];
                    for (int i = 0; i <= degm; i++){
                        MatrixCu[i][k] = alpha * MatrixCu[i][k] + (1.0 - alpha) * MatrixCu[i][k-1];
                    };
                };           
                int cc = degm - j;
                for (int i=save; i<= j+save; i++){
                    MatrixCuNext[i][save] = MatrixCu[cc][degm];
                    cc++;
                };
            };
            nb = nb + 1;
            if (bu_ <= uind + 1){

                au_ = bu_;
                bu_ ++;
            };
        } else {
            nb = nb + 1;
            if (bu_ <= uind + 1){
                au_ = bu_;
            };
        };
    };

    // BEZIER EXTRACTION C V DIRECTION
    au_ = degn;
    bu_ = au_ + 1;
    nb = 0;          //number of Bézier Elements
    for (int i = 0 ; i <= degn; i++){
        MatrixCv[i][i] = 1.;
        MatrixCvNext[i][i] = 1.;
    };
    while (bu_ <= vind + 1){
        for (int i = 0 ; i <= degn; i++){
            for (int j = 0 ; j <= degn; j++){
                MatrixCv[i][j] = MatrixCvNext[i][j];
            };
        };
        for (int i = 0 ; i <= degn; i++){
            MatrixCvNext[i][i] = 1.;
        };
        i_ = bu_;
        while ((bu_ < dimv - 1) && (vknot_[bu_+ 1] == vknot_[bu_])) bu_ ++;
        mult = bu_ - i_ + 1;
        if (mult < degn){
            numer = vknot_[bu_] - vknot_[au_];
            for (int j = degn; j > mult; j--){
                alphasv[j-mult-1] = numer/(vknot_[au_+j] - vknot_[au_]);        
            };
            r = degn - mult; // Insert knot r times
            for (int j = 1; j<= r; j++){
                save = r - j;
                s = mult + j;
                for (int k = degn; k>= s; k--){
                    alpha = alphasv[k-s];
                    for (int i = 0; i <= degn; i++){
                        MatrixCv[i][k] = alpha * MatrixCv[i][k] + (1.0 - alpha) * MatrixCv[i][k-1];
                    };
                };        
                int cc = degn - j;
                for (int i=save; i<= j+save; i++){                      
                    MatrixCvNext[i][save] = MatrixCv[cc][degm];
                    cc++;
                };
            };
            nb = nb + 1;
            if (bu_ <= vind + 1){
                au_ = bu_;
                bu_ ++;
            };
        } else {
            nb = nb + 1;
            if (bu_ <= uind + 1){
                au_ = bu_;
            };
        };
    };

    for (int k = 0; k <= degn ; k++){
        for (int l = 0; l <= degn; l++){
            for (int i = 0; i <= degm ; i++){
                for (int j = 0; j <= degm; j++){
                    MatrixC[i + k*(degn+1)][j + l*(degm+1)] = MatrixCu[i][j] * MatrixCv[k][l];
                };
            };
        };
    }; 
    return;
}


template<>
void Element<2>::getInvMatrixC(double **MatrixCuInv, double **MatrixCvInv){

	// BEZIER EXTRACTION MATRIXC - PARAMETRIC DIRECTION u
    int degm = (*iparameters)[Npatch_] -> getDegree(0); // Degree 
    int npcm = (*iparameters)[Npatch_] -> getNcp(0); // Number of control points 
    int dimu = degm + npcm + 1; //Number of Knots 
    double *uknot_ = (*iparameters)[Npatch_] -> getuKnot();
    int *inc_ = (*nodes_)[connect_[8]] -> getINC(); 
    int uind = inc_[0]; // index of base fuction (parametric direction u) that starts in the index_ element
    double alphasu[3];
    double MatrixCu[3][3] = {};
    double MatrixCuNext[3][3] = {};
    
    // Algorithm 
    int au_ = degm;
    int bu_ = au_ + 1;
    int nb = 0;            //number of Bézier Elements
    int i_,mult,r,s,save;  //auxiliar     
    double numer,alpha;

    for (int i = 0 ; i <= degm; i++){
        MatrixCu[i][i] = 1.;
        MatrixCuNext[i][i] = 1.;
    };
    while (bu_ <= uind + 1){
        for (int i = 0 ; i <= degm; i++){
            for (int j = 0 ; j <= degm; j++){
                MatrixCu[i][j] = MatrixCuNext[i][j];
            };
        };
        for (int i = 0 ; i <= degm; i++){
            MatrixCuNext[i][i] = 1.;
        };
        i_ = bu_;
        while ((bu_ < dimu - 1) && (uknot_[bu_+ 1] == uknot_[bu_])) bu_ ++;
        mult = bu_ - i_ + 1;
        if (mult < degm){
            numer = uknot_[bu_] - uknot_[au_];
            for (int j = degm; j > mult; j--){
                alphasu[j-mult-1] = numer/(uknot_[au_+j] - uknot_[au_]);        
            };
            r = degm - mult; // Insert knot r times
            for (int j = 1; j<= r; j++){
                save = r - j;
                s = mult + j;
                for (int k = degm; k>= s; k--){
                    alpha = alphasu[k-s];
                    for (int i = 0; i <= degm; i++){
                        MatrixCu[i][k] = alpha * MatrixCu[i][k] + (1.0 - alpha) * MatrixCu[i][k-1];
                    };
                } ;              
                int cc = degm - j;
                for (int i=save; i<= j+save; i++){
                    MatrixCuNext[i][save] = MatrixCu[cc][degm];
                    cc++;
                };
            };
            nb = nb + 1;
            if (bu_ <= uind + 1){

                au_ = bu_;
                bu_ ++;
            };
        } else {
            nb = nb + 1;
            if (bu_ <= uind + 1){
                au_ = bu_;
            };
        };
    };


    //Computing Matrix Cu determinant
    double det_ = MatrixCu[0][0] * MatrixCu[1][1] * MatrixCu[2][2] + MatrixCu[0][1] * MatrixCu[1][2] * MatrixCu[2][0] + MatrixCu[0][2] * MatrixCu[1][0] * MatrixCu[2][1] - 
                  MatrixCu[0][1] * MatrixCu[1][0] * MatrixCu[2][2] - MatrixCu[0][0] * MatrixCu[1][2] * MatrixCu[2][1] - MatrixCu[0][2] * MatrixCu[1][1] * MatrixCu[2][0];

    //Inverting Matrix Cu
    MatrixCuInv[0][0] = (MatrixCu[1][1] * MatrixCu[2][2] - MatrixCu[1][2] * MatrixCu[2][1]) / det_;
    MatrixCuInv[1][0] = -(MatrixCu[1][0] * MatrixCu[2][2] - MatrixCu[1][2] * MatrixCu[2][0]) / det_;
    MatrixCuInv[2][0] = (MatrixCu[1][0] * MatrixCu[2][1] - MatrixCu[1][1] * MatrixCu[2][0]) / det_;
    MatrixCuInv[0][1] = -(MatrixCu[0][1] * MatrixCu[2][2] - MatrixCu[0][2] * MatrixCu[2][1]) / det_;
    MatrixCuInv[1][1] = (MatrixCu[0][0] * MatrixCu[2][2] - MatrixCu[0][2] * MatrixCu[2][0]) / det_;
    MatrixCuInv[2][1] = -(MatrixCu[0][0] * MatrixCu[2][1] - MatrixCu[0][1] * MatrixCu[2][0]) / det_;
    MatrixCuInv[0][2] = (MatrixCu[0][1] * MatrixCu[1][2] - MatrixCu[0][2] * MatrixCu[1][1]) / det_;
    MatrixCuInv[1][2] = -(MatrixCu[0][0] * MatrixCu[1][2] - MatrixCu[0][2] * MatrixCu[1][0]) / det_;
    MatrixCuInv[2][2] = (MatrixCu[0][0] * MatrixCu[1][1] - MatrixCu[0][1] * MatrixCu[1][0]) / det_;

    // BEZIER EXTRACTION MATRIXC - PARAMETRIC DIRECTION v
    int degn = (*iparameters)[Npatch_] -> getDegree(1); // Degree 
    int npcn = (*iparameters)[Npatch_] -> getNcp(1); // Number of control points
    int dimv = degn + npcn + 1;
    double *vknot_ = (*iparameters)[Npatch_] -> getvKnot();  
	int vind = inc_[1]; // index of base fuction (parametric direction v) that starts in the index_ element 
    double alphasv[3] = {};
    double MatrixCv[3][3] = {};
    double MatrixCvNext[3][3] = {};
    
    // BEZIER EXTRACTION C V DIRECTION
    au_ = degn;
    bu_ = au_ + 1;
    nb = 0;          //number of Bézier Elements
    for (int i = 0 ; i <= degn; i++){
        MatrixCv[i][i] = 1.;
        MatrixCvNext[i][i] = 1.;
    };    
    while (bu_ <= vind + 1){
        for (int i = 0 ; i <= degn; i++){
            for (int j = 0 ; j <= degn; j++){
                MatrixCv[i][j] = MatrixCvNext[i][j];
            };
        };
        for (int i = 0 ; i <= degn; i++){
            MatrixCvNext[i][i] = 1.;
        };
        i_ = bu_;
        while ((bu_ < dimv - 1) && (vknot_[bu_+ 1] == vknot_[bu_])) bu_ ++;
        mult = bu_ - i_ + 1;
        if (mult < degn){
            numer = vknot_[bu_] - vknot_[au_];
            for (int j = degn; j > mult; j--){
                alphasv[j-mult-1] = numer/(vknot_[au_+j] - vknot_[au_]);        
            };
            r = degn - mult; // Insert knot r times
            for (int j = 1; j<= r; j++){
                save = r - j;
                s = mult + j;
                for (int k = degn; k>= s; k--){
                    alpha = alphasv[k-s];
                    for (int i = 0; i <= degn; i++){
                        MatrixCv[i][k] = alpha * MatrixCv[i][k] + (1.0 - alpha) * MatrixCv[i][k-1];
                    };
                } ;           
                int cc = degn - j;
                for (int i=save; i<= j+save; i++){                        
                    MatrixCvNext[i][save] = MatrixCv[cc][degm];
                    cc++;
                };
            };
            nb = nb + 1;
            if (bu_ <= vind + 1){
                au_ = bu_;
                bu_ ++;
            };
        } else {
            nb = nb + 1;
            if (bu_ <= uind + 1){
                au_ = bu_;
            };
        };
    };

     //Computing Matrix Cv determinant
    det_ = MatrixCv[0][0] * MatrixCv[1][1] * MatrixCv[2][2] + MatrixCv[0][1] * MatrixCv[1][2] * MatrixCv[2][0] + MatrixCv[0][2] * MatrixCv[1][0] * MatrixCv[2][1] - 
           MatrixCv[0][1] * MatrixCv[1][0] * MatrixCv[2][2] - MatrixCv[0][0] * MatrixCv[1][2] * MatrixCv[2][1] - MatrixCv[0][2] * MatrixCv[1][1] * MatrixCv[2][0];

    //Inverting Matrix Cv
    MatrixCvInv[0][0] = (MatrixCv[1][1] * MatrixCv[2][2] - MatrixCv[1][2] * MatrixCv[2][1]) / det_;
    MatrixCvInv[1][0] = -(MatrixCv[1][0] * MatrixCv[2][2] - MatrixCv[1][2] * MatrixCv[2][0]) / det_;
    MatrixCvInv[2][0] = (MatrixCv[1][0] * MatrixCv[2][1] - MatrixCv[1][1] * MatrixCv[2][0]) / det_;
    MatrixCvInv[0][1] = -(MatrixCv[0][1] * MatrixCv[2][2] - MatrixCv[0][2] * MatrixCv[2][1]) / det_;
    MatrixCvInv[1][1] = (MatrixCv[0][0] * MatrixCv[2][2] - MatrixCv[0][2] * MatrixCv[2][0]) / det_;
    MatrixCvInv[2][1] = -(MatrixCv[0][0] * MatrixCv[2][1] - MatrixCv[0][1] * MatrixCv[2][0]) / det_;
    MatrixCvInv[0][2] = (MatrixCv[0][1] * MatrixCv[1][2] - MatrixCv[0][2] * MatrixCv[1][1]) / det_;
    MatrixCvInv[1][2] = -(MatrixCv[0][0] * MatrixCv[1][2] - MatrixCv[0][2] * MatrixCv[1][0]) / det_;
    MatrixCvInv[2][2] = (MatrixCv[0][0] * MatrixCv[1][1] - MatrixCv[0][1] * MatrixCv[1][0]) / det_;

};

//------------------------------------------------------------------------------
//------------------COMPUTES THE SUPG STABILIZATION PARAMETER-------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getNewParameterSUPG_FEM(double **Jac, double *phi_, double **dphi_dx) {
    
    double MatrixD[2][2],MatrixInvD[2][2],MatrixQh[2][2],MatrixInvQh[2][2],MatrixG[2][2],rr[2][2];   
    double r[2] = {}; double ua_[2] = {}; 
    double hrqd,tSUGN1_,tSUGN2_,tSUGN3_;

    double &dTime_ = parameters->getTimeStep();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();

    //Matrix D = polynomial order of the base functions
    MatrixD[0][0] = 2.;
    MatrixD[1][1] = 2.;
    MatrixD[1][0] = 0.;
    MatrixD[0][1] = 0.;

    //Inverse Matrix D
    double detMatrixD = MatrixD[0][0] * MatrixD[1][1] - MatrixD[0][1] * MatrixD[1][0];
    MatrixInvD[0][0] = (1./detMatrixD) * MatrixD[1][1];
    MatrixInvD[1][1] = (1./detMatrixD) * MatrixD[0][0];
    MatrixInvD[0][1] = -(1./detMatrixD) * MatrixD[0][1];
    MatrixInvD[1][0] = -(1./detMatrixD) * MatrixD[1][0];

    // Matrix Q "hat"
    MatrixQh[0][0] = Jac[0][0] * MatrixInvD[0][0] + Jac[0][1] * MatrixInvD[1][0];
    MatrixQh[0][1] = Jac[0][0] * MatrixInvD[0][1] + Jac[0][1] * MatrixInvD[1][1];
    MatrixQh[1][0] = Jac[1][0] * MatrixInvD[0][0] + Jac[1][1] * MatrixInvD[1][0];
    MatrixQh[1][1] = Jac[1][0] * MatrixInvD[0][1] + Jac[1][1] * MatrixInvD[1][1];

    //Matrix Inverse Q "hat"
    double detQh = MatrixQh[0][0] * MatrixQh[1][1] - MatrixQh[0][1] * MatrixQh[1][0];

    MatrixInvQh[0][0] = (1./detQh) *  MatrixQh[1][1];
    MatrixInvQh[1][1] = (1./detQh) *  MatrixQh[0][0];
    MatrixInvQh[0][1] = -(1./detQh) * MatrixQh[0][1];
    MatrixInvQh[1][0] = -(1./detQh) * MatrixQh[1][0];

    // Matrix G
    MatrixG[0][0] = MatrixInvQh[0][0] * MatrixInvQh[0][0] + MatrixInvQh[1][0] * MatrixInvQh[1][0];
    MatrixG[0][1] = MatrixInvQh[0][0] * MatrixInvQh[0][1] + MatrixInvQh[1][0] * MatrixInvQh[1][1];
    MatrixG[1][0] = MatrixInvQh[0][1] * MatrixInvQh[0][0] + MatrixInvQh[1][1] * MatrixInvQh[1][0];
    MatrixG[1][1] = MatrixInvQh[0][1] * MatrixInvQh[0][1] + MatrixInvQh[1][1] * MatrixInvQh[1][1];

    // Calculation hrqd
    for (int i = 0; i < 6; i++){

        double ua = (*nodes_)[connect_[i]] -> getVelocity(0);
        double va = (*nodes_)[connect_[i]] -> getVelocity(1);

        double uma = (*nodes_)[connect_[i]] -> getMeshVelocity(0);
        double vma = (*nodes_)[connect_[i]] -> getMeshVelocity(1);

        ua -= uma;
        va -= vma;

        ua_[0] += ua * phi_[i];
        ua_[1] += va * phi_[i];
        
        r[0] += sqrt(ua * ua + va * va) * dphi_dx[0][i];
        r[1] += sqrt(ua * ua + va * va) * dphi_dx[1][i];
   };

    double rNorm = sqrt(r[0]*r[0] + r[1]*r[1]);;
    double referpar = 0.000000000000000001;

    //physical unitary direction from gradient velocity 
    r[0] /= (rNorm + referpar);
    r[1] /= (rNorm + referpar);
   
    rr[0][0] = r[0] * r[0];
    rr[0][1] = r[0] * r[1];
    rr[1][0] = r[1] * r[0];
    rr[1][1] = r[1] * r[1];

    hrqd = 2./((sqrt(rr[0][0] * MatrixG[0][0] + rr[0][1] * MatrixG[0][1] + rr[1][0] * MatrixG[1][0] + rr[1][1] * MatrixG[1][1])) + referpar);

    //hmin e hmax
    //calculating the eigen-values from G
    double a_ = 1.;
    double b_ = -(MatrixG[0][0] + MatrixG[1][1]);
    double c_ = (MatrixG[0][0] * MatrixG[1][1] - MatrixG[0][1] * MatrixG[1][0]);

    double lambdaMax = (- b_ + sqrt(b_ * b_ - 4. * a_ * c_))/ (2. * a_);
    double lambdaMin = (- b_ - sqrt(b_ * b_ - 4. * a_ * c_))/ (2. * a_);

    double hmin = 2. / sqrt(lambdaMax);
    double hmax = 2. / sqrt(lambdaMin);

    if (hrqd < hmin) hrqd = hmin;
    if (hrqd > hmax) hrqd = hmax;

    //tSUGN1_ (-2)
    tSUGN1_ = (ua_[0] * ua_[0]) * MatrixG[0][0] +
              (ua_[0] * ua_[1]) * MatrixG[0][1] + 
              (ua_[1] * ua_[0]) * MatrixG[1][0] + 
              (ua_[1] * ua_[1]) * MatrixG[1][1];
    
    tSUGN2_ = dTime_/2.;

    //tSUGN3_ (-1)
    tSUGN3_ = (visc_/dens_)*(4./(hrqd*hrqd));

    //Computing tSUPG parameter
    tSUPG_ = 1./(sqrt(tSUGN1_ + (1./(tSUGN2_ * tSUGN2_)) + (tSUGN3_ * tSUGN3_)));

    tPSPG_ = tSUPG_;
   
    tLSIC_ = hrqd * hrqd / tSUPG_;
    
    return;

};

template<>
void Element<2>::getNewParameterSUPG_ISO(double **quadJacMat, double *phi_, double **dphi_dx) {

    double MatrixD[2][2] = {};
    double MatrixInvD[2][2],MatrixQh[2][2],MatrixInvQh[2][2],MatrixG[2][2],rr[2][2];   
    double r[2] = {}; double ua_[2] = {}; 
    double hrqd,tSUGN1_,tSUGN2_,tSUGN3_;
    double &dTime_ = parameters->getTimeStep();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity(); 

    
    //Computation of matrix D
    double **MatrixCuInv;
    MatrixCuInv = new double*[3];
    for (int i = 0; i < 3; ++i) MatrixCuInv[i] = new double[3];

    double **MatrixCvInv;
    MatrixCvInv = new double*[3];
    for (int i = 0; i < 3; ++i) MatrixCvInv[i] = new double[3];

    int degm = (*iparameters)[Npatch_] -> getDegree(0);
    int degn = (*iparameters)[Npatch_] -> getDegree(1);

    // Bezier parametric element lenght
    double dimxsi = 2.;
    double dimeta = 2.;

    getInvMatrixC(MatrixCuInv, MatrixCvInv);

    ublas::vector<double> Dxsi(3), Deta(3);

    Dxsi.clear();
    Deta.clear();

    for (int i = 0; i <= degm; i++){
        Dxsi(1) += i * (MatrixCuInv[i][1] - MatrixCuInv[i][0]);
        Dxsi(2) += i * (MatrixCuInv[i][2] - MatrixCuInv[i][1]);
    }

    Dxsi(1) *= (dimxsi/degm);
    Dxsi(2) *= (dimxsi/degm);

    for (int i = 0; i <= degn; i++){
        Deta(1) += i * (MatrixCvInv[i][1] - MatrixCvInv[i][0]);
        Deta(2) += i * (MatrixCvInv[i][2] - MatrixCvInv[i][1]);
    }

    Deta(1) *= (dimeta/degn);
    Deta(2) *= (dimeta/degn);

    typedef ublas::vector<double>::iterator itt;
    
    std::pair< itt , itt > results = minmax_element(Dxsi.begin(),Dxsi.end());
    std::pair< itt , itt > results2 = minmax_element(Deta.begin(),Deta.end());
    
    //RQD - MAX
    MatrixD[0][0] = *(results.second);
    MatrixD[1][1] = *(results2.second);

    //Inverse Matrix D
    double detMatrixD = MatrixD[0][0] * MatrixD[1][1] - MatrixD[0][1] * MatrixD[1][0];
    MatrixInvD[0][0] = (1./detMatrixD) * MatrixD[1][1];
    MatrixInvD[1][1] = (1./detMatrixD) * MatrixD[0][0];
    MatrixInvD[0][1] = -(1./detMatrixD) * MatrixD[0][1];
    MatrixInvD[1][0] = -(1./detMatrixD) * MatrixD[1][0];

    // Matrix Q "hat"
    MatrixQh[0][0] = quadJacMat[0][0] * MatrixInvD[0][0]+ quadJacMat[0][1] * MatrixInvD[1][0];
    MatrixQh[0][1] = quadJacMat[0][0] * MatrixInvD[0][1] + quadJacMat[0][1] * MatrixInvD[1][1];
    MatrixQh[1][0] = quadJacMat[1][0]* MatrixInvD[0][0]+ quadJacMat[1][1] * MatrixInvD[1][0];
    MatrixQh[1][1] = quadJacMat[1][0]* MatrixInvD[0][1] + quadJacMat[1][1] * MatrixInvD[1][1];

    //Matrix Inverse Q "hat"
    double detQh = MatrixQh[0][0]* MatrixQh[1][1] - MatrixQh[0][1] * MatrixQh[1][0];

    MatrixInvQh[0][0]= (1./detQh) *  MatrixQh[1][1];
    MatrixInvQh[1][1] = (1./detQh) *  MatrixQh[0][0];
    MatrixInvQh[0][1] = -(1./detQh) * MatrixQh[0][1];
    MatrixInvQh[1][0]= -(1./detQh) * MatrixQh[1][0];

    // Matrix G
    MatrixG[0][0]= MatrixInvQh[0][0]* MatrixInvQh[0][0]+ MatrixInvQh[1][0]* MatrixInvQh[1][0];
    MatrixG[0][1] = MatrixInvQh[0][0]* MatrixInvQh[0][1] + MatrixInvQh[1][0]* MatrixInvQh[1][1];
    MatrixG[1][0]= MatrixInvQh[0][1] * MatrixInvQh[0][0]+ MatrixInvQh[1][1] * MatrixInvQh[1][0];
    MatrixG[1][1] = MatrixInvQh[0][1] * MatrixInvQh[0][1] + MatrixInvQh[1][1] * MatrixInvQh[1][1];
   
    // Calculation hrqd
    for (int i = 0; i < 9; i++){

        double ua = (*nodes_)[connect_[i]] -> getVelocity(0);
        double va = (*nodes_)[connect_[i]] -> getVelocity(1);

        double uma = (*nodes_)[connect_[i]] -> getMeshVelocity(0);
        double vma = (*nodes_)[connect_[i]] -> getMeshVelocity(1);

        ua -= uma;
        va -= vma;

        ua_[0]+= ua * phi_[i];
        ua_[1] += va * phi_[i];
        
        r[0]+= sqrt(ua * ua + va * va) * dphi_dx[0][i];
        r[1] += sqrt(ua * ua + va * va) * dphi_dx[1][i];

    };

    double rNorm = sqrt(r[0]*r[0] + r[1]*r[1]);
    double referpar = 0.000000000000000001;

    r[0] /= (rNorm + referpar);
    r[1] /= (rNorm + referpar);

    rr[0][0]= r[0]* r[0];
    rr[0][1] = r[0]* r[1];
    rr[1][0]= r[1] * r[0];
    rr[1][1] = r[1] * r[1];

    hrqd = 2./((sqrt(rr[0][0]* MatrixG[0][0]+ rr[0][1] * MatrixG[0][1] + rr[1][0]* MatrixG[1][0]+ rr[1][1] * MatrixG[1][1])) + referpar);
    
    //hmin e hmax
    //calculating the eigen-values from G
    double a_ = 1.;
    double b_ = -(MatrixG[0][0]+ MatrixG[1][1]);
    double c_ = (MatrixG[0][0]* MatrixG[1][1] - MatrixG[0][1] * MatrixG[1][0]);

    double square = (b_ * b_ - 4. * a_ * c_);

    if (square < 0.000001) square = 0.;

    double lambdaMax = (- b_ + sqrt(square))/ (2. * a_);
    double lambdaMin = (- b_ - sqrt(square))/ (2. * a_);

    double hmin = 2. / sqrt(lambdaMax);
    double hmax = 2. / sqrt(lambdaMin);

    if (hrqd < hmin) hrqd = hmin;
    if (hrqd > hmax) hrqd = hmax;


    //tSUGN1_ (-2)
    tSUGN1_ = (ua_[0]* ua_[0]) * MatrixG[0][0]+
              (ua_[0]* ua_[1]) * MatrixG[0][1] + 
              (ua_[1] * ua_[0]) * MatrixG[1][0]+ 
              (ua_[1] * ua_[1]) * MatrixG[1][1];
    
    tSUGN2_ = dTime_/2.;

    //tSUGN3_ (-1)
    tSUGN3_ = (visc_/dens_)*(4./(hrqd*hrqd));

    //Computing tSUPG parameter
    tSUPG_ = 1./(sqrt(tSUGN1_ + (1./(tSUGN2_ * tSUGN2_)) + (tSUGN3_ * tSUGN3_)));

    tPSPG_ = tSUPG_;
   
    tLSIC_ = hrqd * hrqd / tSUPG_;

    for (int i = 0; i < 3; ++i) delete [] MatrixCuInv[i];
    delete [] MatrixCuInv;
    for (int i = 0; i < 3; ++i) delete [] MatrixCvInv[i];
    delete [] MatrixCvInv;
    
    return;
};


template<>
void Element<2>::getParameterArlequin_FEM(double *phi_, double **dphi_dx){

    double &dTime_ = parameters->getTimeStep();
    double &k1 = parameters->getArlequinK1();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();

    double hRGN_ = 0.;
    double hUGN_ = 0.;
    double tSUGN1_,tSUGN2_,tSUGN3_;  

    double r[2] = {}; 
    double s[2] = {};

    double u_ = 0.;
    double v_ = 0.;
    
    for (int i = 0; i < 6; i++){

        double ua = (*nodes_)[connect_[i]] -> getVelocity(0);
        double va = (*nodes_)[connect_[i]] -> getVelocity(1);

        double uma = (*nodes_)[connect_[i]] -> getMeshVelocity(0);
        double vma = (*nodes_)[connect_[i]] -> getMeshVelocity(1);

        ua -= uma;
        va -= vma;

        u_ += ua * phi_[i];
        v_ += va * phi_[i];
        
    };

    double uNorm = sqrt(u_*u_ + v_*v_);
     if(uNorm > 1.e-10){
        s[0] = u_ / uNorm;
        s[1] = v_ / uNorm;
    }else{
        s[0] = 1. / sqrt(2.);
        s[1] = 1. / sqrt(2.);
    };

    for (int i = 0; i < 6; i++){

        double ua = (*nodes_)[connect_[i]] -> getVelocity(0);
        double va = (*nodes_)[connect_[i]] -> getVelocity(1);

        double uma = (*nodes_)[connect_[i]] -> getMeshVelocity(0);
        double vma = (*nodes_)[connect_[i]] -> getMeshVelocity(1);

        ua -= uma;
        va -= vma;

        r[0] += sqrt(ua * ua + va * va) * dphi_dx[0][i];
        r[1] += sqrt(ua * ua + va * va) * dphi_dx[1][i];
        
    };

    double rNorm = sqrt(r[0]*r[0] + r[1]*r[1]);;
    if (rNorm >= 1.e-10){
        r[0] /= rNorm;
        r[1] /= rNorm;
    }else{
        r[0] = 1. / sqrt(2.);
        r[1] = 1. / sqrt(2.);
    };


    for (int i = 0; i < 6; i++){
        hRGN_ += fabs(r[0] * dphi_dx[0][i] + r[1] * dphi_dx[1][i]);
        hUGN_ += fabs(s[0] * dphi_dx[0][i] + s[1] * dphi_dx[1][i]);        
    };

    if (hRGN_ >= 1.e-10){
        hRGN_ = 2. / hRGN_;
    }else{
        hRGN_ = 2. / 1.e-10;
    };

    if (hUGN_ >= 1.e-10){
        hUGN_ = 2. / hUGN_;
    }else{
        hUGN_ = 2. / 1.e-10;
    };  


    if (uNorm >= 1.e-10){
        tSUGN1_ = hUGN_ / (2. * uNorm);
    }else{
        tSUGN1_ = hUGN_ / 2.e-10;
    };
              
    tSUGN2_ = dTime_ / 2.;

    tSUGN3_ = hRGN_ * hRGN_ / (4. * visc_ / dens_);

   

    if (fabs(tSUGN1_) <= 1.e-10) tSUGN1_ = 1.e-10;
    if (fabs(tSUGN3_) <= 1.e-10) tSUGN3_ = 1.e-10;


 
    //Computing tARLQ parameter
    tSUPG_ = 1. / sqrt(1. / (tSUGN1_ * tSUGN1_) + 
                       1. / (tSUGN2_ * tSUGN2_) + 
                       1. / (tSUGN3_ * tSUGN3_));


    tARLQ_ = -1. * k1 * tSUPG_ * 1.e-2;

    tPSPG_ = tSUPG_;


    return;
}


template<>
void Element<2>::getParameterArlequin2_FEM(double *phi_, double *phiC_, double **dphi_dx, double ***ddphi_dx){

    
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &visc_ = parameters->getViscosity();

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
    double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * vmeshprev_);

    double duna_dx = (alpha_f * du_dx + (1. - alpha_f) * duprev_dx);
    double duna_dy = (alpha_f * du_dy + (1. - alpha_f) * duprev_dy);
    double dvna_dx = (alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx);
    double dvna_dy = (alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy);

    double dumeshna_dx = (alpha_f * dumesh_dx + (1. - alpha_f) * dumeshprev_dx);
    double dumeshna_dy = (alpha_f * dumesh_dy + (1. - alpha_f) * dumeshprev_dy);
    double dvmeshna_dx = (alpha_f * dvmesh_dx + (1. - alpha_f) * dvmeshprev_dx);
    double dvmeshna_dy = (alpha_f * dvmesh_dy + (1. - alpha_f) * dvmeshprev_dy);

    double dduna_dxdx = (alpha_f * ddu_dxdx + (1. - alpha_f) * dduprev_dxdx);
    double dduna_dxdy = (alpha_f * ddu_dxdy + (1. - alpha_f) * dduprev_dxdy);
    double dduna_dydx = (alpha_f * ddu_dydx + (1. - alpha_f) * dduprev_dydx);
    double dduna_dydy = (alpha_f * ddu_dydy + (1. - alpha_f) * dduprev_dydy);
    double ddvna_dxdx = (alpha_f * ddv_dxdx + (1. - alpha_f) * ddvprev_dxdx);
    double ddvna_dxdy = (alpha_f * ddv_dxdy + (1. - alpha_f) * ddvprev_dxdy);
    double ddvna_dydx = (alpha_f * ddv_dydx + (1. - alpha_f) * ddvprev_dydx);
    double ddvna_dydy = (alpha_f * ddv_dydy + (1. - alpha_f) * ddvprev_dydy);

    double daxm_dx = alpha_m * dax_dx + (1. - alpha_m) * daxprev_dx;
    double daxm_dy = alpha_m * dax_dy + (1. - alpha_m) * daxprev_dy;
    double daym_dx = alpha_m * day_dx + (1. - alpha_m) * dayprev_dx;
    double daym_dy = alpha_m * day_dy + (1. - alpha_m) * dayprev_dy;

    double lambda1[12][12] = {};
    double lambda0[12][18] = {};
    double convec[12] = {};
    double iner[12] = {};
    double visc[12] = {};
    double press[12] = {};
    double gradlambda[12] = {};

    double WJ = weight_ * djac_;

    for (int i = 0; i < 6; i++){

            //convection
            convec[2*i]   += (((dphi_dx[0][i] * (dduna_dxdx * una_ + dduna_dxdy * vna_)) + 
                               (dphi_dx[1][i] * (dduna_dydx * una_ + dduna_dydy * vna_))) +
                              ((dphi_dx[0][i] * (duna_dx * (duna_dx - dumeshna_dx) + duna_dy * (dvna_dx - dvmeshna_dx))) + 
                               (dphi_dx[1][i] * (duna_dx * (duna_dy - dumeshna_dy) + duna_dy * (dvna_dy - dvmeshna_dy)))) + 
                               (dphi_dx[0][i] * una_ * duna_dx + dphi_dx[0][i] * vna_ * duna_dy)) * WJ;
            convec[2*i+1] += (((dphi_dx[0][i] * (ddvna_dxdx * una_ + ddvna_dxdy * vna_)) + 
                               (dphi_dx[1][i] * (ddvna_dydx * una_ + ddvna_dydy * vna_))) +
                              ((dphi_dx[0][i] * (dvna_dx * (duna_dx - dumeshna_dx) + dvna_dy * (dvna_dx - dvmeshna_dx))) + 
                               (dphi_dx[1][i] * (dvna_dx * (duna_dy - dumeshna_dy)  + dvna_dy * (dvna_dy - dvmeshna_dy)))) +
                               (dphi_dx[1][i] * una_ * dvna_dx + dphi_dx[1][i] * vna_ * dvna_dy)) * WJ;
            
            //inercia
            iner[2*i]   += (dphi_dx[0][i] * daxm_dx + dphi_dx[1][i] * daxm_dy) * WJ;
            iner[2*i+1] += (dphi_dx[0][i] * daym_dx + dphi_dx[1][i] * daym_dy) * WJ;

            //viscosity
            visc[2*i]   += ((ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * (2. * dduna_dxdx + dduna_dydy) +
                             (ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * ddvna_dxdy) * WJ * visc_;
            visc[2*i+1] += ((ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * dduna_dydx +
                             (ddphi_dx[0][0][i] + ddphi_dx[1][1][i]) * (2. * ddvna_dydy + ddvna_dxdx)) * WJ * visc_;

            //pressure 
            press[2*i]   += (dphi_dx[0][i] * ddp_dxdx + dphi_dx[1][i] * ddp_dydx) * WJ;
            press[2*i+1] += (dphi_dx[0][i] * ddp_dxdy + dphi_dx[1][i] * ddp_dydy) * WJ;

            gradlambda[2*i]   += dphi_dx[0][i] * lamx_dx + dphi_dx[1][i] * lamx_dy;
            gradlambda[2*i+1] += dphi_dx[0][i] * lamy_dx + dphi_dx[1][i] * lamy_dy;
        
        for (int j = 0; j < 6; j++){
            //lagrange multipliers
            lambda1[2*i][2*j]     += phi_[i] * phi_[j] * WJ;
            lambda1[2*i+1][2*j+1] += phi_[i] * phi_[j] * WJ;
        };

    };

    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 9; j++){
            //lagrange multipliers
            lambda0[2*i][2*j] += phi_[i] * phiC_[j] * WJ;
            lambda0[2*i+1][2*j+1] += phi_[i] * phiC_[j] * WJ;
        };
    };

    double normLambda1[12] = {};
    double normLambda0[12] = {};


    double uL[12], uG[18];
    for (int i = 0; i < 6; i++){
        uL[2*i] = (*nodes_)[connect_[i]] -> getVelocity(0);
        uL[2*i+1] = (*nodes_)[connect_[i]] -> getVelocity(1);
    };
    for (int i = 0; i < 9; i++){
        uG[2*i] = (*nodesC_)[connectC_[i]] -> getVelocity(0);
        uG[2*i+1] = (*nodesC_)[connectC_[i]] -> getVelocity(1);
    };

    for (int i = 0; i < 12 ; i++){
        for (int j = 0; j < 12; j++){
            normLambda1[i] += lambda1[i][j] * uL[j];
        };
    };

    for (int i = 0; i < 12 ; i++){
        for (int j = 0; j < 18; j++){
            normLambda0[i] += lambda0[i][j] * uG[j];
        };
    };


    double vecnormLambda0 = 0.;
    double vecnormLambda1 = 0.;
    for (int i = 0; i < 12; i++){
        vecnormLambda0 += normLambda0[i] * normLambda0[i];
        vecnormLambda1 += normLambda1[i] * normLambda1[i];
    };


    double normConvec = 0.;
    double normIner = 0.;
    double normVisc = 0.;
    double normPress = 0.;
    double normGradLambda = 0.;

    for (int i = 0; i < 12 ; i++){

        normConvec += convec[i] * convec[i];
        normIner += iner[i] * iner[i];
        normVisc += visc[i] * visc[i];
        normPress += press[i] * press[i];
        normGradLambda += gradlambda[i] * gradlambda[i];
    };

    normConvec = sqrt(normConvec);
    normIner = sqrt(normIner);
    normVisc = sqrt(normVisc);
    normPress = sqrt(normPress);
    normGradLambda = sqrt(normGradLambda);
    vecnormLambda0 = sqrt(vecnormLambda0);
    vecnormLambda1 = sqrt(vecnormLambda1);

    if (fabs(normIner) <= 1.e-10) normIner = 1.e-10;
    if (fabs(normVisc) <= 1.e-10) normVisc = 1.e-10;
    if (fabs(normPress) <= 1.e-10) normPress = 1.e-10;
    if (fabs(normConvec) <= 1.e-10) normConvec = 1.e-10;
    if (fabs(normGradLambda) <= 1.e-10) normGradLambda = 1.e-10;
    
    double tA1 = std::fabs(vecnormLambda0)/std::fabs(normConvec);
    double tA2 = std::fabs(vecnormLambda0)/std::fabs(normConvec);
    double tB1 = std::fabs(vecnormLambda1)/std::fabs(normIner);
    double tB2 = std::fabs(vecnormLambda0)/std::fabs(normIner);
    double tC1 = std::fabs(vecnormLambda1)/std::fabs(normVisc);
    double tC2 = std::fabs(vecnormLambda0)/std::fabs(normVisc);
    double tD1 = std::fabs(vecnormLambda1)/std::fabs(normPress);
    double tD2 = std::fabs(vecnormLambda0)/std::fabs(normPress);
    double tE1 = std::fabs(vecnormLambda1)/std::fabs(normGradLambda);
    double tE2 = std::fabs(vecnormLambda0)/std::fabs(normGradLambda);

    double tA =  1. / sqrt(1. / (tA1 * tA1) + 
                       1. / (tA2 * tA2));

    double tB =  1. / sqrt(1. / (tB1 * tB1) + 
                       1. / (tB2 * tB2));

    double tC =  1. / sqrt(1. / (tC1 * tC1) + 
                       1. / (tC2 * tC2));

    double tD =  1. / sqrt(1. / (tD1 * tD1) + 
                       1. / (tD2 * tD2));

    double tE =  1. / sqrt(1. / (tE1 * tE1) + 
                       1. / (tE2 * tE2));


    tARLQ_ = -1. / sqrt(1. / (tA * tA) + 
                       1. / (tB * tB) +
                       1. / (tC * tC) +
                       1. / (tE * tE) + 
                       1. / (tD * tD));


};



template<>
void Element<2>::getParameterArlequin_ISO(double *phi_, double **dphi_dx){



    double &dTime_ = parameters->getTimeStep();
    double &k1 = parameters->getArlequinK1();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();

    double hRGN_ = 0.;
    double hUGN_ = 0.;
    double tSUGN1_,tSUGN2_,tSUGN3_;  

    double r[2] = {}; 
    double s[2] = {};

    double u_ = 0.;
    double v_ = 0.;
    
    for (int i = 0; i < 9; i++){

        double ua = (*nodes_)[connect_[i]] -> getVelocity(0);
        double va = (*nodes_)[connect_[i]] -> getVelocity(1);

        double uma = (*nodes_)[connect_[i]] -> getMeshVelocity(0);
        double vma = (*nodes_)[connect_[i]] -> getMeshVelocity(1);

        ua -= uma;
        va -= vma;

        u_ += ua * phi_[i];
        v_ += va * phi_[i];
        
    };

    double uNorm = sqrt(u_*u_ + v_*v_);
     if(uNorm > 1.e-10){
        s[0] = u_ / uNorm;
        s[1] = v_ / uNorm;
    }else{
        s[0] = 1. / sqrt(2.);
        s[1] = 1. / sqrt(2.);
    };

    for (int i = 0; i < 9; i++){

        double ua = (*nodes_)[connect_[i]] -> getVelocity(0);
        double va = (*nodes_)[connect_[i]] -> getVelocity(1);

        double uma = (*nodes_)[connect_[i]] -> getMeshVelocity(0);
        double vma = (*nodes_)[connect_[i]] -> getMeshVelocity(1);

        ua -= uma;
        va -= vma;

        r[0] += sqrt(ua * ua + va * va) * dphi_dx[0][i];
        r[1] += sqrt(ua * ua + va * va) * dphi_dx[1][i];
        
    };

    double rNorm = sqrt(r[0]*r[0] + r[1]*r[1]);;
    if (rNorm >= 1.e-10){
        r[0] /= rNorm;
        r[1] /= rNorm;
    }else{
        r[0] = 1. / sqrt(2.);
        r[1] = 1. / sqrt(2.);
    };


    for (int i = 0; i < 9; i++){
        hRGN_ += fabs(r[0] * dphi_dx[0][i] + r[1] * dphi_dx[1][i]);
        hUGN_ += fabs(s[0] * dphi_dx[0][i] + s[1] * dphi_dx[1][i]);        
    };

    if (hRGN_ >= 1.e-10){
        hRGN_ = 2. / hRGN_;
    }else{
        hRGN_ = 2. / 1.e-10;
    };

    if (hUGN_ >= 1.e-10){
        hUGN_ = 2. / hUGN_;
    }else{
        hUGN_ = 2. / 1.e-10;
    };  


    if (uNorm >= 1.e-10){
        tSUGN1_ = hUGN_ / (2. * uNorm);
    }else{
        tSUGN1_ = hUGN_ / 2.e-10;
    };
              
    tSUGN2_ = dTime_ / 2.;

    tSUGN3_ = hRGN_ * hRGN_ / (4. * visc_ / dens_);

   

    if (fabs(tSUGN1_) <= 1.e-10) tSUGN1_ = 1.e-10;
    if (fabs(tSUGN3_) <= 1.e-10) tSUGN3_ = 1.e-10;


 
    //Computing tARLQ parameter
    double tsupg = 1. / sqrt(1. / (tSUGN1_ * tSUGN1_) + 
                       1. / (tSUGN2_ * tSUGN2_) + 
                       1. / (tSUGN3_ * tSUGN3_));


    tARLQ_ = -1. * k1 * tsupg * 1.e-2;


    return;
}



//------------------------------------------------------------------------------
//-----------------------------ELEMENT LOCAL MATRIX-----------------------------
//------------------------------------------------------------------------------

template<>
void Element<2>::getElemMatrix_FEM(int &index,double *phi_, double **dphi_dx,double **jacobianNRMatrix){
    
    double &dTime_ = parameters->getTimeStep();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &gamma = parameters->getGamma();

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
    double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * vmeshprev_);

    double duna_dx = (alpha_f * du_dx + (1. - alpha_f) * duprev_dx);
    double duna_dy = (alpha_f * du_dy + (1. - alpha_f) * duprev_dy);
    double dvna_dx = (alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx);
    double dvna_dy = (alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy);

    double wna_ = alpha_f * intPointWeightFunction_FEM[index] + (1. - alpha_f) * intPointWeightFunctionPrev_FEM[index];
    double WJ = weight_ * djac_ * wna_;

    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){

            double wSUPGi = una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i];
            double wSUPGj = una_ * dphi_dx[0][j] + vna_ * dphi_dx[1][j];
            
            //Mass matrix (use for both directions) 
            double mM =  phi_[i] * phi_[j] * dens_ * alpha_m + 
                         wSUPGi * phi_[j] * tSUPG_ * dens_ * alpha_m;
                                 
            //Difusion matrix (viscosity)
            double Kxx = (2. * dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) 
                         * visc_* alpha_f * gamma * dTime_;
            double Kxy = dphi_dx[1][i] * dphi_dx[0][j] * visc_* alpha_f * gamma * dTime_;
            double Kyx = dphi_dx[0][i] * dphi_dx[1][j] * visc_ * alpha_f * gamma * dTime_;
            double Kyy = (2. * dphi_dx[1][i] * dphi_dx[1][j] + dphi_dx[0][i] * dphi_dx[0][j])
                         * visc_* alpha_f * gamma * dTime_;
            
            //Convection matrixes
            double Cxx = (dphi_dx[0][j] * una_ + dphi_dx[1][j] * vna_) 
                        * phi_[i] * dens_ * alpha_f * gamma * dTime_
                        + wSUPGi * wSUPGj * tSUPG_ * dens_* alpha_f * gamma * dTime_;
            double Cyy  = Cxx;
            
            double Cuu = phi_[i] * duna_dx * phi_[j] * dens_* alpha_f * gamma * dTime_;
            double Cuv = phi_[i] * duna_dy * phi_[j] * dens_ * alpha_f * gamma * dTime_;
            double Cvu = phi_[i] * dvna_dx * phi_[j] * dens_ * alpha_f * gamma * dTime_;
            double Cvv = phi_[i] * dvna_dy * phi_[j] * dens_ * alpha_f * gamma * dTime_;

            //Stabilization LSIC matrix          
            double KLSxx = dphi_dx[0][i] * dphi_dx[0][j] * tLSIC_ * dens_ * alpha_f * gamma * dTime_;
            double KLSxy = dphi_dx[0][i] * dphi_dx[1][j] * tLSIC_ * dens_ * alpha_f * gamma * dTime_;
            double KLSyx = dphi_dx[1][i] * dphi_dx[0][j] * tLSIC_ * dens_ * alpha_f * gamma * dTime_;
            double KLSyy = dphi_dx[1][i] * dphi_dx[1][j] * tLSIC_ * dens_ * alpha_f * gamma * dTime_;

            jacobianNRMatrix[2*i  ][2*j  ] += (mM + Kxx + KLSxx + Cxx + Cuu) * WJ;
            jacobianNRMatrix[2*i+1][2*j+1] += (mM + Kyy + KLSyy + Cyy + Cvv) * WJ;
            jacobianNRMatrix[2*i  ][2*j+1] += (Kxy + Cuv + KLSxy) * WJ;
            jacobianNRMatrix[2*i+1][2*j  ] += (Kyx + Cvu + KLSyx) * WJ; 
            
            //multipy pressure direction x and y
            double QSUPGx = - (dphi_dx[0][i] * phi_[j]) + 
                              wSUPGi * dphi_dx[0][j] * tSUPG_;
            double QSUPGy = - (dphi_dx[1][i] * phi_[j]) +
                               wSUPGi *dphi_dx[1][j] * tSUPG_;
            
            //multiply velocity direction x and y
            double Qx = dphi_dx[0][j] * phi_[i] * alpha_f * gamma * dTime_;
            double Qy = dphi_dx[1][j] * phi_[i] * alpha_f * gamma * dTime_;                

            jacobianNRMatrix[12+i][2*j  ] += Qx * WJ;
            jacobianNRMatrix[12+i][2*j+1] += Qy * WJ;
            jacobianNRMatrix[2*i  ][12+j] += QSUPGx * WJ;
            jacobianNRMatrix[2*i+1][12+j] += QSUPGy * WJ;

            //PSPG stabilization matrixes
            double Hx = dphi_dx[0][i] * phi_[j] * tPSPG_ * alpha_m;
            double Hy = dphi_dx[1][i] * phi_[j] * tPSPG_* alpha_m;
            
            double Gx = dphi_dx[0][i] * wSUPGj * tPSPG_ * alpha_f * gamma * dTime_;
            double Gy = dphi_dx[1][i] * wSUPGj * tPSPG_ * alpha_f * gamma * dTime_;

            double Q = (dphi_dx[0][i] * dphi_dx[0][j] + 
                        dphi_dx[1][i] * dphi_dx[1][j]) * tPSPG_ / (dens_);         

            jacobianNRMatrix[12+i][2*j  ] += (Hx + Gx) * WJ;
            jacobianNRMatrix[12+i][2*j+1] += (Hy + Gy) * WJ;
            jacobianNRMatrix[12+i][12+j]  +=  Q * WJ;              
        };
    };


    // //Stokes Problem

    // for (int i = 0; i < 6; i++){
    //     for (int j = 0; j < 6; j++){

    //         //Difusion matrix (viscosity)
    //         double Kxx = (2. * dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) 
    //                      * visc_* alpha_f * gamma * dTime_;
    //         double Kxy = dphi_dx[1][i] * dphi_dx[0][j] * visc_* alpha_f * gamma * dTime_;
    //         double Kyx = dphi_dx[0][i] * dphi_dx[1][j] * visc_ * alpha_f * gamma * dTime_;
    //         double Kyy = (2. * dphi_dx[1][i] * dphi_dx[1][j] + dphi_dx[0][i] * dphi_dx[0][j])
    //                      * visc_* alpha_f * gamma * dTime_;
            
    //         jacobianNRMatrix[2*i  ][2*j  ] += (Kxx) * WJ;
    //         jacobianNRMatrix[2*i+1][2*j+1] += (Kyy) * WJ;
    //         jacobianNRMatrix[2*i  ][2*j+1] += (Kxy) * WJ;
    //         jacobianNRMatrix[2*i+1][2*j  ] += (Kyx) * WJ; 
            
    //         //multipy pressure direction x and y
    //         double QSUPGx = - (dphi_dx[0][i] * phi_[j]);
    //         double QSUPGy = - (dphi_dx[1][i] * phi_[j]);
            
    //         //multiply velocity direction x and y
    //         double Qx = dphi_dx[0][j] * phi_[i] * alpha_f * gamma * dTime_;
    //         double Qy = dphi_dx[1][j] * phi_[i] * alpha_f * gamma * dTime_;                

    //         jacobianNRMatrix[12+i][2*j  ] += Qx * WJ;
    //         jacobianNRMatrix[12+i][2*j+1] += Qy * WJ;
    //         jacobianNRMatrix[2*i  ][12+j] += QSUPGx * WJ;
    //         jacobianNRMatrix[2*i+1][12+j] += QSUPGy * WJ;

    //         double Q = (dphi_dx[0][i] * dphi_dx[0][j] + 
    //                     dphi_dx[1][i] * dphi_dx[1][j]) * tPSPG_ / (dens_);         

    //         jacobianNRMatrix[12+i][12+j]  +=  Q * WJ;              
    //     };
    // };

    return;
};

template<>
void Element<2>::getElemMatrix_ISO(int &index, double *phi_, double **dphi_dx, double **jacobianNRMatrix){
    
    double &dTime_ = parameters->getTimeStep();
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &gamma = parameters->getGamma();

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
    double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * umeshprev_);

    double duna_dx = (alpha_f * du_dx + (1. - alpha_f) * duprev_dx);
    double duna_dy = (alpha_f * du_dy + (1. - alpha_f) * duprev_dy);
    double dvna_dx = (alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx);
    double dvna_dy = (alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy);

    
    double wna_ = alpha_f * intPointWeightFunction_ISO[index] + (1. - alpha_f) * intPointWeightFunctionPrev_ISO[index];

    double WJ = weight_ * djac_ * wna_;

    for (int i = 0; i < 9; i++){
        for (int j = 0; j < 9; j++){


            double wSUPGi = una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i];
            double wSUPGj = una_ * dphi_dx[0][j] + vna_ * dphi_dx[1][j];
            

            //Mass matrix (use for both directions)
            double mM =  phi_[i] * phi_[j] * dens_ * alpha_m + 
                         wSUPGi * phi_[j] * tSUPG_ * dens_ * alpha_m;
                       
            //Difusion matrix (viscosity)
            double Kxx = (2. * dphi_dx[0][i] * dphi_dx[0][j] + 
                         dphi_dx[1][i] * dphi_dx[1][j]) * visc_* alpha_f * gamma * dTime_;
            double Kxy = dphi_dx[1][i] * dphi_dx[0][j] * visc_ * alpha_f * gamma * dTime_;
            double Kyx = dphi_dx[0][i] * dphi_dx[1][j] * visc_ * alpha_f * gamma * dTime_;
            double Kyy = (2. * dphi_dx[1][i] * dphi_dx[1][j] + 
                         dphi_dx[0][i] * dphi_dx[0][j]) * visc_* alpha_f * gamma * dTime_;
            
            // //Convection matrix
            double Cxx = (dphi_dx[0][j] * una_ + 
                          dphi_dx[1][j] * vna_) * phi_[i] * dens_* alpha_f * gamma * dTime_
                        + wSUPGi * wSUPGj * tSUPG_ * dens_* alpha_f * gamma * dTime_;
            
            double Cyy  = Cxx;
            
            double Cuu = phi_[i] * duna_dx * phi_[j] * dens_* alpha_f * gamma * dTime_;           
            double Cuv = phi_[i] * duna_dy * phi_[j] * dens_* alpha_f * gamma * dTime_; 
            double Cvu = phi_[i] * dvna_dx * phi_[j] * dens_* alpha_f * gamma * dTime_;     
            double Cvv = phi_[i] * dvna_dy * phi_[j] * dens_* alpha_f * gamma * dTime_; 
                        
            double KLSxx = dphi_dx[0][i] * dphi_dx[0][j] * tLSIC_ * dens_* alpha_f * gamma * dTime_;
            double KLSxy = dphi_dx[0][i] * dphi_dx[1][j] * tLSIC_ * dens_* alpha_f * gamma * dTime_;
            double KLSyx = dphi_dx[1][i] * dphi_dx[0][j] * tLSIC_ * dens_* alpha_f * gamma * dTime_;
            double KLSyy = dphi_dx[1][i] * dphi_dx[1][j] * tLSIC_ * dens_* alpha_f * gamma * dTime_;


            jacobianNRMatrix[2*i  ][2*j  ] += (mM + Kxx + KLSxx + Cxx + Cuu) * WJ;
            jacobianNRMatrix[2*i+1][2*j+1] += (mM + Kyy + KLSyy + Cyy + Cvv) * WJ;
            jacobianNRMatrix[2*i  ][2*j+1] += (Kxy + Cuv + KLSxy) * WJ;
            jacobianNRMatrix[2*i+1][2*j  ] += (Kyx + Cvu + KLSyx) * WJ; 


            double QSUPGx = - (dphi_dx[0][i] * phi_[j]) + 
                            ((dphi_dx[0][i] * una_ + 
                            dphi_dx[1][i] * vna_) * 
                            dphi_dx[0][j] * tSUPG_);
            //multiply pressure direction y
            double QSUPGy = - (dphi_dx[1][i] * phi_[j]) +
                            ((dphi_dx[0][i] * una_ + 
                            dphi_dx[1][i] * vna_) *   
                            dphi_dx[1][j] * tSUPG_);

            //multiply velocity direction x
            double Qx = dphi_dx[0][j] * phi_[i] * alpha_f * gamma * dTime_;
            //multiply velocity direction y
            double Qy = dphi_dx[1][j] * phi_[i] * alpha_f * gamma * dTime_;

            jacobianNRMatrix[18+i][2*j  ] += Qx * WJ;
            jacobianNRMatrix[18+i][2*j+1] += Qy * WJ;
            jacobianNRMatrix[2*i  ][18+j] += QSUPGx * WJ;
            jacobianNRMatrix[2*i+1][18+j] += QSUPGy * WJ;


            double Hx = dphi_dx[0][i] * phi_[j] * tPSPG_ * alpha_m;
            double Hy = dphi_dx[1][i] * phi_[j] * tPSPG_ * alpha_m;
            
            double Gx = dphi_dx[0][i] * wSUPGj * tPSPG_ * alpha_f * gamma * dTime_;
            double Gy = dphi_dx[1][i] * wSUPGj * tPSPG_ * alpha_f * gamma * dTime_;

            double Q = (dphi_dx[0][i] * dphi_dx[0][j] + 
                        dphi_dx[1][i] * dphi_dx[1][j]) * tPSPG_ / (dens_);
           
            jacobianNRMatrix[18+i][2*j  ] += (Hx + Gx) * WJ;
            jacobianNRMatrix[18+i][2*j+1] += (Hy + Gy)* WJ;
            jacobianNRMatrix[18+i][18+j] +=  Q * WJ;

        };
    };

    //Stokes Problem   
   //  for (int i = 0; i < 9; i++){
   //      for (int j = 0; j < 9; j++){
                   
   //          //Difusion matrix (viscosity)
   //          double Kxx = (2. * dphi_dx[0][i] * dphi_dx[0][j] + 
   //                       dphi_dx[1][i] * dphi_dx[1][j]) * visc_* alpha_f * gamma * dTime_;
   //          double Kxy = dphi_dx[1][i] * dphi_dx[0][j] * visc_ * alpha_f * gamma * dTime_;
   //          double Kyx = dphi_dx[0][i] * dphi_dx[1][j] * visc_ * alpha_f * gamma * dTime_;
   //          double Kyy = (2. * dphi_dx[1][i] * dphi_dx[1][j] + 
   //                       dphi_dx[0][i] * dphi_dx[0][j]) * visc_* alpha_f * gamma * dTime_;
            
   //          jacobianNRMatrix[2*i  ][2*j  ] += Kxx * weight_ * djac_ * wna_;
   //          jacobianNRMatrix[2*i+1][2*j+1] += Kyy * weight_ * djac_* wna_;
   //          jacobianNRMatrix[2*i  ][2*j+1] += (Kxy) * WJ;
   //          jacobianNRMatrix[2*i+1][2*j  ] += (Kyx) * WJ;

			// //multiply pressure 
   //          double QSUPGx = - dphi_dx[0][i] * phi_[j];
   //          double QSUPGy = - dphi_dx[1][i] * phi_[j];

   //          //multiply velocity 
   //          double Qx = dphi_dx[0][j] * phi_[i] * alpha_f * gamma * dTime_;
   //          double Qy = dphi_dx[1][j] * phi_[i] * alpha_f * gamma * dTime_;

   //          jacobianNRMatrix[18+i][2*j  ] += Qx * weight_ * djac_* wna_;
   //          jacobianNRMatrix[18+i][2*j+1] += Qy * weight_ * djac_* wna_;
   //          jacobianNRMatrix[2*i  ][18+j] += QSUPGx * weight_ * djac_* wna_;
   //          jacobianNRMatrix[2*i+1][18+j] += QSUPGy * weight_ * djac_* wna_;

   //          double Q = (dphi_dx[0][i] * dphi_dx[0][j] + 
   //                      dphi_dx[1][i] * dphi_dx[1][j]) * tPSPG_ / (dens_);
           
   //          jacobianNRMatrix[18+i][18+j] +=  Q * weight_ * djac_* wna_;

   //      };
   //  };

    return;
};


template<>
void Element<2>::getMatrixAndVectorsSameMesh_FEM(double *phi_, double **dphi_dx, double **lagrMultMatrix, 
                                                 double *rhsVector1, double *rhsVector2){


    //Fluid Data
    double &alpha_f = parameters->getAlphaF();
    double &dens_ = parameters->getDensity();
    double &k1 = parameters->getArlequinK1();
    double &k2 = parameters->getArlequinK2();

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
    double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * umeshprev_);

    double duna_dx = (alpha_f * du_dx + (1. - alpha_f) * duprev_dx);
    double duna_dy = (alpha_f * du_dy + (1. - alpha_f) * duprev_dy);
    double dvna_dx = (alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx);
    double dvna_dy = (alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy);

    double WJ = weight_ * djac_;

    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){

            // //L2 operator
            double L2 = phi_[i] * phi_[j] * k1;

            lagrMultMatrix[2*i][2*j] -= L2 *  WJ ;
            lagrMultMatrix[2*i+1][2*j+1] -= L2 *  WJ ;

            //H1 operator
            double H1xx = (2. * dphi_dx[0][i] * dphi_dx[0][j] + 
                          dphi_dx[1][i] * dphi_dx[1][j]) * k2;
            double H1xy = dphi_dx[1][i] * dphi_dx[0][j]* k2;
            double H1yx = dphi_dx[0][i] * dphi_dx[1][j]* k2;
            double H1yy = (2. * dphi_dx[1][i] * dphi_dx[1][j] + 
                           dphi_dx[0][i] * dphi_dx[0][j])* k2;

            lagrMultMatrix[2*i][2*j] -= H1xx * WJ;
            lagrMultMatrix[2*i+1][2*j] -= H1yx  * WJ;
            lagrMultMatrix[2*i][2*j+1] -= H1xy  * WJ;
            lagrMultMatrix[2*i+1][2*j+1] -= H1yy  * WJ;

        };

        //L2 operator - Lagrange
        double l2x_ = phi_[i] * lamx_ * k1;
        double l2y_ = phi_[i] * lamy_ * k1;

        //H1 operator - Lagrange
        double h1x_ = (2. * dphi_dx[0][i] * lamx_dx + 
                      dphi_dx[1][i] * lamx_dy + 
                      dphi_dx[1][i] * lamy_dx) * k2 ;
        
        double h1y_ = (dphi_dx[0][i] * lamx_dy + 
                      2. * dphi_dx[1][i] * lamy_dy + 
                      dphi_dx[0][i] * lamy_dx) * k2;  

        rhsVector1[2*i] += (l2x_+ h1x_) * WJ;
        rhsVector1[2*i+1] += (l2y_+ h1y_) * WJ;

        //L2 operator - Velocity
        double l2ux_ = phi_[i] * una_ * k1;
        double l2uy_ = phi_[i] * vna_ * k1;

        //H1 operator - Velocity
        double h1ux_ = (2. * dphi_dx[0][i] * duna_dx+ 
                       dphi_dx[1][i] * duna_dy + 
                       dphi_dx[1][i] * dvna_dx) * k2;
        
        double h1uy_= (dphi_dx[0][i] * duna_dy + 
                      2. * dphi_dx[1][i] * dvna_dy + 
                      dphi_dx[0][i] * dvna_dx) * k2; 

        rhsVector2[2*i] += (l2ux_ + h1ux_) * WJ;
        rhsVector2[2*i+1] += (l2uy_ + h1uy_)  * WJ;

    };

};


template<>
void Element<2>::getMatrixAndVectorsSameMesh_tSUPG_tPSPG_FEM(double *phi_, double **dphi_dx, double **jacobianNRMatrix,double *rhsVector){

	double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
	double &k1 = parameters->getArlequinK1();

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
	double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * umeshprev_);


	double WJ = weight_ * djac_;

	for (int i = 0; i < 6; i++){
		for (int j = 0; j < 6; j++){

			double msupg = -(una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) * phi_[j] * tSUPG_;
			jacobianNRMatrix[2*i][2*j] += msupg * WJ * k1;
			jacobianNRMatrix[2*i+1][2*j+1] += msupg * WJ * k1;

			double mpspgx = -(dphi_dx[0][i] * phi_[j]) * tPSPG_ / dens_;
			double mpspgy = -(dphi_dx[1][i] * phi_[j]) * tPSPG_ / dens_;

			jacobianNRMatrix[12+i][2*j] += mpspgx * WJ * k1;
			jacobianNRMatrix[12+i][2*j+1] += mpspgy * WJ * k1;
		};

		double vsupgx = (una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) * lamx_ * tSUPG_;
		double vsupgy = (una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) * lamy_ * tSUPG_;

		double vpspg = (dphi_dx[0][i] * lamx_ + dphi_dx[1][i] * lamy_) * tPSPG_ / dens_;

		rhsVector[2*i] += vsupgx * WJ * k1;
		rhsVector[2*i+1] += vsupgy * WJ * k1;
		rhsVector[12+i] += vpspg * WJ * k1;

	};
};


template<>
void Element<2>::getMatrixAndVectorsSameMeshArlqStab_FEM(int &index, double *phi_, double **dphi_dx, double ***ddphi_dx,
                                                         double **arlequinStabD, double *arlequinStabVectorD,
                                                         double **arlequinStab1, double *arlequinStabVector1){

    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &gamma = parameters->getGamma();
    double &dTime_ = parameters->getTimeStep();

    double daxm_dx = alpha_m * dax_dx + (1. - alpha_m) * daxprev_dx;
    double daxm_dy = alpha_m * dax_dy + (1. - alpha_m) * daxprev_dy;
    double daym_dx = alpha_m * day_dx + (1. - alpha_m) * dayprev_dx;
    double daym_dy = alpha_m * day_dy + (1. - alpha_m) * dayprev_dy;

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
    double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * umeshprev_);

    double duna_dx = (alpha_f * du_dx + (1. - alpha_f) * duprev_dx);
    double duna_dy = (alpha_f * du_dy + (1. - alpha_f) * duprev_dy);
    double dvna_dx = (alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx);
    double dvna_dy = (alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy);

    double dumeshna_dx = (alpha_f * dumesh_dx + (1. - alpha_f) * dumeshprev_dx);
    double dumeshna_dy = (alpha_f * dumesh_dy + (1. - alpha_f) * dumeshprev_dy);
    double dvmeshna_dx = (alpha_f * dvmesh_dx + (1. - alpha_f) * dvmeshprev_dx);
    double dvmeshna_dy = (alpha_f * dvmesh_dy + (1. - alpha_f) * dvmeshprev_dy);

    double dduna_dxdx = (alpha_f * ddu_dxdx + (1. - alpha_f) * dduprev_dxdx);
    double dduna_dxdy = (alpha_f * ddu_dxdy + (1. - alpha_f) * dduprev_dxdy);
    double dduna_dydx = (alpha_f * ddu_dydx + (1. - alpha_f) * dduprev_dydx);
    double dduna_dydy = (alpha_f * ddu_dydy + (1. - alpha_f) * dduprev_dydy);
    double ddvna_dxdx = (alpha_f * ddv_dxdx + (1. - alpha_f) * ddvprev_dxdx);
    double ddvna_dxdy = (alpha_f * ddv_dxdy + (1. - alpha_f) * ddvprev_dxdy);
    double ddvna_dydx = (alpha_f * ddv_dydx + (1. - alpha_f) * ddvprev_dydx);
    double ddvna_dydy = (alpha_f * ddv_dydy + (1. - alpha_f) * ddvprev_dydy);

    double wna_ = alpha_f * intPointWeightFunction_FEM[index] + (1. - alpha_f) * intPointWeightFunctionPrev_FEM[index];
    double WJ = weight_ * djac_;

    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){

            //tarlq x mass matrix
            double mass = tARLQ_ * wna_ * (dphi_dx[0][i]* dphi_dx[0][j] + dphi_dx[1][i]* dphi_dx[1][j]);

            arlequinStab1[2*i][2*j] -= mass * WJ * alpha_m;
            arlequinStab1[2*i+1][2*j+1] -= mass * WJ * alpha_m;

            // tarlq x convecction matrix
            double convec1 = tARLQ_ * wna_ * ((dphi_dx[0][i] * (ddphi_dx[0][0][j] * una_ + ddphi_dx[0][1][j] * vna_)) + 
                                              (dphi_dx[1][i] * (ddphi_dx[1][0][j] * una_ + ddphi_dx[1][1][j] * vna_)));
            double convec2 = tARLQ_ * wna_ * ((dphi_dx[0][i] * (dphi_dx[0][j] * (duna_dx - dumeshna_dx) + dphi_dx[1][j] * (dvna_dx - dvmeshna_dx))) + 
                                              (dphi_dx[1][i] * (dphi_dx[0][j] * (duna_dy - dumeshna_dy)  + dphi_dx[1][j] * (dvna_dy - dvmeshna_dy))));
            double convec00 = tARLQ_ * wna_ * (dphi_dx[0][i] * phi_[j] * duna_dx);
            double convec01 = tARLQ_ * wna_ * (dphi_dx[0][i] * phi_[j] * duna_dy);
            double convec10 = tARLQ_ * wna_ * (dphi_dx[1][i] * phi_[j] * dvna_dx);
            double convec11 = tARLQ_ * wna_ * (dphi_dx[1][i] * phi_[j] * dvna_dy);

            arlequinStab1[2*i][2*j] -= (convec1 + convec2 + convec00) * WJ * alpha_f * gamma * dTime_;
            arlequinStab1[2*i][2*j+1] -= (convec01) * WJ * alpha_f * gamma * dTime_;
        	arlequinStab1[2*i+1][2*j] -= (convec10) * WJ * alpha_f * gamma * dTime_;
        	arlequinStab1[2*i+1][2*j+1] -= (convec1 + convec2 + convec11) * WJ * alpha_f * gamma * dTime_;

            // // //tarlq x pressure term
            // double pressx = tARLQ_/dens_ * wna_ * (dphi_dx[0][i] * ddphi_dx[0][0][j] + dphi_dx[1][i] * ddphi_dx[0][1][j]);
            // double pressy = tARLQ_/dens_ * wna_ * (dphi_dx[0][i] * ddphi_dx[1][0][j] + dphi_dx[1][i] * ddphi_dx[1][1][j]);
            // arlequinStab1[2*i][12+j] -= pressx * WJ;
            // arlequinStab1[2*i+1][12+j] -= pressy * WJ;

            //Diagonal matrix
            double LL = (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * tARLQ_ / dens_;

            arlequinStabD[2*i][2*j] += LL * WJ;
            arlequinStabD[2*i+1][2*j+1] += LL * WJ; 
        };


        //tarlq x mass matrix
        double massx = tARLQ_ * wna_ * (dphi_dx[0][i] * daxm_dx + dphi_dx[1][i] * daxm_dy);
        double massy = tARLQ_ * wna_ * (dphi_dx[0][i] * daym_dx + dphi_dx[1][i] * daym_dy);

        arlequinStabVector1[2*i] += massx * WJ;
        arlequinStabVector1[2*i+1] += massy * WJ;

        // tarlq x convecction matrix
        double convec1x = tARLQ_ * wna_ * ((dphi_dx[0][i] * (dduna_dxdx * una_ + dduna_dxdy * vna_)) + 
                                           (dphi_dx[1][i] * (dduna_dydx * una_ + dduna_dydy * vna_)));
        double convec2x = tARLQ_ * wna_ * ((dphi_dx[0][i] * (duna_dx * (duna_dx - dumeshna_dx) + duna_dy * (dvna_dx - dvmeshna_dx))) + 
                                           (dphi_dx[1][i] * (duna_dx * (duna_dy - dumeshna_dy)  + duna_dy * (dvna_dy - dvmeshna_dy))));
        double convec0j = tARLQ_ * wna_ * (dphi_dx[0][i] * una_ * duna_dx + dphi_dx[0][i] * vna_ * duna_dy);
        double convec1y = tARLQ_ * wna_ * ((dphi_dx[0][i] * (ddvna_dxdx * una_ + ddvna_dxdy * vna_)) + 
                                           (dphi_dx[1][i] * (ddvna_dydx * una_ + ddvna_dydy * vna_)));
        double convec2y = tARLQ_ * wna_ * ((dphi_dx[0][i] * (dvna_dx * (duna_dx - dumeshna_dx) + dvna_dy * (dvna_dx - dvmeshna_dx))) + 
                                           (dphi_dx[1][i] * (dvna_dx * (duna_dy - dumeshna_dy)  + dvna_dy * (dvna_dy - dvmeshna_dy))));
        double convec1j = tARLQ_ * wna_ * (dphi_dx[1][i] * una_ * dvna_dx + dphi_dx[1][i] * vna_ * dvna_dy);


        arlequinStabVector1[2*i] += (convec1x + convec2x + convec0j) * WJ;
        arlequinStabVector1[2*i+1] += (convec1y + convec2y + convec1j) * WJ;

        //tarlq x pressure term
        // double pressx = tARLQ_/dens_ * wna_ * (dphi_dx[0][i] * ddp_dxdx + dphi_dx[1][i] * ddp_dydx);
        // double pressy = tARLQ_/dens_ * wna_ * (dphi_dx[0][i] * ddp_dxdy + dphi_dx[1][i] * ddp_dydy);
        
        // arlequinStabVector1[2*i] += pressx * WJ;
        // arlequinStabVector1[2*i+1] += pressy * WJ;

        //Diagonal matrix
        double LLx = (dphi_dx[0][i] * lamx_dx + dphi_dx[1][i] * lamx_dy) * tARLQ_ / dens_;
        double LLy = (dphi_dx[0][i] * lamy_dx + dphi_dx[1][i] * lamy_dy) * tARLQ_ / dens_;

        arlequinStabVectorD[2*i] -= LLx * WJ;
        arlequinStabVectorD[2*i+1] -= LLy * WJ;

    };


};


template<>
void Element<2>::getMatrixAndVectorsSameMesh_ISO(double *phi_, double **dphi_dx, double **lagrMultMatrix, 
        	                                     double *rhsVector1, double *rhsVector2){

	//Fluid Data
    double &alpha_f = parameters->getAlphaF();
    double &dens_ = parameters->getDensity();
    double &k1 = parameters->getArlequinK1();
    double &k2 = parameters->getArlequinK2();

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
	double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * vmeshprev_);

	double duna_dx = (alpha_f * du_dx + (1. - alpha_f) * duprev_dx);
	double duna_dy = (alpha_f * du_dy + (1. - alpha_f) * duprev_dy);
	double dvna_dx = (alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx);
	double dvna_dy = (alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy);
    
    double WJ = weight_ * djac_;

	for (int i = 0; i < 9; i++){
		for (int j = 0; j < 9; j++){

			//L2 operator
			double L2 = phi_[i] * phi_[j] * k1;

			lagrMultMatrix[2*i][2*j] -= L2 * WJ;
			lagrMultMatrix[2*i+1][2*j+1] -= L2 * WJ;

			//H1 operator
			double H1xx = (2. * dphi_dx[0][i] * dphi_dx[0][j] + 
                          dphi_dx[1][i] * dphi_dx[1][j]) * k2;
			double H1xy = dphi_dx[1][i] * dphi_dx[0][j]* k2;
			double H1yx = dphi_dx[0][i] * dphi_dx[1][j]* k2;
			double H1yy = (2. * dphi_dx[1][i] * dphi_dx[1][j] + 
                           dphi_dx[0][i] * dphi_dx[0][j])* k2;

			lagrMultMatrix[2*i][2*j] -= H1xx * WJ;
			lagrMultMatrix[2*i+1][2*j] -= H1yx  * WJ;
			lagrMultMatrix[2*i][2*j+1] -= H1xy  * WJ;
			lagrMultMatrix[2*i+1][2*j+1] -= H1yy  * WJ;

		};

		//L2 operator - Lagrange
		double l2x_ = phi_[i] * lamx_ * k1;
		double l2y_ = phi_[i] * lamy_ * k1;

		//H1 operator - Lagrange
		double h1x_ = (2. * dphi_dx[0][i] * lamx_dx + 
                      dphi_dx[1][i] * lamx_dy + 
                      dphi_dx[1][i] * lamy_dx) * k2 ;
        
        double h1y_ = (dphi_dx[0][i] * lamx_dy + 
                      2. * dphi_dx[1][i] * lamy_dy + 
                      dphi_dx[0][i] * lamy_dx) * k2;  

		rhsVector1[2*i] += (l2x_+ h1x_) * WJ;
		rhsVector1[2*i+1] += (l2y_+ h1y_) * WJ;

		//L2 operator - Velocity
		double l2ux_ = phi_[i] * una_ * k1;
		double l2uy_ = phi_[i] * vna_ * k1;

		//H1 operator - Velocity
		double h1ux_ = (2. * dphi_dx[0][i] * duna_dx+ 
                       dphi_dx[1][i] * duna_dy + 
                       dphi_dx[1][i] * dvna_dx) * k2;
        
        double h1uy_= (dphi_dx[0][i] * duna_dy + 
                      2. * dphi_dx[1][i] * dvna_dy + 
                      dphi_dx[0][i] * dvna_dx) * k2; 

		rhsVector2[2*i] += (l2ux_ + h1ux_) * WJ;
		rhsVector2[2*i+1] += (l2uy_ + h1uy_)  * WJ;

	};

};


template<>
void Element<2>::getMatrixAndVectorsSameMesh_tSUPG_tPSPG_ISO(double *phi_, double **dphi_dx, double **jacobianNRMatrix,double *rhsVector){

    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &k1 = parameters->getArlequinK1();

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
    double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * umeshprev_);

    double WJ = weight_ * djac_;

    for (int i = 0; i < 9; i++){
        for (int j = 0; j < 9; j++){

            double msupg = -(una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) * phi_[j] * tSUPG_;
            jacobianNRMatrix[2*i][2*j] += msupg * WJ * k1;
            jacobianNRMatrix[2*i+1][2*j+1] += msupg * WJ * k1;

            double mpspgx = -(dphi_dx[0][i] * phi_[j]) * tPSPG_ / dens_;
            double mpspgy = -(dphi_dx[1][i] * phi_[j]) * tPSPG_ / dens_;

            jacobianNRMatrix[18+i][2*j] += mpspgx * WJ * k1;
            jacobianNRMatrix[18+i][2*j+1] += mpspgy * WJ * k1;
        };

        double vsupgx = (una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) * lamx_ * tSUPG_;
        double vsupgy = (una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) * lamy_ * tSUPG_;

        double vpspg = (dphi_dx[0][i] * lamx_ + dphi_dx[1][i] * lamy_) * tPSPG_ / dens_;

        rhsVector[2*i] += vsupgx * WJ * k1;
        rhsVector[2*i+1] += vsupgy * WJ * k1;
        rhsVector[18+i] += vpspg * WJ * k1;

    };
};


template<>
void Element<2>::getMatrixAndVectorsSameMeshArlqStab_ISO(int &index, double *phi_, double **dphi_dx, double ***ddphi_dx,
                                                         double **arlequinStabD, double *arlequinStabVectorD,
                                                         double **arlequinStab1, double *arlequinStabVector1){

    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &gamma = parameters->getGamma();
    double &dTime_ = parameters->getTimeStep();

    double daxm_dx = alpha_m * dax_dx + (1. - alpha_m) * daxprev_dx;
    double daxm_dy = alpha_m * dax_dy + (1. - alpha_m) * daxprev_dy;
    double daym_dx = alpha_m * day_dx + (1. - alpha_m) * dayprev_dx;
    double daym_dy = alpha_m * day_dy + (1. - alpha_m) * dayprev_dy;

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
    double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * umeshprev_);

    double duna_dx = (alpha_f * du_dx + (1. - alpha_f) * duprev_dx);
    double duna_dy = (alpha_f * du_dy + (1. - alpha_f) * duprev_dy);
    double dvna_dx = (alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx);
    double dvna_dy = (alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy);


    double dumeshna_dx = (alpha_f * dumesh_dx + (1. - alpha_f) * dumeshprev_dx);
    double dumeshna_dy = (alpha_f * dumesh_dy + (1. - alpha_f) * dumeshprev_dy);
    double dvmeshna_dx = (alpha_f * dvmesh_dx + (1. - alpha_f) * dvmeshprev_dx);
    double dvmeshna_dy = (alpha_f * dvmesh_dy + (1. - alpha_f) * dvmeshprev_dy);

    double dduna_dxdx = (alpha_f * ddu_dxdx + (1. - alpha_f) * dduprev_dxdx);
    double dduna_dxdy = (alpha_f * ddu_dxdy + (1. - alpha_f) * dduprev_dxdy);
    double dduna_dydx = (alpha_f * ddu_dydx + (1. - alpha_f) * dduprev_dydx);
    double dduna_dydy = (alpha_f * ddu_dydy + (1. - alpha_f) * dduprev_dydy);
    double ddvna_dxdx = (alpha_f * ddv_dxdx + (1. - alpha_f) * ddvprev_dxdx);
    double ddvna_dxdy = (alpha_f * ddv_dxdy + (1. - alpha_f) * ddvprev_dxdy);
    double ddvna_dydx = (alpha_f * ddv_dydx + (1. - alpha_f) * ddvprev_dydx);
    double ddvna_dydy = (alpha_f * ddv_dydy + (1. - alpha_f) * ddvprev_dydy);

    double wna_ = alpha_f * intPointWeightFunction_ISO[index] + (1. - alpha_f) * intPointWeightFunctionPrev_ISO[index];
    double WJ = weight_ * djac_;

    for (int i = 0; i < 9; i++){
        for (int j = 0; j < 9; j++){

         //    //tarlq x mass matrix
         //    double mass = tARLQ_ * wna_ * (dphi_dx[0][i]* dphi_dx[0][j] + dphi_dx[1][i]* dphi_dx[1][j]);

         //    arlequinStab1[2*i][2*j] -= mass * WJ * alpha_m;
         //    arlequinStab1[2*i+1][2*j+1] -= mass * WJ * alpha_m;

         //    // tarlq x convecction matrix
            // double convec1 = tARLQ_ * wna_ * ((dphi_dx[0][i] * (ddphi_dx[0][0][j] * una_ + ddphi_dx[0][1][j] * vna_)) + 
            //                                   (dphi_dx[1][i] * (ddphi_dx[1][0][j] * una_ + ddphi_dx[1][1][j] * vna_)));
            // double convec2 = tARLQ_ * wna_ * ((dphi_dx[0][i] * (dphi_dx[0][j] * (duna_dx - dumeshna_dx) + dphi_dx[1][j] * (dvna_dx - dvmeshna_dx))) + 
            //                                   (dphi_dx[1][i] * (dphi_dx[0][j] * (duna_dy - dumeshna_dy)  + dphi_dx[1][j] * (dvna_dy - dvmeshna_dy))));
            // double convec00 = tARLQ_ * wna_ * (dphi_dx[0][i] * phi_[j] * duna_dx);
            // double convec01 = tARLQ_ * wna_ * (dphi_dx[0][i] * phi_[j] * duna_dy);
            // double convec10 = tARLQ_ * wna_ * (dphi_dx[1][i] * phi_[j] * dvna_dx);
            // double convec11 = tARLQ_ * wna_ * (dphi_dx[1][i] * phi_[j] * dvna_dy);
            
         //    arlequinStab1[2*i][2*j] -= (convec1 + convec2 + convec00) * WJ * alpha_f * gamma * dTime_;
         //    arlequinStab1[2*i][2*j+1] -= convec01 * WJ * alpha_f * gamma * dTime_;
        	// arlequinStab1[2*i+1][2*j] -= convec10 * WJ * alpha_f * gamma * dTime_;
        	// arlequinStab1[2*i+1][2*j+1] -= (convec1 + convec2 + convec11) * WJ * alpha_f * gamma * dTime_;

         //    //tarlq x pressure term
         //    double pressx = tARLQ_/dens_ * wna_ * (dphi_dx[0][i] * ddphi_dx[0][0][j] + dphi_dx[1][i] * ddphi_dx[0][1][j]);
         //    double pressy = tARLQ_/dens_ * wna_ * (dphi_dx[0][i] * ddphi_dx[1][0][j] + dphi_dx[1][i] * ddphi_dx[1][1][j]);
         //    arlequinStab1[2*i][18+j] -= pressx * WJ;
         //    arlequinStab1[2*i+1][18+j] -= pressy * WJ;

            //Diagonal matrix
            double LL = (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * tARLQ_ / dens_;

            arlequinStabD[2*i][2*j] += LL * WJ;
            arlequinStabD[2*i+1][2*j+1] += LL * WJ; 
        };

        // //tarlq x mass matrix
        // double massx = tARLQ_ * wna_ * (dphi_dx[0][i]* daxm_dx + dphi_dx[1][i]* daxm_dy);
        // double massy = tARLQ_ * wna_ * (dphi_dx[0][i]* daym_dx + dphi_dx[1][i]* daym_dy);

        // arlequinStabVector1[2*i] += massx * WJ;
        // arlequinStabVector1[2*i+1] += massy * WJ;

        // // tarlq x convecction matrix
        // double convec1x = tARLQ_ * wna_ * ((dphi_dx[0][i] * (dduna_dxdx * una_ + dduna_dxdy * vna_)) + 
        //                                    (dphi_dx[1][i] * (dduna_dydx * una_ + dduna_dydy * vna_)));
        // double convec2x = tARLQ_ * wna_ * ((dphi_dx[0][i] * (duna_dx * (duna_dx - dumeshna_dx) + duna_dy * (dvna_dx - dvmeshna_dx))) + 
        //                                    (dphi_dx[1][i] * (duna_dx * (duna_dy - dumeshna_dy)  + duna_dy * (dvna_dy - dvmeshna_dy))));
        // double convec0j = tARLQ_ * wna_ * (dphi_dx[0][i] * una_ * duna_dx + dphi_dx[0][i] * vna_ * duna_dy);
        // double convec1y = tARLQ_ * wna_ * ((dphi_dx[0][i] * (ddvna_dxdx * una_ + ddvna_dxdy * vna_)) + 
        //                                    (dphi_dx[1][i] * (ddvna_dydx * una_ + ddvna_dydy * vna_)));
        // double convec2y = tARLQ_ * wna_ * ((dphi_dx[0][i] * (dvna_dx * (duna_dx - dumeshna_dx) + dvna_dy * (dvna_dx - dvmeshna_dx))) + 
        //                                    (dphi_dx[1][i] * (dvna_dx * (duna_dy - dumeshna_dy)  + dvna_dy * (dvna_dy - dvmeshna_dy))));
        // double convec1j = tARLQ_ * wna_ * (dphi_dx[1][i] * una_ * dvna_dx + dphi_dx[1][i] * vna_ * dvna_dy);


        // arlequinStabVector1[2*i] += (convec1x + convec2x + convec0j) * WJ;
        // arlequinStabVector1[2*i+1] += (convec1y + convec2y + convec1j) * WJ;

        // // tarlq x pressure term
        // double pressx = tARLQ_/dens_ * wna_ * (dphi_dx[0][i] * ddp_dxdx + dphi_dx[1][i] * ddp_dydx);
        // double pressy = tARLQ_/dens_ * wna_ * (dphi_dx[0][i] * ddp_dxdy + dphi_dx[1][i] * ddp_dydy);
        
        // arlequinStabVector1[2*i] += pressx * WJ;
        // arlequinStabVector1[2*i+1] += pressy * WJ;

        //Diagonal matrix
        double LLx = (dphi_dx[0][i] * lamx_dx + dphi_dx[1][i] * lamx_dy) * tARLQ_ / dens_;
        double LLy = (dphi_dx[0][i] * lamy_dx + dphi_dx[1][i] * lamy_dy) * tARLQ_ / dens_;

        arlequinStabVectorD[2*i] -= LLx * WJ;
        arlequinStabVectorD[2*i+1] -= LLy * WJ;

    };


};


template<>
void Element<2>::getMatrixAndVectorsDifferentMesh_FEM_FEM(double *phi_, double **dphi_dx, double *phiC_, double **dphiC_dx,
                                                          double **lagrMultMatrix, double *rhsVector1, double *rhsVector2,
                                                          double **arlequinStab, double *arlequinStabVector){

    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &k1 = parameters->getArlequinK1();
    double &k2 = parameters->getArlequinK2();

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
    double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * vmeshprev_);

    double duna_dx = (alpha_f * du_dx + (1. - alpha_f) * duprev_dx);
    double duna_dy = (alpha_f * du_dy + (1. - alpha_f) * duprev_dy);
    double dvna_dx = (alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx);
    double dvna_dy = (alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy);
    
    double WJ = weight_ * djac_;
    
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){

            //L2 operator
            double L2 = phi_[i] * phiC_[j] * k1;

            lagrMultMatrix[2*i][2*j] += L2 * WJ;
            lagrMultMatrix[2*i+1][2*j+1] += L2 * WJ;

            //H1 operator
            double H1xx = (2. * dphi_dx[0][i] * dphiC_dx[0][j] + 
                           dphi_dx[1][i] * dphiC_dx[1][j]) * k2;
            double H1xy = dphi_dx[1][i] * dphiC_dx[0][j]* k2;
            double H1yx = dphi_dx[0][i] * dphiC_dx[1][j]* k2;
            double H1yy = (2. * dphi_dx[1][i] * dphiC_dx[1][j] + 
                           dphi_dx[0][i] * dphiC_dx[0][j])* k2;

            lagrMultMatrix[2*i][2*j] += H1xx * WJ;
            lagrMultMatrix[2*i+1][2*j] += H1yx * WJ;
            lagrMultMatrix[2*i][2*j+1] += H1xy * WJ;
            lagrMultMatrix[2*i+1][2*j+1] += H1yy * WJ;

             // Arlequin Stabilization term
            double LL = (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * tARLQ_ / dens_;

            arlequinStab[2*i][2*j] += LL * WJ;
            arlequinStab[2*i+1][2*j+1] += LL * WJ; 

        };


        //L2 operator - Lagrange
        double l2x_ = phiC_[i] * lamx_ * k1;
        double l2y_ = phiC_[i] * lamy_ * k1;

        //H1 operator - Lagrange
        double h1x_ = (2. * dphiC_dx[0][i] * lamx_dx + 
                      dphiC_dx[1][i] * lamx_dy + 
                      dphiC_dx[1][i] * lamy_dx) * k2 ;
        
        double h1y_ = (dphiC_dx[0][i] * lamx_dy + 
                      2. * dphiC_dx[1][i] * lamy_dy + 
                      dphiC_dx[0][i] * lamy_dx) * k2;  

        rhsVector1[2*i] -= (l2x_+ h1x_)  * WJ;
        rhsVector1[2*i+1] -= (l2y_+ h1y_)  * WJ;


        //L2 operator - Velocity
        double l2ux_ = phi_[i] * una_ * k1;
        double l2uy_ = phi_[i] * vna_ * k1;

        //H1 operator - Velocity
        double h1ux_ = (2. * dphi_dx[0][i] * duna_dx+ 
                       dphi_dx[1][i] * duna_dy + 
                       dphi_dx[1][i] * dvna_dx) * k2;
        
        double h1uy_= (dphi_dx[0][i] * duna_dy + 
                      2. * dphi_dx[1][i] * dvna_dy + 
                      dphi_dx[0][i] * dvna_dx) * k2; 

        rhsVector2[2*i] -= (l2ux_ + h1ux_) * WJ;
        rhsVector2[2*i+1] -= (l2uy_ + h1uy_) * WJ;

        //Arlequin Stabilization vector
        double LLx = (dphi_dx[0][i] * lamx_dx + dphi_dx[1][i] * lamx_dy) * tARLQ_ / dens_;
        double LLy = (dphi_dx[0][i] * lamy_dx + dphi_dx[1][i] * lamy_dy) * tARLQ_ / dens_;

        arlequinStabVector[2*i] -= LLx * WJ;
        arlequinStabVector[2*i+1] -= LLy * WJ;

    };
};

template<>
void Element<2>::getMatrixAndVectorsDifferentMesh_FEM_ISO(double *phi_, double **dphi_dx, double *phiC_, double **dphiC_dx,
                                                         double **lagrMultMatrix, double *rhsVector1, double *rhsVector2){

    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &k1 = parameters->getArlequinK1();
    double &k2 = parameters->getArlequinK2();

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
    double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * vmeshprev_);

    double duna_dx = (alpha_f * du_dx + (1. - alpha_f) * duprev_dx);
    double duna_dy = (alpha_f * du_dy + (1. - alpha_f) * duprev_dy);
    double dvna_dx = (alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx);
    double dvna_dy = (alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy);
    
    double WJ = weight_ * djac_;
    
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 9; j++){

            //L2 operator
            double L2 = phi_[i] * phiC_[j] * k1;

            lagrMultMatrix[2*i][2*j] += L2 * WJ;
            lagrMultMatrix[2*i+1][2*j+1] += L2 * WJ;

            //H1 operator
            double H1xx = (2. * dphi_dx[0][i] * dphiC_dx[0][j] + 
                           dphi_dx[1][i] * dphiC_dx[1][j]) * k2;
            double H1xy = dphi_dx[1][i] * dphiC_dx[0][j]* k2;
            double H1yx = dphi_dx[0][i] * dphiC_dx[1][j]* k2;
            double H1yy = (2. * dphi_dx[1][i] * dphiC_dx[1][j] + 
                           dphi_dx[0][i] * dphiC_dx[0][j])* k2;

            lagrMultMatrix[2*i][2*j] += H1xx * WJ;
            lagrMultMatrix[2*i+1][2*j] += H1yx * WJ;
            lagrMultMatrix[2*i][2*j+1] += H1xy * WJ;
            lagrMultMatrix[2*i+1][2*j+1] += H1yy * WJ;

        };

        //L2 operator - Velocity
        double l2ux_ = phi_[i] * una_ * k1;
        double l2uy_ = phi_[i] * vna_ * k1;

        //H1 operator - Velocity
        double h1ux_ = (2. * dphi_dx[0][i] * duna_dx+ 
                       dphi_dx[1][i] * duna_dy + 
                       dphi_dx[1][i] * dvna_dx) * k2;
        
        double h1uy_= (dphi_dx[0][i] * duna_dy + 
                      2. * dphi_dx[1][i] * dvna_dy + 
                      dphi_dx[0][i] * dvna_dx) * k2; 

        rhsVector2[2*i] -= (l2ux_ + h1ux_) * WJ;
        rhsVector2[2*i+1] -= (l2uy_ + h1uy_) * WJ;

    };


    for (int i = 0; i < 9; i ++){

        //L2 operator - Lagrange
        double l2x_ = phiC_[i] * lamx_ * k1;
        double l2y_ = phiC_[i] * lamy_ * k1;

        //H1 operator - Lagrange
        double h1x_ = (2. * dphiC_dx[0][i] * lamx_dx + 
                      dphiC_dx[1][i] * lamx_dy + 
                      dphiC_dx[1][i] * lamy_dx) * k2 ;
        
        double h1y_ = (dphiC_dx[0][i] * lamx_dy + 
                      2. * dphiC_dx[1][i] * lamy_dy + 
                      dphiC_dx[0][i] * lamy_dx) * k2;  

        rhsVector1[2*i] -= (l2x_+ h1x_)  * WJ;
        rhsVector1[2*i+1] -= (l2y_+ h1y_)  * WJ;

    };

};

template<>
void Element<2>::getMatrixAndVectorsDifferentMesh_ISO_FEM(double *phi_, double **dphi_dx, double *phiC_, double **dphiC_dx,
                                                          double **lagrMultMatrix, double *rhsVector1, double *rhsVector2){

    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &k1 = parameters->getArlequinK1();
    double &k2 = parameters->getArlequinK2();

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
    double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * vmeshprev_);

    double duna_dx = (alpha_f * du_dx + (1. - alpha_f) * duprev_dx);
    double duna_dy = (alpha_f * du_dy + (1. - alpha_f) * duprev_dy);
    double dvna_dx = (alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx);
    double dvna_dy = (alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy);
    
    double WJ = weight_ * djac_;
    
    for (int i = 0; i < 9; i++){
        for (int j = 0; j < 6; j++){

            //L2 operator
            double L2 = phi_[i] * phiC_[j] * k1;

            lagrMultMatrix[2*i][2*j] += L2 * WJ;
            lagrMultMatrix[2*i+1][2*j+1] += L2 * WJ;

            //H1 operator
            double H1xx = (2. * dphi_dx[0][i] * dphiC_dx[0][j] + 
                           dphi_dx[1][i] * dphiC_dx[1][j]) * k2;
            double H1xy = dphi_dx[1][i] * dphiC_dx[0][j]* k2;
            double H1yx = dphi_dx[0][i] * dphiC_dx[1][j]* k2;
            double H1yy = (2. * dphi_dx[1][i] * dphiC_dx[1][j] + 
                           dphi_dx[0][i] * dphiC_dx[0][j])* k2;

            lagrMultMatrix[2*i][2*j] += H1xx * WJ;
            lagrMultMatrix[2*i+1][2*j] += H1yx * WJ;
            lagrMultMatrix[2*i][2*j+1] += H1xy * WJ;
            lagrMultMatrix[2*i+1][2*j+1] += H1yy * WJ;

        };

        //L2 operator - Velocity
        double l2ux_ = phi_[i] * una_ * k1;
        double l2uy_ = phi_[i] * vna_ * k1;

        //H1 operator - Velocity
        double h1ux_ = (2. * dphi_dx[0][i] * duna_dx+ 
                       dphi_dx[1][i] * duna_dy + 
                       dphi_dx[1][i] * dvna_dx) * k2;
        
        double h1uy_= (dphi_dx[0][i] * duna_dy + 
                      2. * dphi_dx[1][i] * dvna_dy + 
                      dphi_dx[0][i] * dvna_dx) * k2; 

        rhsVector2[2*i] -= (l2ux_ + h1ux_) * WJ;
        rhsVector2[2*i+1] -= (l2uy_ + h1uy_) * WJ;

    };


    for (int i = 0; i < 6; i ++){

        //L2 operator - Lagrange
        double l2x_ = phiC_[i] * lamx_ * k1;
        double l2y_ = phiC_[i] * lamy_ * k1;

        //H1 operator - Lagrange
        double h1x_ = (2. * dphiC_dx[0][i] * lamx_dx + 
                      dphiC_dx[1][i] * lamx_dy + 
                      dphiC_dx[1][i] * lamy_dx) * k2 ;
        
        double h1y_ = (dphiC_dx[0][i] * lamx_dy + 
                      2. * dphiC_dx[1][i] * lamy_dy + 
                      dphiC_dx[0][i] * lamy_dx) * k2;  

        rhsVector1[2*i] -= (l2x_+ h1x_)  * WJ;
        rhsVector1[2*i+1] -= (l2y_+ h1y_)  * WJ;

    };

};

template<>
void Element<2>::getMatrixAndVectorsDifferentMesh_tSUPG_tPSPG_FEM_ISO(double *phi_,double **dphiC_dx,double **jacobianNRMatrix,
																	  double *rhsVector){

	double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
	double &k1 = parameters->getArlequinK1();

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
	double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * vmeshprev_);

	double WJ = weight_ * djac_;

	for (int i = 0; i < 9; i++){
		for (int j = 0; j < 6; j++){

			double msupg = -(una_ * dphiC_dx[0][i] + vna_ * dphiC_dx[1][i]) * phi_[j] * tSUPG_;
			jacobianNRMatrix[2*i][2*j] += msupg * WJ * k1;
			jacobianNRMatrix[2*i+1][2*j+1] += msupg * WJ * k1;

			double mpspgx = -(dphiC_dx[0][i] * phi_[j]) * tPSPG_ / dens_;
			double mpspgy = -(dphiC_dx[1][i] * phi_[j]) * tPSPG_ / dens_;

			jacobianNRMatrix[18+i][2*j] += mpspgx * WJ * k1;
			jacobianNRMatrix[18+i][2*j+1] += mpspgy * WJ * k1;
		};

		double vsupgx = (una_ * dphiC_dx[0][i] + vna_ * dphiC_dx[1][i]) * lamx_ * tSUPG_;
		double vsupgy = (una_ * dphiC_dx[0][i] + vna_ * dphiC_dx[1][i]) * lamy_ * tSUPG_;

		double vpspg = (dphiC_dx[0][i] * lamx_ + dphiC_dx[1][i] * lamy_) * tPSPG_ / dens_;

		rhsVector[2*i] += vsupgx * WJ * k1;
		rhsVector[2*i+1] += vsupgy * WJ * k1;
		rhsVector[18+i] += vpspg * WJ * k1;

	};

};

template<>
void Element<2>::getMatrixAndVectorsDifferentMesh_tSUPG_tPSPG_ISO_FEM(double *phi_,double **dphiC_dx,double **jacobianNRMatrix,
                                                                      double *rhsVector){

    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &k1 = parameters->getArlequinK1();

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
	double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * vmeshprev_);

    double WJ = weight_ * djac_;

    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 9; j++){

            double msupg = -(una_ * dphiC_dx[0][i] + vna_ * dphiC_dx[1][i]) * phi_[j] * tSUPG_;
            jacobianNRMatrix[2*i][2*j] += msupg * WJ * k1;
            jacobianNRMatrix[2*i+1][2*j+1] += msupg * WJ * k1;

            double mpspgx = -(dphiC_dx[0][i] * phi_[j]) * tPSPG_ / dens_;
            double mpspgy = -(dphiC_dx[1][i] * phi_[j]) * tPSPG_ / dens_;

            jacobianNRMatrix[12+i][2*j] += mpspgx * WJ * k1;
            jacobianNRMatrix[12+i][2*j+1] += mpspgy * WJ * k1;
        };

        double vsupgx = (una_ * dphiC_dx[0][i] + vna_ * dphiC_dx[1][i]) * lamx_ * tSUPG_;
        double vsupgy = (una_ * dphiC_dx[0][i] + vna_ * dphiC_dx[1][i]) * lamy_ * tSUPG_;

        double vpspg = (dphiC_dx[0][i] * lamx_ + dphiC_dx[1][i] * lamy_) * tPSPG_ / dens_;

        rhsVector[2*i] += vsupgx * WJ * k1;
        rhsVector[2*i+1] += vsupgy * WJ * k1;
        rhsVector[12+i] += vpspg * WJ * k1;

    };

};

template<>
void Element<2>::getMatrixAndVectorsDifferentMeshArlqStab_FEM_ISO(int &index, double **dphi_dx, double *phiC_, double **dphiC_dx,double ***ddphiC_dx,
                                                                  double **arlequinStabD, double *arlequinStabVectorD,
                                                                  double **arlequinStab0, double *arlequinStabVector0){
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &gamma = parameters->getGamma();
    double &dTime_ = parameters->getTimeStep();

    double daxm_dx = alpha_m * dax_dx + (1. - alpha_m) * daxprev_dx;
    double daxm_dy = alpha_m * dax_dy + (1. - alpha_m) * daxprev_dy;
    double daym_dx = alpha_m * day_dx + (1. - alpha_m) * dayprev_dx;
    double daym_dy = alpha_m * day_dy + (1. - alpha_m) * dayprev_dy;

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
    double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * umeshprev_);

    double duna_dx = (alpha_f * du_dx + (1. - alpha_f) * duprev_dx);
    double duna_dy = (alpha_f * du_dy + (1. - alpha_f) * duprev_dy);
    double dvna_dx = (alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx);
    double dvna_dy = (alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy);

    double dumeshna_dx = (alpha_f * dumesh_dx + (1. - alpha_f) * dumeshprev_dx);
    double dumeshna_dy = (alpha_f * dumesh_dy + (1. - alpha_f) * dumeshprev_dy);
    double dvmeshna_dx = (alpha_f * dvmesh_dx + (1. - alpha_f) * dvmeshprev_dx);
    double dvmeshna_dy = (alpha_f * dvmesh_dy + (1. - alpha_f) * dvmeshprev_dy);

    double dduna_dxdx = (alpha_f * ddu_dxdx + (1. - alpha_f) * dduprev_dxdx);
    double dduna_dxdy = (alpha_f * ddu_dxdy + (1. - alpha_f) * dduprev_dxdy);
    double dduna_dydx = (alpha_f * ddu_dydx + (1. - alpha_f) * dduprev_dydx);
    double dduna_dydy = (alpha_f * ddu_dydy + (1. - alpha_f) * dduprev_dydy);
    double ddvna_dxdx = (alpha_f * ddv_dxdx + (1. - alpha_f) * ddvprev_dxdx);
    double ddvna_dxdy = (alpha_f * ddv_dxdy + (1. - alpha_f) * ddvprev_dxdy);
    double ddvna_dydx = (alpha_f * ddv_dydx + (1. - alpha_f) * ddvprev_dydx);
    double ddvna_dydy = (alpha_f * ddv_dydy + (1. - alpha_f) * ddvprev_dydy);

    double wna_ = alpha_f * intPointWeightFunction_FEM[index] + (1. - alpha_f) * intPointWeightFunctionPrev_FEM[index];
    double WJ = weight_ * djac_;  

    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 9; j++){

            //tarlq x mass matrix
            double mass = tARLQ_ * wna_ * (dphi_dx[0][i]* dphiC_dx[0][j] + dphi_dx[1][i]* dphiC_dx[1][j]);

            arlequinStab0[2*i][2*j] += mass * WJ * alpha_m;
            arlequinStab0[2*i+1][2*j+1] += mass * WJ * alpha_m;;

            //tarlq x convecction matrix
            double convec1 = tARLQ_ * wna_ * ((dphi_dx[0][i] * (ddphiC_dx[0][0][j] * una_ + ddphiC_dx[0][1][j] * vna_)) + 
                                              (dphi_dx[1][i] * (ddphiC_dx[1][0][j] * una_ + ddphiC_dx[1][1][j] * vna_)));
            double convec2 = tARLQ_ * wna_ * ((dphi_dx[0][i] * (dphiC_dx[0][j] * (duna_dx - dumeshna_dx) + dphiC_dx[1][j] * (dvna_dx - dvmeshna_dx))) + 
                                              (dphi_dx[1][i] * (dphiC_dx[0][j] * (duna_dy - dumeshna_dy)  + dphiC_dx[1][j] * (dvna_dy - dvmeshna_dy))));
    		double convec00 = tARLQ_ * wna_ * (dphi_dx[0][i] * phiC_[j] * duna_dx);
            double convec01 = tARLQ_ * wna_ * (dphi_dx[0][i] * phiC_[j] * duna_dy);
            double convec10 = tARLQ_ * wna_ * (dphi_dx[1][i] * phiC_[j] * dvna_dx);
            double convec11 = tARLQ_ * wna_ * (dphi_dx[1][i] * phiC_[j] * dvna_dy);

            arlequinStab0[2*i][2*j] += (convec1 + convec2 + convec00) * WJ * alpha_f * gamma * dTime_;
            arlequinStab0[2*i][2*j+1] += convec01 * WJ * alpha_f * gamma * dTime_;
        	arlequinStab0[2*i+1][2*j] += convec10 * WJ * alpha_f * gamma * dTime_;
        	arlequinStab0[2*i+1][2*j+1] += (convec1 + convec2 + convec11) * WJ * alpha_f * gamma * dTime_;


            // tarlq x pressure term
            // double pressx = tARLQ_/dens_ * wna_ * (dphi_dx[0][i]* ddphiC_dx[0][0][j] + dphi_dx[1][i]* ddphiC_dx[0][1][j]);
            // double pressy = tARLQ_/dens_ * wna_ * (dphi_dx[0][i]* ddphiC_dx[1][0][j] + dphi_dx[1][i]* ddphiC_dx[1][1][j]);
            // arlequinStab0[2*i][18+j] += pressx * WJ;
            // arlequinStab0[2*i+1][18+j] += pressy * WJ;

        };


        //tarlq x mass matrix
        double massx = tARLQ_ * wna_ * (dphi_dx[0][i]* daxm_dx + dphi_dx[1][i]* daxm_dy);
        double massy = tARLQ_ * wna_ * (dphi_dx[0][i]* daym_dx + dphi_dx[1][i]* daym_dy);

        arlequinStabVector0[2*i] -= massx * WJ;
        arlequinStabVector0[2*i+1] -= massy * WJ;

        //tarlq x convecction matrix
        double convec1x = tARLQ_ * wna_ * ((dphi_dx[0][i] * (dduna_dxdx * una_ + dduna_dxdy * vna_)) + 
                                           (dphi_dx[1][i] * (dduna_dydx * una_ + dduna_dydy * vna_)));
        double convec2x = tARLQ_ * wna_ * ((dphi_dx[0][i] * (duna_dx * (duna_dx - dumeshna_dx) + duna_dy * (dvna_dx - dvmeshna_dx))) + 
                                           (dphi_dx[1][i] * (duna_dx * (duna_dy - dumeshna_dy)  + duna_dy * (dvna_dy - dvmeshna_dy))));
        double convec0j = tARLQ_ * wna_ * (dphi_dx[0][i] * una_ * duna_dx + dphi_dx[0][i] * vna_ * duna_dy);
        double convec1y = tARLQ_ * wna_ * ((dphi_dx[0][i] * (ddvna_dxdx * una_ + ddvna_dxdy * vna_)) + 
                                           (dphi_dx[1][i] * (ddvna_dydx * una_ + ddvna_dydy * vna_)));
        double convec2y = tARLQ_ * wna_ * ((dphi_dx[0][i] * (dvna_dx * (duna_dx - dumeshna_dx) + dvna_dy * (dvna_dx - dvmeshna_dx))) + 
                                           (dphi_dx[1][i] * (dvna_dx * (duna_dy - dumeshna_dy)  + dvna_dy * (dvna_dy - dvmeshna_dy))));
        double convec1j = tARLQ_ * wna_ * (dphi_dx[1][i] * una_ * dvna_dx + dphi_dx[1][i] * vna_ * dvna_dy);

        arlequinStabVector0[2*i] -= (convec1x + convec2x + convec0j) * WJ;
        arlequinStabVector0[2*i+1] -= (convec1y + convec2y + convec1j) * WJ;

        // //tarlq x pressure term
        // double pressx = tARLQ_/dens_ * wna_ * (dphi_dx[0][i]* ddp_dxdx + dphi_dx[1][i]* ddp_dydx);
        // double pressy = tARLQ_/dens_ * wna_ * (dphi_dx[0][i]* ddp_dxdy + dphi_dx[1][i]* ddp_dydy);
        // arlequinStabVector0[2*i] -= pressx * WJ;
        // arlequinStabVector0[2*i+1] -= pressy * WJ;

    };

    //Arlequin Diagonal Stabilization terms
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 6; j++){

            double LL = (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * tARLQ_ / dens_;

            arlequinStabD[2*i][2*j] += LL * WJ;
            arlequinStabD[2*i+1][2*j+1] += LL * WJ; 

        };

        //Stabilization term
        double LLx = (dphi_dx[0][i] * lamx_dx + dphi_dx[1][i] * lamx_dy) * tARLQ_ / dens_;
        double LLy = (dphi_dx[0][i] * lamy_dx + dphi_dx[1][i] * lamy_dy) * tARLQ_ / dens_;

        arlequinStabVectorD[2*i] -= LLx * WJ;
        arlequinStabVectorD[2*i+1] -= LLy * WJ;
    };


};

template<>
void Element<2>::getMatrixAndVectorsDifferentMeshArlqStab_ISO_FEM(int &index, double **dphi_dx, double *phiC_, double **dphiC_dx,double ***ddphiC_dx,
                                                                  double **arlequinStabD, double *arlequinStabVectorD,
                                                                  double **arlequinStab0, double *arlequinStabVector0){
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &gamma = parameters->getGamma();
    double &dTime_ = parameters->getTimeStep();

    double daxm_dx = alpha_m * dax_dx + (1. - alpha_m) * daxprev_dx;
    double daxm_dy = alpha_m * dax_dy + (1. - alpha_m) * daxprev_dy;
    double daym_dx = alpha_m * day_dx + (1. - alpha_m) * dayprev_dx;
    double daym_dy = alpha_m * day_dy + (1. - alpha_m) * dayprev_dy;

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
    double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * umeshprev_);

    double duna_dx = (alpha_f * du_dx + (1. - alpha_f) * duprev_dx);
    double duna_dy = (alpha_f * du_dy + (1. - alpha_f) * duprev_dy);
    double dvna_dx = (alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx);
    double dvna_dy = (alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy);

    double dumeshna_dx = (alpha_f * dumesh_dx + (1. - alpha_f) * dumeshprev_dx);
    double dumeshna_dy = (alpha_f * dumesh_dy + (1. - alpha_f) * dumeshprev_dy);
    double dvmeshna_dx = (alpha_f * dvmesh_dx + (1. - alpha_f) * dvmeshprev_dx);
    double dvmeshna_dy = (alpha_f * dvmesh_dy + (1. - alpha_f) * dvmeshprev_dy);

    double dduna_dxdx = (alpha_f * ddu_dxdx + (1. - alpha_f) * dduprev_dxdx);
    double dduna_dxdy = (alpha_f * ddu_dxdy + (1. - alpha_f) * dduprev_dxdy);
    double dduna_dydx = (alpha_f * ddu_dydx + (1. - alpha_f) * dduprev_dydx);
    double dduna_dydy = (alpha_f * ddu_dydy + (1. - alpha_f) * dduprev_dydy);
    double ddvna_dxdx = (alpha_f * ddv_dxdx + (1. - alpha_f) * ddvprev_dxdx);
    double ddvna_dxdy = (alpha_f * ddv_dxdy + (1. - alpha_f) * ddvprev_dxdy);
    double ddvna_dydx = (alpha_f * ddv_dydx + (1. - alpha_f) * ddvprev_dydx);
    double ddvna_dydy = (alpha_f * ddv_dydy + (1. - alpha_f) * ddvprev_dydy);

    double wna_ = alpha_f * intPointWeightFunction_ISO[index] + (1. - alpha_f) * intPointWeightFunctionPrev_ISO[index];
    double WJ = weight_ * djac_;
   

    // for (int i = 0; i < 9; i++){
    //     for (int j = 0; j < 6; j++){

    //         // tarlq x mass matrix
    //         double mass = tARLQ_ * wna_ * (dphi_dx[0][i]* dphiC_dx[0][j] + dphi_dx[1][i]* dphiC_dx[1][j]);

    //         arlequinStab0[2*i][2*j] += mass * WJ * alpha_m;
    //         arlequinStab0[2*i+1][2*j+1] += mass * WJ * alpha_m;;

    //        //tarlq x convecction matrix
    //         double convec1 = tARLQ_ * wna_ * ((dphi_dx[0][i] * (ddphiC_dx[0][0][j] * una_ + ddphiC_dx[0][1][j] * vna_)) + 
    //                                           (dphi_dx[1][i] * (ddphiC_dx[1][0][j] * una_ + ddphiC_dx[1][1][j] * vna_)));
    // double convec2 = tARLQ_ * wna_ * ((dphi_dx[0][i] * (dphiC_dx[0][j] * (duna_dx - dumeshna_dx) + dphiC_dx[1][j] * (dvna_dx - dvmeshna_dx))) + 
    //                                           (dphi_dx[1][i] * (dphiC_dx[0][j] * (duna_dy - dumeshna_dy)  + dphiC_dx[1][j] * (dvna_dy - dvmeshna_dy))));
    // 		double convec00 = tARLQ_ * wna_ * (dphi_dx[0][i] * phiC_[j] * duna_dx);
    //         double convec01 = tARLQ_ * wna_ * (dphi_dx[0][i] * phiC_[j] * duna_dy);
    //         double convec10 = tARLQ_ * wna_ * (dphi_dx[1][i] * phiC_[j] * dvna_dx);
    //         double convec11 = tARLQ_ * wna_ * (dphi_dx[1][i] * phiC_[j] * dvna_dy);
    //         arlequinStab0[2*i][2*j] += (convec1 + convec2 + convec00) * WJ * alpha_f * gamma * dTime_;
    //         arlequinStab0[2*i][2*j+1] += convec01 * WJ * alpha_f * gamma * dTime_;
    //     	arlequinStab0[2*i+1][2*j] += convec10 * WJ * alpha_f * gamma * dTime_;
    //     	arlequinStab0[2*i+1][2*j+1] += (convec1 + convec2 + convec11) * WJ * alpha_f * gamma * dTime_;

    //         // tarlq x pressure term
    //         double pressx = tARLQ_/dens_ * wna_ * (dphi_dx[0][i]* ddphiC_dx[0][0][j] + dphi_dx[1][i]* ddphiC_dx[0][1][j]);
    //         double pressy = tARLQ_/dens_ * wna_ * (dphi_dx[0][i]* ddphiC_dx[1][0][j] + dphi_dx[1][i]* ddphiC_dx[1][1][j]);
    //         arlequinStab0[2*i][12+j] += pressx * WJ;
    //         arlequinStab0[2*i+1][12+j] += pressy * WJ;

    //     };


    //     // //tarlq x mass matrix
    //     double massx = tARLQ_ * wna_ * (dphi_dx[0][i]* daxm_dx + dphi_dx[1][i]* daxm_dy);
    //     double massy = tARLQ_ * wna_ * (dphi_dx[0][i]* daym_dx + dphi_dx[1][i]* daym_dy);

    //     arlequinStabVector0[2*i] -= massx * WJ;
    //     arlequinStabVector0[2*i+1] -= massy * WJ;

    //     // tarlq x convecction matrix
    //     double convec1x = tARLQ_ * wna_ * ((dphi_dx[0][i] * (dduna_dxdx * una_ + dduna_dxdy * vna_)) + 
    //                                        (dphi_dx[1][i] * (dduna_dydx * una_ + dduna_dydy * vna_)));
// double convec2x = tARLQ_ * wna_ * ((dphi_dx[0][i] * (duna_dx * (duna_dx - dumeshna_dx) + duna_dy * (dvna_dx - dvmeshna_dx))) + 
//                                            (dphi_dx[1][i] * (duna_dx * (duna_dy - dumeshna_dy)  + duna_dy * (dvna_dy - dvmeshna_dy))));
    //     double convec0j = tARLQ_ * wna_ * (dphi_dx[0][i] * una_ * duna_dx + dphi_dx[0][i] * vna_ * duna_dy);
    //     double convec1y = tARLQ_ * wna_ * ((dphi_dx[0][i] * (ddvna_dxdx * una_ + ddvna_dxdy * vna_)) + 
    //                                        (dphi_dx[1][i] * (ddvna_dydx * una_ + ddvna_dydy * vna_)));
// double convec2y = tARLQ_ * wna_ * ((dphi_dx[0][i] * (dvna_dx * (duna_dx - dumeshna_dx) + dvna_dy * (dvna_dx - dvmeshna_dx))) + 
//                                            (dphi_dx[1][i] * (dvna_dx * (duna_dy - dumeshna_dy)  + dvna_dy * (dvna_dy - dvmeshna_dy))));
    //     double convec1j = tARLQ_ * wna_ * (dphi_dx[1][i] * una_ * dvna_dx + dphi_dx[1][i] * vna_ * dvna_dy);


    //     arlequinStabVector0[2*i] -= (convec1x + convec2x + convec0j) * WJ;
    //     arlequinStabVector0[2*i+1] -= (convec1y + convec2y + convec1j) * WJ;

    //     // tarlq x pressure term
    //     double pressx = tARLQ_/dens_ * wna_ * (dphi_dx[0][i]* ddp_dxdx + dphi_dx[1][i]* ddp_dydx);
    //     double pressy = tARLQ_/dens_ * wna_ * (dphi_dx[0][i]* ddp_dxdy + dphi_dx[1][i]* ddp_dydy);
    //     arlequinStabVector0[2*i] -= pressx * WJ;
    //     arlequinStabVector0[2*i+1] -= pressy * WJ;

    // };


    //Arlequin Diagonal Stabilization terms
    for (int i = 0; i < 9; i++){
        for (int j = 0; j < 9; j++){

            double LL = (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * tARLQ_ / dens_;

            arlequinStabD[2*i][2*j] += LL * WJ;
            arlequinStabD[2*i+1][2*j+1] += LL * WJ; 

        };

        //Stabilization term
        double LLx = (dphi_dx[0][i] * lamx_dx + dphi_dx[1][i] * lamx_dy) * tARLQ_ / dens_;
        double LLy = (dphi_dx[0][i] * lamy_dx + dphi_dx[1][i] * lamy_dy) * tARLQ_ / dens_;

        arlequinStabVectorD[2*i] -= LLx * WJ;
        arlequinStabVectorD[2*i+1] -= LLy * WJ;
    };


};

template<>
void Element<2>::getMatrixAndVectorsDifferentMesh_ISO_ISO(double *phi_, double **dphi_dx, double *phiC_, double **dphiC_dx,
        										          double **lagrMultMatrix, double *rhsVector1, double *rhsVector2,
                                                          double **arlequinStab, double *arlequinStabVector){

	double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();
    double &k1 = parameters->getArlequinK1();
    double &k2 = parameters->getArlequinK2();

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
	double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * vmeshprev_) ;

	double duna_dx = (alpha_f * du_dx + (1. - alpha_f) * duprev_dx);
	double duna_dy = (alpha_f * du_dy + (1. - alpha_f) * duprev_dy);
	double dvna_dx = (alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx);
	double dvna_dy = (alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy);
	
    double WJ = weight_ * djac_;
    
	for (int i = 0; i < 9; i++){
		for (int j = 0; j < 9; j++){

			//L2 operator
			double L2 = phi_[i] * phiC_[j] * k1;

			lagrMultMatrix[2*i][2*j] += L2 * WJ;
			lagrMultMatrix[2*i+1][2*j+1] += L2 * WJ;

			//H1 operator
			double H1xx = (2. * dphi_dx[0][i] * dphiC_dx[0][j] + 
                           dphi_dx[1][i] * dphiC_dx[1][j]) * k2;
			double H1xy = dphi_dx[1][i] * dphiC_dx[0][j]* k2;
			double H1yx = dphi_dx[0][i] * dphiC_dx[1][j]* k2;
			double H1yy = (2. * dphi_dx[1][i] * dphiC_dx[1][j] + 
                           dphi_dx[0][i] * dphiC_dx[0][j])* k2;

			lagrMultMatrix[2*i][2*j] += H1xx * WJ;
			lagrMultMatrix[2*i+1][2*j] += H1yx  * WJ;
			lagrMultMatrix[2*i][2*j+1] += H1xy  * WJ;
			lagrMultMatrix[2*i+1][2*j+1] += H1yy  * WJ;

            //Stabilization terms
            double LL = (dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) * tARLQ_ / dens_;

            arlequinStab[2*i][2*j] += LL * WJ;
            arlequinStab[2*i+1][2*j+1] += LL * WJ; 

		};

		//L2 operator - Lagrange
		double l2x_ = phiC_[i] * lamx_ * k1;
		double l2y_ = phiC_[i] * lamy_ * k1;

		//H1 operator - Lagrange
		double h1x_ = (2. * dphiC_dx[0][i] * lamx_dx + 
                      dphiC_dx[1][i] * lamx_dy + 
                      dphiC_dx[1][i] * lamy_dx) * k2 ;
        
        double h1y_ = (dphiC_dx[0][i] * lamx_dy + 
                      2. * dphiC_dx[1][i] * lamy_dy + 
                      dphiC_dx[0][i] * lamy_dx) * k2;  

		rhsVector1[2*i] -= (l2x_+ h1x_)  * WJ;
		rhsVector1[2*i+1] -= (l2y_+ h1y_)  * WJ;

		//L2 operator - Velocity
		double l2ux_ = phi_[i] * una_ * k1;
		double l2uy_ = phi_[i] * vna_ * k1;

		//H1 operator - Velocity
		double h1ux_ = (2. * dphi_dx[0][i] * duna_dx+ 
                       dphi_dx[1][i] * duna_dy + 
                       dphi_dx[1][i] * dvna_dx) * k2;
        
        double h1uy_= (dphi_dx[0][i] * duna_dy + 
                      2. * dphi_dx[1][i] * dvna_dy + 
                      dphi_dx[0][i] * dvna_dx) * k2; 

		rhsVector2[2*i] -= (l2ux_ + h1ux_) * WJ;
		rhsVector2[2*i+1] -= (l2uy_ + h1uy_) * WJ;


        //Stabilization terms
        double LLx = (dphi_dx[0][i] * lamx_dx + dphi_dx[1][i] * lamx_dy) * tARLQ_ / dens_;
        double LLy = (dphi_dx[0][i] * lamy_dx + dphi_dx[1][i] * lamy_dy) * tARLQ_ / dens_;

        arlequinStabVector[2*i] -= LLx * WJ;
        arlequinStabVector[2*i+1] -= LLy * WJ;

	};
};


//------------------------------------------------------------------------------
//-----------------------------RESIDUAL - RHS VECTOR----------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::getResidualVector_FEM(int &index,double *phi_, double **dphi_dx, double *rhsVector){
    
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();

    double ff[2];
    ff[0] = parameters->getFieldForce(0);
    ff[1] = parameters->getFieldForce(1);

    double axm_ = alpha_m * ax_ + (1. - alpha_m) * axprev_;
    double aym_ = alpha_m * ay_ + (1. - alpha_m) * ayprev_;

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
    double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * vmeshprev_);

    double duna_dx = (alpha_f * du_dx + (1. - alpha_f) * duprev_dx);
    double duna_dy = (alpha_f * du_dy + (1. - alpha_f) * duprev_dy);
    double dvna_dx = (alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx);
    double dvna_dy = (alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy);


    double wna_ = alpha_f * intPointWeightFunction_FEM[index] + (1. - alpha_f) * intPointWeightFunctionPrev_FEM[index];
    double WJ = weight_ * djac_ * wna_;

    
    for (int i = 0; i < 6; i++){
       
        double mx = phi_[i] * (axm_) * dens_ + 
                    (una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) * 
                    axm_ * tSUPG_ * dens_;
        double my = phi_[i] * (aym_) * dens_ + 
                    (una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) * 
                    aym_ * tSUPG_ * dens_;

        double Kx = 2. * dphi_dx[0][i] * duna_dx * visc_ + 
                    dphi_dx[1][i] * duna_dy * visc_ + 
                    dphi_dx[1][i] * dvna_dx * visc_ ;
        
        double Ky= dphi_dx[0][i] * duna_dy * visc_ + 
                   2. * dphi_dx[1][i] * dvna_dy * visc_ + 
                   dphi_dx[0][i] * dvna_dx * visc_;            

        double KLSx = dphi_dx[0][i] * (duna_dx + dvna_dy) * tLSIC_ * dens_;
        double KLSy = dphi_dx[1][i] * (duna_dx + dvna_dy) * tLSIC_ * dens_;

        double Cx = (duna_dx * una_ + 
                    duna_dy * vna_) * phi_[i] * dens_ +
                    (una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) *
                    (una_ * duna_dx + vna_ * duna_dy) * tSUPG_ *dens_;
        double Cy = (dvna_dx * una_ + 
                    dvna_dy * vna_) * phi_[i] * dens_ +
                    (una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) *
                    (una_ * dvna_dx + vna_ * dvna_dy) * tSUPG_ *dens_;
 
        double Px = - (dphi_dx[0][i] * p_) + 
                      ((dphi_dx[0][i] * una_ + dphi_dx[1][i] * vna_) * dp_dx * tSUPG_);
        double Py = - (dphi_dx[1][i] * p_) + 
                      ((dphi_dx[0][i] * una_ + dphi_dx[1][i] * vna_) * dp_dy * tSUPG_);
           
        double Q = ((duna_dx + dvna_dy) * phi_[i]) +
                    (dphi_dx[0][i] * dp_dx + dphi_dx[1][i] * dp_dy) * tPSPG_ / dens_ +
                     dphi_dx[0][i] * (una_ * duna_dx +
                                    vna_ * duna_dy) * tPSPG_ +
                    dphi_dx[1][i] * (una_ * dvna_dx +
                                    vna_ * dvna_dy) * tPSPG_ +
                    dphi_dx[0][i] * axm_* tPSPG_ +
                    dphi_dx[1][i] * aym_ * tPSPG_;
                   
        double Ffvx = phi_[i]* dens_ * ff[0] + 
                      tSUPG_* dens_ * (una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) * ff[0];
        double Ffvy = phi_[i]* dens_ * ff[1] + 
                      tSUPG_* dens_ * (una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) * ff[1];

        double Ffp = tPSPG_ * (dphi_dx[0][i] * ff[0] + dphi_dx[1][i] * ff[1]);
        

        rhsVector[2*i  ] += (-mx + Ffvx -Kx - Px - Cx - KLSx) * WJ;
        rhsVector[2*i+1] += (-my + Ffvy -Ky - Py - Cy - KLSy) * WJ;
        rhsVector[12+i] +=  (Ffp -Q) * WJ;

    };


    // //Stokes Problem
    // for (int i = 0; i < 6; i++){
       

    //     double Kx = 2. * dphi_dx[0][i] * duna_dx * visc_ + 
    //                 dphi_dx[1][i] * duna_dy * visc_ + 
    //                 dphi_dx[1][i] * dvna_dx * visc_ ;
        
    //     double Ky= dphi_dx[0][i] * duna_dy * visc_ + 
    //                2. * dphi_dx[1][i] * dvna_dy * visc_ + 
    //                dphi_dx[0][i] * dvna_dx * visc_;            


    //     double Px = - (dphi_dx[0][i] * p_);
    //     double Py = - (dphi_dx[1][i] * p_);
           
    //     double Q = ((duna_dx + dvna_dy) * phi_[i]) +
    //                 (dphi_dx[0][i] * dp_dx + dphi_dx[1][i] * dp_dy) * tPSPG_ / dens_;
                   
    //     double Ffvx = phi_[i]* dens_ * ff[0];
    //     double Ffvy = phi_[i]* dens_ * ff[1];

    //     double Ffp = tPSPG_ * (dphi_dx[0][i] * ff[0] + dphi_dx[1][i] * ff[1]);
        

    //     rhsVector[2*i  ] += (Ffvx -Kx - Px) * WJ;
    //     rhsVector[2*i+1] += (Ffvy -Ky - Py) * WJ;
    //     rhsVector[12+i] +=  (Ffp -Q) * WJ;

    // };

    return;
};

template<>
void Element<2>::getResidualVector_ISO(int &index, double *phi_, double **dphi_dx, double *rhsVector){
    
    double &visc_ = parameters->getViscosity();
    double &dens_ = parameters->getDensity();
    double &alpha_f = parameters->getAlphaF();
    double &alpha_m = parameters->getAlphaM();

    double ff[2];
    ff[0] = parameters->getFieldForce(0);
    ff[1] = parameters->getFieldForce(1);

    double axm_ = alpha_m * ax_ + (1. - alpha_m) * axprev_;
    double aym_ = alpha_m * ay_ + (1. - alpha_m) * ayprev_;

    double una_ = (alpha_f * u_ + (1. - alpha_f) * uprev_) - (alpha_f * umesh_ + (1. - alpha_f) * umeshprev_);
    double vna_ = (alpha_f * v_ + (1. - alpha_f) * vprev_) - (alpha_f * vmesh_ + (1. - alpha_f) * vmeshprev_);

    double duna_dx = (alpha_f * du_dx + (1. - alpha_f) * duprev_dx);
    double duna_dy = (alpha_f * du_dy + (1. - alpha_f) * duprev_dy);
    double dvna_dx = (alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx);
    double dvna_dy = (alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy);
   

    double wna_ = alpha_f * intPointWeightFunction_ISO[index] + (1. - alpha_f) * intPointWeightFunctionPrev_ISO[index];

    double WJ = weight_ * djac_ * wna_;


    for (int i = 0; i < 9; i++){


        double mx = phi_[i] * axm_ * dens_ + 
                    (una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) * 
                    (axm_) * tSUPG_ * dens_;
        double my = phi_[i] * aym_ * dens_ + 
                    (una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) * 
                    (aym_) * tSUPG_ * dens_;

        double Kx = 2. * dphi_dx[0][i] * duna_dx * visc_ + 
                    dphi_dx[1][i] * duna_dy * visc_ + 
                    dphi_dx[1][i] * dvna_dx * visc_;
        
        double Ky= dphi_dx[0][i] * duna_dy * visc_ + 
                   2. * dphi_dx[1][i] * dvna_dy * visc_ + 
                   dphi_dx[0][i] * dvna_dx * visc_;                
  

        double KLSx = dphi_dx[0][i] * (duna_dx + dvna_dy) * tLSIC_ * dens_;
        double KLSy = dphi_dx[1][i] * (duna_dx + dvna_dy) * tLSIC_ * dens_;


        double Cx = (duna_dx * una_ + duna_dy * vna_) * phi_[i] * dens_ +
                    (una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) *
                    (una_ * duna_dx + vna_ * duna_dy) * tSUPG_ *dens_;
        double Cy = (dvna_dx * una_ + dvna_dy * vna_) * phi_[i] * dens_ +
                    (una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) *
                    (una_ * dvna_dx + vna_ * dvna_dy) * tSUPG_ *dens_;
 
        double Px = - (dphi_dx[0][i] * p_) 
                    + ((dphi_dx[0][i] * una_ + dphi_dx[1][i] * vna_)* dp_dx * tSUPG_);
        double Py = - (dphi_dx[1][i] * p_) 
                    + ((dphi_dx[0][i] * una_ + dphi_dx[1][i] * vna_) * dp_dy * tSUPG_);
           

        double Q = ((duna_dx + dvna_dy) * phi_[i]) +
                    (dphi_dx[0][i] * dp_dx + dphi_dx[1][i] * dp_dy) * tPSPG_ / dens_ +
                    dphi_dx[0][i] * (una_ * duna_dx +
                                     vna_ * duna_dy) * tPSPG_ +
                    dphi_dx[1][i] * (una_ * dvna_dx +
                                     vna_ * dvna_dy) * tPSPG_ +    
                    dphi_dx[0][i] * axm_ * tPSPG_ +
                    dphi_dx[1][i] * aym_ * tPSPG_;
             
        
        double Ffvx = phi_[i]* dens_ * ff[0] + 
                      tSUPG_* dens_ * (una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) * ff[0];
        double Ffvy = phi_[i]* dens_ * ff[1] + 
                      tSUPG_* dens_ * (una_ * dphi_dx[0][i] + vna_ * dphi_dx[1][i]) * ff[1];

        double Ffp = tPSPG_ * (dphi_dx[0][i] * ff[0] + dphi_dx[1][i] * ff[1]);

    

        rhsVector[2*i  ] += (-mx + Ffvx - Kx - Px - Cx - KLSx) * WJ;
        rhsVector[2*i+1] += (-my + Ffvy - Ky - Py - Cy - KLSy) * WJ;
       
        rhsVector[18+i] +=  (-Q + Ffp) * WJ;

    };
    // //Stokes problem
    // for (int i = 0; i < 9; i++){

    //     double Kx = 2. * dphi_dx[0][i] * duna_dx * visc_ + 
    //                 dphi_dx[1][i] * duna_dy * visc_ + 
    //                 dphi_dx[1][i] * dvna_dx * visc_;
        
    //     double Ky = dphi_dx[0][i] * duna_dy * visc_ + 
    //                2. * dphi_dx[1][i] * dvna_dy * visc_ + 
    //                dphi_dx[0][i] * dvna_dx * visc_;                
  
    //     double Px = - (dphi_dx[0][i] * p_);
    //     double Py = - (dphi_dx[1][i] * p_);
           
    //     double Q = ((duna_dx + dvna_dy) * phi_[i]) +
    //                 (dphi_dx[0][i] * dp_dx + dphi_dx[1][i] * dp_dy) * tPSPG_ / dens_;
             
    //     double Ffvx = phi_[i]* dens_ * ff[0];
    //     double Ffvy = phi_[i]* dens_ * ff[1];

    //     double Ffp = tPSPG_ * (dphi_dx[0][i] * ff[0] + dphi_dx[1][i] * ff[1]);


    //     rhsVector[2*i  ] += (Ffvx - Kx - Px) * weight_ * djac_ * wna_;
    //     rhsVector[2*i+1] += (Ffvy - Ky - Py) * weight_ * djac_ * wna_;
       
    //     rhsVector[18+i] +=  (-Q + Ffp) * weight_ * djac_ * wna_;

    // };
   
    return;
};


//------------------------------------------------------------------------------
//----------------------ARLEQUIN ELEMENT LOCAL MATRIX/VECTOR--------------------
//------------------------------------------------------------------------------
// template<>
// void Element<2>::getLagrangeMultipliersSameMesh_FEM(double **lagrMultMatrix, double *rhsVector1, double *rhsVector2){


//     //quadrature and functions local classes
//     NormalQuad              nQuad = NormalQuad(); 
//     QuadShapeFunction<2>    shapeQuad;

//     //data for computation of IGA basis functions
//     double xsi[2], phi_[6];

//     double **dphi_dx;
//     dphi_dx = new double*[2];
//     for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[6];

//     double **ainv_;
//     ainv_ = new double*[2];
//     for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

//     double **Jac;
//     Jac = new double*[2];
//     for (int i = 0; i < 2; ++i) Jac[i] = new double[2];
    
//     int index = 0;
//     for(double* it = nQuad.beginFem(); it != nQuad.endFem(); it++){

//         //Defines the integration points adimentional coordinates
//         xsi[0] = nQuad.PointListFem(index,0);
//         xsi[1] = nQuad.PointListFem(index,1);
//         //Returns the quadrature integration weight
//         weight_ = nQuad.WeightListFem(index);

//         //Computes the velocity shape functions
//         shapeQuad.evaluateFem(xsi,phi_);
       
//         //Computes the jacobian matrix
//         getJacobianMatrix_FEM(xsi,Jac,ainv_);

//         //Computes spatial derivatives
//         getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

//         //Interpolates variables
//         getInterpolatedVariablesSameMesh_FEM(phi_,dphi_dx);

//         //Computes matrixes and vectors
//         getMatrixAndVectorsSameMesh_FEM(phi_,dphi_dx,lagrMultMatrix, 
//                                         rhsVector1,rhsVector2);
      
//         index++; 

  
//     };

//     for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
//     delete [] dphi_dx;
//     for (int i = 0; i < 2; ++i) delete [] ainv_[i];
//     delete [] ainv_;
//     for (int i = 0; i < 2; ++i) delete [] Jac[i];
//     delete [] Jac;

//     return;
      
// };

template<>
void Element<2>::getLagrangeMultipliersSameMesh_FEM(int &index, double **lagrMultMatrix, double *rhsVector1, double *rhsVector2){


    //quadrature and functions local classes
    NormalQuad              nQuad = NormalQuad(); 
    QuadShapeFunction<2>    shapeQuad;

    //data for computation of IGA basis functions
    double xsi[2], phi_[6];

    double **dphi_dx;
    dphi_dx = new double*[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[6];

    double **ainv_;
    ainv_ = new double*[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **Jac;
    Jac = new double*[2];
    for (int i = 0; i < 2; ++i) Jac[i] = new double[2];
    

    //Defines the integration points adimentional coordinates
    xsi[0] = nQuad.PointListFem(index,0);
    xsi[1] = nQuad.PointListFem(index,1);
    //Returns the quadrature integration weight
    weight_ = nQuad.WeightListFem(index);

    //Computes the velocity shape functions
    shapeQuad.evaluateFem(xsi,phi_);
   
    //Computes the jacobian matrix
    getJacobianMatrix_FEM(xsi,Jac,ainv_);

    //Computes spatial derivatives
    getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

    //Interpolates variables
    getInterpolatedVariablesSameMesh_FEM(phi_,dphi_dx);

    //Computes matrixes and vectors
    getMatrixAndVectorsSameMesh_FEM(phi_,dphi_dx,lagrMultMatrix, 
                                    rhsVector1,rhsVector2);
      

    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] Jac[i];
    delete [] Jac;

    return;
      
};

template<>
void Element<2>::getLagrangeMultipliersSameMesh_tSUPG_tPSPG_FEM(double **jacobianNRMatrix, double *rhsVector){


	//quadrature and functions local classes
    NormalQuad              nQuad = NormalQuad(); 
    QuadShapeFunction<2>    shapeQuad;

    //data for computation of IGA basis functions
    double xsi[2], phi_[6];

    double **dphi_dx;
    dphi_dx = new double*[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[6];

    double **ainv_;
    ainv_ = new double*[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **Jac;
    Jac = new double*[2];
    for (int i = 0; i < 2; ++i) Jac[i] = new double[2];
    
    int index = 0;
    for(double* it = nQuad.beginFem(); it != nQuad.endFem(); it++){

        //Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointListFem(index,0);
        xsi[1] = nQuad.PointListFem(index,1);
        //Returns the quadrature integration weight
        weight_ = nQuad.WeightListFem(index);

        //Computes the velocity shape functions
        shapeQuad.evaluateFem(xsi,phi_);
       
        //Computes the jacobian matrix
        getJacobianMatrix_FEM(xsi,Jac,ainv_);

        //Computes spatial derivatives
        getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

        //Interpolates variables
        getInterpolatedVariablesSameMesh_FEM(phi_,dphi_dx);

        //Computes stabilization parameters
        getNewParameterSUPG_FEM(Jac,phi_,dphi_dx);

        //Computes matrixes and vectors
        getMatrixAndVectorsSameMesh_tSUPG_tPSPG_FEM(phi_,dphi_dx,jacobianNRMatrix,rhsVector);
      
        index++; 

  
    };

    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] Jac[i];
    delete [] Jac;

    return;
}

template<>
void Element<2>::getLagrangeMultipliersSameMeshArlqStab_FEM(int &index, int &ipatchC, std::vector<Nodes *> &nodesCoarse_,int *connecC,
                                                            std::vector<IParameters *> &iparamC,double **arlequinStabD, double *arlequinStabVectorD,
                                                            double **arlequinStab1, double *arlequinStabVector1){

    //quadrature and functions local classes
    NormalQuad              nQuad = NormalQuad(); 
    QuadShapeFunction<2>    shapeQuad;


    //Data for IGA coarse mesh computation
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[9],xsiC[2], phiC_[9];
    for (int i = 0; i < 9; i ++) wpcC[i] = (*nodesC_)[connectC_[i]] -> getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[8]] -> getINC(); 
  

    //data for computation of FEM fine mesh computations
    double xsi[2], phi_[6];

    double **dphi_dx;
    dphi_dx = new double*[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[6];

    double ***ddphi_dx;
    ddphi_dx = new double**[2];
    for (int i = 0; i < 2; ++i) {
        ddphi_dx[i] = new double*[2];
        for (int j = 0; j < 2; j++) ddphi_dx[i][j] = new double[6];
    };

    double **ainv_;
    ainv_ = new double*[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **Jac;
    Jac = new double*[2];
    for (int i = 0; i < 2; ++i) Jac[i] = new double[2];


    //Defines the integration points adimentional coordinates
    xsi[0] = nQuad.PointListFem(index,0);
    xsi[1] = nQuad.PointListFem(index,1);
    //Returns the quadrature integration weight
    weight_ = nQuad.WeightListFem(index);

    //Computes the shape functions
    shapeQuad.evaluateFem(xsi,phi_);
   
    //Computes the jacobian matrix
    getJacobianMatrix_FEM(xsi,Jac,ainv_);

    //Computes spatial derivatives
    getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

    //Computes spatial second derivatives
    getSecondSpatialDerivatives_FEM(xsi,ainv_,ddphi_dx);


    
    //Integration point in the coarse mesh element
    for (int k = 0; k < 2; k++) xsiC[k] = intPointCorrespXsi_FEM[index][k];

    //Computes the coqrse shape functions
    shapeQuad.evaluateIso(xsiC,phiC_,wpcC,incC_,(*iparametersC),NpatchC_);

    

    //Interpolates variables
    getInterpolatedVariablesSameMeshArlqStab_FEM(phi_,dphi_dx,ddphi_dx);


    getParameterArlequin2_FEM(phi_,phiC_,dphi_dx,ddphi_dx);
    //get Arlequin stabilization parameter
    if (iTimeStep < 1) getParameterArlequin_FEM(phi_,dphi_dx);

    // getParameterArlequin_FEM(phi_,dphi_dx);

    //Computes matrixes and vectors
    getMatrixAndVectorsSameMeshArlqStab_FEM(index,phi_,dphi_dx,ddphi_dx,arlequinStabD,arlequinStabVectorD,
                                            arlequinStab1,arlequinStabVector1);
      
  

    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) delete [] ddphi_dx[i][j];
        delete [] ddphi_dx[i];
    };
    delete [] ddphi_dx;

    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] Jac[i];
    delete [] Jac;
};



// template<>
// void Element<2>::getLagrangeMultipliersSameMeshArlqStab_FEM(double **arlequinStabD, double *arlequinStabVectorD,
//                                                             double **arlequinStab1, double *arlequinStabVector1){

    //     //quadrature and functions local classes
    // NormalQuad              nQuad = NormalQuad(); 
    // QuadShapeFunction<2>    shapeQuad;

    // //data for computation of IGA basis functions
    // double xsi[2], phi_[6];

    // double **dphi_dx;
    // dphi_dx = new double*[2];
    // for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[6];

    // double ***ddphi_dx;
    // ddphi_dx = new double**[2];
    // for (int i = 0; i < 2; ++i) {
    //     ddphi_dx[i] = new double*[2];
    //     for (int j = 0; j < 2; j++) ddphi_dx[i][j] = new double[6];
    // };

    // double **ainv_;
    // ainv_ = new double*[2];
    // for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    // double **Jac;
    // Jac = new double*[2];
    // for (int i = 0; i < 2; ++i) Jac[i] = new double[2];
    
    // int index = 0;
    // for(double* it = nQuad.beginFem(); it != nQuad.endFem(); it++){

    //     //Defines the integration points adimentional coordinates
    //     xsi[0] = nQuad.PointListFem(index,0);
    //     xsi[1] = nQuad.PointListFem(index,1);
    //     //Returns the quadrature integration weight
    //     weight_ = nQuad.WeightListFem(index);

    //     //Computes the velocity shape functions
    //     shapeQuad.evaluateFem(xsi,phi_);
       
    //     //Computes the jacobian matrix
    //     getJacobianMatrix_FEM(xsi,Jac,ainv_);

    //     //Computes spatial derivatives
    //     getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

    //     //Computes spatial second derivatives
    //     getSecondSpatialDerivatives_FEM(xsi,ainv_,ddphi_dx);

    //     //Interpolates variables
    //     getInterpolatedVariablesSameMeshArlqStab_FEM(phi_,dphi_dx,ddphi_dx);

    //     // getParameterArlequin2_FEM(phi_,dphi_dx,ddphi_dx);
    //     // //get Arlequin stabilization parameter
    //     // if (iTimeStep < 1) getParameterArlequin_FEM(phi_,dphi_dx);

    //     getParameterArlequin_FEM(phi_,dphi_dx);

    //     //Computes matrixes and vectors
    //     getMatrixAndVectorsSameMeshArlqStab_FEM(index,phi_,dphi_dx,ddphi_dx,arlequinStabD,arlequinStabVectorD,
    //                                             arlequinStab1,arlequinStabVector1);
      
    //     index++; 

  
    // };

    // for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    // delete [] dphi_dx;
    
    // for (int i = 0; i < 2; ++i) {
    //     for (int j = 0; j < 2; ++j) delete [] ddphi_dx[i][j];
    //     delete [] ddphi_dx[i];
    // };
    // delete [] ddphi_dx;

    // for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    // delete [] ainv_;
    // for (int i = 0; i < 2; ++i) delete [] Jac[i];
    // delete [] Jac;

    // return;

// };


template<>
void Element<2>::getLagrangeMultipliersSameMesh_ISO(double **lagrMultMatrix, double *rhsVector1, double *rhsVector2){

	//quadrature and functions local classes
    NormalQuad              nQuad = NormalQuad(); 
    QuadShapeFunction<2>    shapeQuad;
    
    //data for computation of IGA basis functions
    double wpc[9],xsi[2], phi_[9];
    for (int i = 0; i < 9; i ++) wpc[i] = (*nodes_)[connect_[i]] -> getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]] -> getINC();  

    double **dphi_dx;
    dphi_dx = new double*[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double*[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **quadJacMat;
    quadJacMat = new double*[2];
    for (int i = 0; i < 2; ++i) quadJacMat[i] = new double[2];

    int index = 0;
    for(double* it = nQuad.beginIso(); it != nQuad.endIso(); it++){

        //Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointListIso(index,0);
        xsi[1] = nQuad.PointListIso(index,1);
        //Returns the quadrature integration weight
        weight_ = nQuad.WeightListIso(index);

        //Computes the velocity shape functions
        shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,(*iparameters),Npatch_);
       
        //Computes the jacobian matrix
        getJacobianMatrix_ISO(xsi,quadJacMat,ainv_);

        //Computes spatial derivatives
        getSpatialDerivatives_ISO(xsi,ainv_,dphi_dx);

        //Interpolates variables
        getInterpolatedVariablesSameMesh_ISO(phi_,dphi_dx);

        //Computes matrixes and vectors
        getMatrixAndVectorsSameMesh_ISO(phi_,dphi_dx,lagrMultMatrix, 
        	                            rhsVector1,rhsVector2);

       
        index++; 
    };


    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] quadJacMat[i];
    delete [] quadJacMat;

    return;
      
};

template<>
void Element<2>::getLagrangeMultipliersSameMesh_tSUPG_tPSPG_ISO(double **jacobianNRMatrix, double *rhsVector){


    //quadrature and functions local classes
    NormalQuad              nQuad = NormalQuad(); 
    QuadShapeFunction<2>    shapeQuad;

    //data for computation of IGA basis functions
    double wpc[9],xsi[2], phi_[9];
    for (int i = 0; i < 9; i ++) wpc[i] = (*nodes_)[connect_[i]] -> getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]] -> getINC();  

    double **dphi_dx;
    dphi_dx = new double*[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double*[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **quadJacMat;
    quadJacMat = new double*[2];
    for (int i = 0; i < 2; ++i) quadJacMat[i] = new double[2];
    

    int index = 0;
    for(double* it = nQuad.beginIso(); it != nQuad.endIso(); it++){

        //Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointListIso(index,0);
        xsi[1] = nQuad.PointListIso(index,1);
        //Returns the quadrature integration weight
        weight_ = nQuad.WeightListIso(index);

        //Computes the velocity shape functions
        shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,(*iparameters),Npatch_);
       
        //Computes the jacobian matrix
        getJacobianMatrix_ISO(xsi,quadJacMat,ainv_);

        //Computes spatial derivatives
        getSpatialDerivatives_ISO(xsi,ainv_,dphi_dx);

        //Interpolates variables
        getInterpolatedVariablesSameMesh_ISO(phi_,dphi_dx);
       
        //Computes stabilization parameters
        getNewParameterSUPG_ISO(quadJacMat,phi_,dphi_dx);

        //Computes matrixes and vectors
        getMatrixAndVectorsSameMesh_tSUPG_tPSPG_ISO(phi_,dphi_dx,jacobianNRMatrix,rhsVector);
      
        index++; 

    };

    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] quadJacMat[i];
    delete [] quadJacMat;

    return;
};



template<>
void Element<2>::getLagrangeMultipliersSameMeshArlqStab_ISO(double **arlequinStabD, double *arlequinStabVectorD,
                                                            double **arlequinStab1, double *arlequinStabVector1){

        //quadrature and functions local classes
    NormalQuad              nQuad = NormalQuad(); 
    QuadShapeFunction<2>    shapeQuad;


    //data for computation of IGA basis functions
    double wpc[9],xsi[2], phi_[9];
    for (int i = 0; i < 9; i ++) wpc[i] = (*nodes_)[connect_[i]] -> getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]] -> getINC();  

    double **dphi_dx;
    dphi_dx = new double*[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double*[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **quadJacMat;
    quadJacMat = new double*[2];
    for (int i = 0; i < 2; ++i) quadJacMat[i] = new double[2];

    double ***ddphi_dx;
    ddphi_dx = new double**[2];
    for (int i = 0; i < 2; ++i) {
        ddphi_dx[i] = new double*[2];
        for (int j = 0; j < 2; j++) ddphi_dx[i][j] = new double[9];
    };
  
    int index = 0;
    for(double* it = nQuad.beginIso(); it != nQuad.endIso(); it++){

       
        //Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointListIso(index,0);
        xsi[1] = nQuad.PointListIso(index,1);
        //Returns the quadrature integration weight
        weight_ = nQuad.WeightListIso(index);

        //Computes the velocity shape functions
        shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,(*iparameters),Npatch_);
       
        //Computes the jacobian matrix
        getJacobianMatrix_ISO(xsi,quadJacMat,ainv_);

        //Computes spatial derivatives
        getSpatialDerivatives_ISO(xsi,ainv_,dphi_dx);

        //Computes spatial second derivatives
        getSecondSpatialDerivatives_ISO(xsi,ainv_,ddphi_dx);

        //get Arlequin stabilization parameter
        getParameterArlequin_ISO(phi_,dphi_dx);

        //Interpolates variables
        getInterpolatedVariablesSameMeshArlqStab_ISO(phi_,dphi_dx,ddphi_dx);

        //Computes matrixes and vectors
        getMatrixAndVectorsSameMeshArlqStab_ISO(index,phi_,dphi_dx,ddphi_dx,arlequinStabD,arlequinStabVectorD,
                                                arlequinStab1,arlequinStabVector1);
      
        index++; 

  
    };


    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] quadJacMat[i];
    delete [] quadJacMat;
    
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j) delete [] ddphi_dx[i][j];
        delete [] ddphi_dx[i];
    };
    delete [] ddphi_dx;

    return;

};



template<>
void Element<2>::getLagrangeMultipliersDifferentMesh_FEM_FEM(std::vector<Nodes *> &nodesCoarse_,int *connecC,int &ielem,
                                                             double **lagrMultMatrix, double *rhsVector1, double *rhsVector2,
                                                             double **arlequinStab, double *arlequinStabVector){


    int dim = 2;
    //quadrature and functions local classes
    SpecialQuad              sQuad = SpecialQuad(); 
    QuadShapeFunction<2>     shapeQuad;
    
    //Data for FEM coarse mesh computation (Velocity field)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    
    double xsiC[dim], phiC_[6];

    double **dphiC_dx;
    dphiC_dx = new double*[dim];
    for (int i = 0; i < dim; ++i) dphiC_dx[i] = new double[6];

    double **ainvC_;
    ainvC_ = new double*[dim];
    for (int i = 0; i < dim; ++i) ainvC_[i] = new double[dim];

    double **JacC;
    JacC = new double*[dim];
    for (int i = 0; i < dim; ++i) JacC[i] = new double[dim];

   
   //Data for FEM fine mesh computation (Lagrange field)
    double xsi[dim], phi_[6];

    double **dphi_dx;
    dphi_dx = new double*[dim];
    for (int i = 0; i < dim; ++i) dphi_dx[i] = new double[6];

    double **ainv_;
    ainv_ = new double*[dim];
    for (int i = 0; i < dim; ++i) ainv_[i] = new double[dim];

    double **Jac;
    Jac = new double*[dim];
    for (int i = 0; i < dim; ++i) Jac[i] = new double[dim];

    int index = 0;    
    for(double* it = sQuad.beginFem(); it != sQuad.endFem(); it++){
        
        if ((intPointCorrespElem_FEM[index] == ielem)){

            //Fine mesh computation
            //Defines the integration points adimentional coordinates
            for (int k = 0; k < dim; k ++) xsi[k] = sQuad.PointListFem(index,k);
            
            //Returns the quadrature integration weight
            weight_ = sQuad.WeightListFem(index);

            //Computes the velocity shape functions
            shapeQuad.evaluateFem(xsi,phi_);
        
            //Computes the jacobian matrix
            getJacobianMatrix_FEM(xsi,Jac,ainv_);

            //Computes spatial derivatives
            getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

            //Computes Arlequin Stabilization term
            getParameterArlequin_FEM(phi_,dphi_dx);

           
            //Coarse mesh computatation
            //Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < dim; k++) xsiC[k] = intPointCorrespXsi_FEM[index][k];

            //Computes the velocity shape functions
            shapeQuad.evaluateFem(xsiC,phiC_);

            //Computes the jacobian matrix
            getJacobianMatrix_COARSE_FEM(xsiC,JacC,ainvC_);

            //Computes spatial derivatives
            getSpatialDerivatives_FEM(xsiC,ainvC_,dphiC_dx);



            //Interpolates Lagrange multipliers  and velocity field
            getInterpolatedVariablesDifferentMesh_FEM_FEM(phi_,dphi_dx,phiC_,dphiC_dx);

            //Computes Matrixes and vectors
            getMatrixAndVectorsDifferentMesh_FEM_FEM(phi_,dphi_dx,phiC_,dphiC_dx,
                                                     lagrMultMatrix,rhsVector1,rhsVector2,
                                                     arlequinStab,arlequinStabVector);


        };
        index++;        
    };  

    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] Jac[i];
    delete [] Jac;

    for (int i = 0; i < 2; ++i) delete [] dphiC_dx[i];
    delete [] dphiC_dx;
    for (int i = 0; i < 2; ++i) delete [] ainvC_[i];
    delete [] ainvC_;
    for (int i = 0; i < 2; ++i) delete [] JacC[i];
    delete [] JacC;  

    return;

};

template<>
void Element<2>::getLagrangeMultipliersDifferentMesh_FEM_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_,int *connecC,
                                                            std::vector<IParameters *> &iparamC, int &ielem,
                                                            double **lagrMultMatrix, double *rhsVector1, double *rhsVector2){


    int dim = 2;
    //quadrature and functions local classes
    SpecialQuad              sQuad = SpecialQuad(); 
    QuadShapeFunction<2>     shapeQuad;
    
    //Data for IGA coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[9],xsiC[dim], phiC_[9];
    for (int i = 0; i < 9; i ++) wpcC[i] = (*nodesC_)[connectC_[i]] -> getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[8]] -> getINC(); 

    double **dphiC_dx;
    dphiC_dx = new double*[dim];
    for (int i = 0; i < dim; ++i) dphiC_dx[i] = new double[9];

    double **ainvC_;
    ainvC_ = new double*[dim];
    for (int i = 0; i < dim; ++i) ainvC_[i] = new double[dim];

    double **quadJacMatC;
    quadJacMatC = new double*[dim];
    for (int i = 0; i < dim; ++i) quadJacMatC[i] = new double[dim];

   
   //Data for FEM fine mesh computation (Lagrange field)
    double xsi[dim], phi_[6];

    double **dphi_dx;
    dphi_dx = new double*[dim];
    for (int i = 0; i < dim; ++i) dphi_dx[i] = new double[6];

    double **ainv_;
    ainv_ = new double*[dim];
    for (int i = 0; i < dim; ++i) ainv_[i] = new double[dim];

    double **Jac;
    Jac = new double*[dim];
    for (int i = 0; i < dim; ++i) Jac[i] = new double[dim];

    int index = 0;    
    for(double* it = sQuad.beginFem(); it != sQuad.endFem(); it++){
        
        if ((intPointCorrespElem_FEM[index] == ielem)){

            //Fine mesh computation
            //Defines the integration points adimentional coordinates
            for (int k = 0; k < dim; k ++) xsi[k] = sQuad.PointListFem(index,k);
            
            //Returns the quadrature integration weight
            weight_ = sQuad.WeightListFem(index);

            //Computes the velocity shape functions
            shapeQuad.evaluateFem(xsi,phi_);
        
            //Computes the jacobian matrix
            getJacobianMatrix_FEM(xsi,Jac,ainv_);

            //Computes spatial derivatives
            getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

           
            //Coarse mesh computatation
            //Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < dim; k++) xsiC[k] = intPointCorrespXsi_FEM[index][k];

            //Computes the velocity shape functions
            shapeQuad.evaluateIso(xsiC,phiC_,wpcC,incC_,(*iparametersC),NpatchC_);

            //Computes the jacobian matrix
            getJacobianMatrix_COARSE_ISO(xsiC,quadJacMatC,ainvC_);

            //Computes spatial derivatives
            getSpatialDerivatives_COARSE_ISO(xsiC,ainvC_,dphiC_dx);

            

            //Interpolates Lagrange multiplier and velocity
            getInterpolatedVariablesDifferentMesh_FEM_ISO(phi_,dphi_dx,phiC_,dphiC_dx);

            //Computes Matrix and vectors
            getMatrixAndVectorsDifferentMesh_FEM_ISO(phi_,dphi_dx,phiC_,dphiC_dx,
                                                     lagrMultMatrix,rhsVector1,rhsVector2);


        };
        index++;        
    };  

    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] Jac[i];
    delete [] Jac;

    for (int i = 0; i < 2; ++i) delete [] dphiC_dx[i];
    delete [] dphiC_dx;
    for (int i = 0; i < 2; ++i) delete [] ainvC_[i];
    delete [] ainvC_;
    for (int i = 0; i < 2; ++i) delete [] quadJacMatC[i];
    delete [] quadJacMatC;  

    return;

};


template<>
void Element<2>::getLagrangeMultipliersDifferentMesh_tSUPG_tPSPG_FEM_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_,int *connecC,
                                                    					 std::vector<IParameters *> &iparamC, int &ielem, 
                                                    		    		 double **jacobianNRMatrix, double *rhsVector){

	int dim = 2;
    //quadrature and functions local classes
    SpecialQuad              sQuad = SpecialQuad(); 
    QuadShapeFunction<2>     shapeQuad;
    
    //Data for IGA coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[9],xsiC[dim], phiC_[9];
    for (int i = 0; i < 9; i ++) wpcC[i] = (*nodesC_)[connectC_[i]] -> getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[8]] -> getINC(); 

    double **dphiC_dx;
    dphiC_dx = new double*[dim];
    for (int i = 0; i < dim; ++i) dphiC_dx[i] = new double[9];

    double **ainvC_;
    ainvC_ = new double*[dim];
    for (int i = 0; i < dim; ++i) ainvC_[i] = new double[dim];

    double **quadJacMatC;
    quadJacMatC = new double*[dim];
    for (int i = 0; i < dim; ++i) quadJacMatC[i] = new double[dim];

   
   //Data for FEM fine mesh computation (Lagrange field)
    double xsi[dim], phi_[6];

    double **dphi_dx;
    dphi_dx = new double*[dim];
    for (int i = 0; i < dim; ++i) dphi_dx[i] = new double[6];

    double **ainv_;
    ainv_ = new double*[dim];
    for (int i = 0; i < dim; ++i) ainv_[i] = new double[dim];

    double **Jac;
    Jac = new double*[dim];
    for (int i = 0; i < dim; ++i) Jac[i] = new double[dim];

    int index = 0;    
    for(double* it = sQuad.beginFem(); it != sQuad.endFem(); it++){
        
        if ((intPointCorrespElem_FEM[index] == ielem)){

            //Fine mesh computation
            //Defines the integration points adimentional coordinates
            for (int k = 0; k < dim; k ++) xsi[k] = sQuad.PointListFem(index,k);
            
            //Returns the quadrature integration weight
            weight_ = sQuad.WeightListFem(index);

            //Computes the velocity shape functions
            shapeQuad.evaluateFem(xsi,phi_);
        
            //Computes the jacobian matrix
            getJacobianMatrix_FEM(xsi,Jac,ainv_);

            //Computes spatial derivatives
            getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

            //Computes stabilization parameters
            getNewParameterSUPG_FEM(Jac,phi_,dphi_dx);
           
            //Coarse mesh computatation
            //Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < dim; k++) xsiC[k] = intPointCorrespXsi_FEM[index][k];

            //Computes the jacobian matrix
            getJacobianMatrix_COARSE_ISO(xsiC,quadJacMatC,ainvC_);

            //Computes spatial derivatives
            getSpatialDerivatives_COARSE_ISO(xsiC,ainvC_,dphiC_dx);

            


            //Interpolates Lagrange multiplier and velocity
            getInterpolatedVariablesDifferentMesh_FEM_ISO(phi_,dphi_dx,phiC_,dphiC_dx);

            //Computes Matrix and vectors
            getMatrixAndVectorsDifferentMesh_tSUPG_tPSPG_FEM_ISO(phi_,dphiC_dx,
                                                                jacobianNRMatrix,rhsVector);


        };
        index++;        
    };  

    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] Jac[i];
    delete [] Jac;

    for (int i = 0; i < 2; ++i) delete [] dphiC_dx[i];
    delete [] dphiC_dx;
    for (int i = 0; i < 2; ++i) delete [] ainvC_[i];
    delete [] ainvC_;
    for (int i = 0; i < 2; ++i) delete [] quadJacMatC[i];
    delete [] quadJacMatC;  

    return;


};


template<>
void Element<2>::getLagrangeMultipliersDifferentMeshArlqStab_FEM_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_,int *connecC,
                                                                     std::vector<IParameters *> &iparamC, int &ielem,
                                                                     double **arlequinStabD, double *arlequinStabVectorD,
                                                                     double **arlequinStab0, double *arlequinStabVector0){


    int dim = 2;
    //quadrature and functions local classes
    SpecialQuad              sQuad = SpecialQuad(); 
    QuadShapeFunction<2>     shapeQuad;
    
    //Data for IGA coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[9],xsiC[dim], phiC_[9];
    for (int i = 0; i < 9; i ++) wpcC[i] = (*nodesC_)[connectC_[i]] -> getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[8]] -> getINC(); 

    double **dphiC_dx;
    dphiC_dx = new double*[dim];
    for (int i = 0; i < dim; ++i) dphiC_dx[i] = new double[9];

    double ***ddphiC_dx;
    ddphiC_dx = new double**[dim];
    for (int i = 0; i < dim; ++i) {
        ddphiC_dx[i] = new double*[2];
        for (int j = 0; j < dim; j++) ddphiC_dx[i][j] = new double[9];
    };    

    double **ainvC_;
    ainvC_ = new double*[dim];
    for (int i = 0; i < dim; ++i) ainvC_[i] = new double[dim];

    double **quadJacMatC;
    quadJacMatC = new double*[dim];
    for (int i = 0; i < dim; ++i) quadJacMatC[i] = new double[dim];

   
   //Data for FEM fine mesh computation (Lagrange field)
    double xsi[dim], phi_[6];

    double **dphi_dx;
    dphi_dx = new double*[dim];
    for (int i = 0; i < dim; ++i) dphi_dx[i] = new double[6];

     double ***ddphi_dx;
    ddphi_dx = new double**[dim];
    for (int i = 0; i < dim; ++i) {
        ddphi_dx[i] = new double*[2];
        for (int j = 0; j < dim; j++) ddphi_dx[i][j] = new double[6];
    };  

    double **ainv_;
    ainv_ = new double*[dim];
    for (int i = 0; i < dim; ++i) ainv_[i] = new double[dim];

    double **Jac;
    Jac = new double*[dim];
    for (int i = 0; i < dim; ++i) Jac[i] = new double[dim];

    int index = 0;    
    for(double* it = sQuad.beginFem(); it != sQuad.endFem(); it++){
        
        if ((intPointCorrespElem_FEM[index] == ielem)){

            //Fine mesh computation
            //Defines the integration points adimentional coordinates
            for (int k = 0; k < dim; k ++) xsi[k] = sQuad.PointListFem(index,k);
            
            //Returns the quadrature integration weight
            weight_ = sQuad.WeightListFem(index);

            //Computes the velocity shape functions
            shapeQuad.evaluateFem(xsi,phi_);
        
            //Computes the jacobian matrix
            getJacobianMatrix_FEM(xsi,Jac,ainv_);

            //Computes spatial derivatives
            getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

	        //Computes spatial second derivatives
	        getSecondSpatialDerivatives_FEM(xsi,ainv_,ddphi_dx);

	        //Interpolates variables
	        getInterpolatedVariablesSameMeshArlqStab_FEM(phi_,dphi_dx,ddphi_dx);

	       
            // //get Arlequin stabilization parameter
            // getParameterArlequin_FEM(phi_,dphi_dx);            
           
            //Coarse mesh computatation
            //Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < dim; k++) xsiC[k] = intPointCorrespXsi_FEM[index][k];

            //Computes the velocity shape functions
            shapeQuad.evaluateIso(xsiC,phiC_,wpcC,incC_,(*iparametersC),NpatchC_);

            //Computes the jacobian matrix
            getJacobianMatrix_COARSE_ISO(xsiC,quadJacMatC,ainvC_);

            //Computes spatial derivatives
            getSpatialDerivatives_COARSE_ISO(xsiC,ainvC_,dphiC_dx);

            //Computes second spatial derivatives
            getSecondSpatialDerivatives_COARSE_ISO(xsiC,ainvC_,ddphiC_dx);

             //get Arlequin stabilization parameter (Usa velocidades interpoladas da malha FINE/Unico parametro para malha fine e para malha coarse)
            getParameterArlequin2_FEM(phi_,phiC_,dphi_dx,ddphi_dx);
            if (iTimeStep < 1) getParameterArlequin_FEM(phi_,dphi_dx);


            //Interpolates Lagrange multiplier and velocity
            getInterpolatedVariablesDifferentMeshArlqStab_FEM_ISO(dphi_dx,phiC_,dphiC_dx,ddphiC_dx);




            //Computes Matrix and vectors
            getMatrixAndVectorsDifferentMeshArlqStab_FEM_ISO(index,dphi_dx,phiC_,dphiC_dx,ddphiC_dx,
                                                             arlequinStabD,arlequinStabVectorD,
                                                             arlequinStab0,arlequinStabVector0);


        };
        index++;        
    };  

    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] Jac[i];
    delete [] Jac;

    for (int i = 0; i < 2; ++i) delete [] dphiC_dx[i];
    delete [] dphiC_dx;

    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < 2; ++j){
        	ddphiC_dx[i][j];
        	ddphi_dx[i][j];
        } 
        delete [] ddphiC_dx[i];
        delete [] ddphi_dx[i];
    }
    delete [] ddphiC_dx;
    delete [] ddphi_dx;



    for (int i = 0; i < 2; ++i) delete [] ainvC_[i];
    delete [] ainvC_;
    for (int i = 0; i < 2; ++i) delete [] quadJacMatC[i];
    delete [] quadJacMatC;  

    return;

};

template<>
void Element<2>::getLagrangeMultipliersDifferentMesh_ISO_FEM(std::vector<Nodes *> &nodesCoarse_,int *connecC,int &ielem, 
                                                             double **lagrMultMatrix,double *rhsVector1, double *rhsVector2){
    



    int dim = 2;
    //quadrature and functions local classes
    SpecialQuad              sQuad = SpecialQuad(); 
    QuadShapeFunction<2>     shapeQuad;
    
    //Data for FEM coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;

    double xsiC[dim], phiC_[6];

    double **dphiC_dx;
    dphiC_dx = new double*[dim];
    for (int i = 0; i < dim; ++i) dphiC_dx[i] = new double[6];

    double **ainvC_;
    ainvC_ = new double*[dim];
    for (int i = 0; i < dim; ++i) ainvC_[i] = new double[dim];

    double **JacC;
    JacC = new double*[dim];
    for (int i = 0; i < dim; ++i) JacC[i] = new double[dim];

   
    //Data for IGA fine mesh computation (Lagrange field)
    double wpc[9],xsi[dim], phi_[9];
    for (int i = 0; i < 9; i ++) wpc[i] = (*nodes_)[connect_[i]] -> getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]] -> getINC(); 

    double **dphi_dx;
    dphi_dx = new double*[dim];
    for (int i = 0; i < dim; ++i) dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double*[dim];
    for (int i = 0; i < dim; ++i) ainv_[i] = new double[dim];

    double **quadJacMat;
    quadJacMat = new double*[dim];
    for (int i = 0; i < dim; ++i) quadJacMat[i] = new double[dim];

    int index = 0;    
    for(double* it = sQuad.beginIso(); it != sQuad.endIso(); it++){
        
        if ((intPointCorrespElem_ISO[index] == ielem)){

            
            //Fine mesh computation
            //Defines the integration points adimentional coordinates
            for (int k = 0; k < dim; k ++) xsi[k] = sQuad.PointListIso(index,k);
            
            //Returns the quadrature integration weight
            weight_ = sQuad.WeightListIso(index);

            //Computes the velocity shape functions
            shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,(*iparameters),Npatch_);
        
            //Computes the jacobian matrix
            getJacobianMatrix_ISO(xsi,quadJacMat,ainv_);

            //Computes spatial derivatives
            getSpatialDerivatives_ISO(xsi,ainv_,dphi_dx);

           
            
            //Coarse mesh computatation
            //Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < dim; k++) xsiC[k] = intPointCorrespXsi_ISO[index][k];

            //Computes the velocity shape function
            shapeQuad.evaluateFem(xsiC,phiC_);

            //Computes the jacobian matrix
            getJacobianMatrix_COARSE_FEM(xsiC,JacC,ainvC_);

            //Computes spatial derivatives
            getSpatialDerivatives_FEM(xsiC, ainvC_, dphiC_dx);

            

            //Interpolates Lagrange multiplier and velocity
            getInterpolatedVariablesDifferentMesh_ISO_FEM(phi_,dphi_dx,phiC_,dphiC_dx);

            //Computes Matrix and vectors
            getMatrixAndVectorsDifferentMesh_ISO_FEM(phi_,dphi_dx,phiC_,dphiC_dx,
                                                     lagrMultMatrix,rhsVector1,rhsVector2);


        };
        index++;        
    };  

    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] quadJacMat[i];
    delete [] quadJacMat;

    for (int i = 0; i < 2; ++i) delete [] dphiC_dx[i];
    delete [] dphiC_dx;
    for (int i = 0; i < 2; ++i) delete [] ainvC_[i];
    delete [] ainvC_;
    for (int i = 0; i < 2; ++i) delete [] JacC[i];
    delete [] JacC;  

    return;
};

template<>
void Element<2>::getLagrangeMultipliersDifferentMesh_tSUPG_tPSPG_ISO_FEM(std::vector<Nodes *> &nodesCoarse_,int *connecC,int &ielem, 
                                                                         double **jacobianNRMatrix, double *rhsVector){



    int dim = 2;
    //quadrature and functions local classes
    SpecialQuad              sQuad = SpecialQuad(); 
    QuadShapeFunction<2>     shapeQuad;
    
    //Data for FEM coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;

    double xsiC[dim], phiC_[6];

    double **dphiC_dx;
    dphiC_dx = new double*[dim];
    for (int i = 0; i < dim; ++i) dphiC_dx[i] = new double[6];

    double **ainvC_;
    ainvC_ = new double*[dim];
    for (int i = 0; i < dim; ++i) ainvC_[i] = new double[dim];

    double **JacC;
    JacC = new double*[dim];
    for (int i = 0; i < dim; ++i) JacC[i] = new double[dim];

   
    //Data for IGA fine mesh computation (Lagrange field)
    double wpc[9],xsi[dim], phi_[9];
    for (int i = 0; i < 9; i ++) wpc[i] = (*nodes_)[connect_[i]] -> getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]] -> getINC(); 

    double **dphi_dx;
    dphi_dx = new double*[dim];
    for (int i = 0; i < dim; ++i) dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double*[dim];
    for (int i = 0; i < dim; ++i) ainv_[i] = new double[dim];

    double **quadJacMat;
    quadJacMat = new double*[dim];
    for (int i = 0; i < dim; ++i) quadJacMat[i] = new double[dim];

    int index = 0;    
    for(double* it = sQuad.beginIso(); it != sQuad.endIso(); it++){
        
        if ((intPointCorrespElem_ISO[index] == ielem)){

            
            //Fine mesh computation
            //Defines the integration points adimentional coordinates
            for (int k = 0; k < dim; k ++) xsi[k] = sQuad.PointListIso(index,k);
            
            //Returns the quadrature integration weight
            weight_ = sQuad.WeightListIso(index);

            //Computes the velocity shape functions
            shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,(*iparameters),Npatch_);
        
            //Computes the jacobian matrix
            getJacobianMatrix_ISO(xsi,quadJacMat,ainv_);

            //Computes spatial derivatives
            getSpatialDerivatives_ISO(xsi,ainv_,dphi_dx);

            //Computes stabilization parameters
            getNewParameterSUPG_ISO(quadJacMat,phi_,dphi_dx);
            
            //Coarse mesh computatation
            //Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < dim; k++) xsiC[k] = intPointCorrespXsi_ISO[index][k];

            //Computes the velocity shape function
            shapeQuad.evaluateFem(xsiC,phiC_);

            //Computes the jacobian matrix
            getJacobianMatrix_COARSE_FEM(xsiC,JacC,ainvC_);

            //Computes spatial derivatives
            getSpatialDerivatives_FEM(xsiC, ainvC_, dphiC_dx);

            

            //Interpolates Lagrange multiplier and velocity
            getInterpolatedVariablesDifferentMesh_ISO_FEM(phi_,dphi_dx,phiC_,dphiC_dx);

            //Computes Matrix and vectors
            getMatrixAndVectorsDifferentMesh_tSUPG_tPSPG_ISO_FEM(phi_,dphiC_dx,
                                                                jacobianNRMatrix,rhsVector);


        };
        index++;        
    };  
     


    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] quadJacMat[i];
    delete [] quadJacMat;

    for (int i = 0; i < 2; ++i) delete [] dphiC_dx[i];
    delete [] dphiC_dx;
    for (int i = 0; i < 2; ++i) delete [] ainvC_[i];
    delete [] ainvC_;
    for (int i = 0; i < 2; ++i) delete [] JacC[i];
    delete [] JacC;  

    return;
};

template<>
void Element<2>::getLagrangeMultipliersDifferentMeshArlqStab_ISO_FEM(std::vector<Nodes *> &nodesCoarse_,int *connecC,int &ielem, 
                                                                    double **arlequinStabD, double *arlequinStabVectorD,
                                                                    double **arlequinStab0, double *arlequinStabVector0){
        int dim = 2;
    //quadrature and functions local classes
    SpecialQuad              sQuad = SpecialQuad(); 
    QuadShapeFunction<2>     shapeQuad;
    
    //Data for FEM coarse mesh computation (w function)
    nodesC_ = &nodesCoarse_;
    connectC_ = connecC;

    double xsiC[dim], phiC_[6];

    double **dphiC_dx;
    dphiC_dx = new double*[dim];
    for (int i = 0; i < dim; ++i) dphiC_dx[i] = new double[6];

    double **ainvC_;
    ainvC_ = new double*[dim];
    for (int i = 0; i < dim; ++i) ainvC_[i] = new double[dim];

    double **JacC;
    JacC = new double*[dim];
    for (int i = 0; i < dim; ++i) JacC[i] = new double[dim];

    double ***ddphiC_dx;
    ddphiC_dx = new double**[dim];
    for (int i = 0; i < dim; ++i) {
        ddphiC_dx[i] = new double*[2];
        for (int j = 0; j < dim; j++) ddphiC_dx[i][j] = new double[6];
    }; 

   
    //Data for IGA fine mesh computation (Lagrange field)
    double wpc[9],xsi[dim], phi_[9];
    for (int i = 0; i < 9; i ++) wpc[i] = (*nodes_)[connect_[i]] -> getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]] -> getINC(); 

    double **dphi_dx;
    dphi_dx = new double*[dim];
    for (int i = 0; i < dim; ++i) dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double*[dim];
    for (int i = 0; i < dim; ++i) ainv_[i] = new double[dim];

    double **quadJacMat;
    quadJacMat = new double*[dim];
    for (int i = 0; i < dim; ++i) quadJacMat[i] = new double[dim];

    int index = 0;    
    for(double* it = sQuad.beginIso(); it != sQuad.endIso(); it++){
        
        if ((intPointCorrespElem_ISO[index] == ielem)){

            
            //Fine mesh computation
            //Defines the integration points adimentional coordinates
            for (int k = 0; k < dim; k ++) xsi[k] = sQuad.PointListIso(index,k);
            
            //Returns the quadrature integration weight
            weight_ = sQuad.WeightListIso(index);

            //Computes the velocity shape functions
            shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,(*iparameters),Npatch_);
        
            //Computes the jacobian matrix
            getJacobianMatrix_ISO(xsi,quadJacMat,ainv_);

            //Computes spatial derivatives
            getSpatialDerivatives_ISO(xsi,ainv_,dphi_dx);

            //get Arlequin stabilization parameter
            getParameterArlequin_ISO(phi_,dphi_dx);

            
            //Coarse mesh computatation
            //Defines the equivalent integration point in coarse mesh
            for (int k = 0; k < dim; k++) xsiC[k] = intPointCorrespXsi_ISO[index][k];

            //Computes the velocity shape function
            shapeQuad.evaluateFem(xsiC,phiC_);

            //Computes the jacobian matrix
            getJacobianMatrix_COARSE_FEM(xsiC,JacC,ainvC_);

            //Computes spatial derivatives
            getSpatialDerivatives_FEM(xsiC, ainvC_, dphiC_dx);

            //Computes second spatial derivatives
            getSecondSpatialDerivatives_FEM(xsiC,ainvC_,ddphiC_dx);


            //Interpolates Lagrange multiplier and velocity
            getInterpolatedVariablesDifferentMeshArlqStab_ISO_FEM(dphi_dx,phiC_,dphiC_dx,ddphiC_dx);

            //Computes Matrix and vectors
            getMatrixAndVectorsDifferentMeshArlqStab_ISO_FEM(index,dphi_dx,phiC_,dphiC_dx,ddphiC_dx,
                                                             arlequinStabD,arlequinStabVectorD,
                                                             arlequinStab0,arlequinStabVector0);

        };
        index++;        
    };  


    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] quadJacMat[i];
    delete [] quadJacMat;

    for (int i = 0; i < 2; ++i) delete [] dphiC_dx[i];
    delete [] dphiC_dx;
    for (int i = 0; i < 2; ++i) delete [] ainvC_[i];
    delete [] ainvC_;
    for (int i = 0; i < 2; ++i) delete [] JacC[i];
    delete [] JacC;  

    return;
};



template<>
void Element<2>::getLagrangeMultipliersDifferentMesh_ISO_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_,int *connecC,
												 	        std::vector<IParameters *> &iparamC, int &ielem,
												 	        double **lagrMultMatrix, double *rhsVector1, double *rhsVector2,
                                                            double **arlequinStab, double *arlequinStabVector){
    int dim = 2;
    //quadrature and functions local classes
    SpecialQuad              sQuad = SpecialQuad(); 
    QuadShapeFunction<2>     shapeQuad;
    
    //Data for IGA coarse mesh computation (w function)
	nodesC_ = &nodesCoarse_;
    connectC_ = connecC;
    NpatchC_ = ipatchC;
    iparametersC = &iparamC;
    double wpcC[9],xsiC[dim], phiC_[9];
    for (int i = 0; i < 9; i ++) wpcC[i] = (*nodesC_)[connectC_[i]] -> getWeightPC();
    int *incC_ = (*nodesC_)[connectC_[8]] -> getINC(); 

    double **dphiC_dx;
    dphiC_dx = new double*[dim];
    for (int i = 0; i < dim; ++i) dphiC_dx[i] = new double[9];

    double **ainvC_;
    ainvC_ = new double*[dim];
    for (int i = 0; i < dim; ++i) ainvC_[i] = new double[dim];

    double **quadJacMatC;
    quadJacMatC = new double*[dim];
    for (int i = 0; i < dim; ++i) quadJacMatC[i] = new double[dim];

   
   //Data for IGA fine mesh computation (Lagrange field)
    double wpc[9],xsi[dim], phi_[9];
    for (int i = 0; i < 9; i ++) wpc[i] = (*nodes_)[connect_[i]] -> getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]] -> getINC();  

    double **dphi_dx;
    dphi_dx = new double*[dim];
    for (int i = 0; i < dim; ++i) dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double*[dim];
    for (int i = 0; i < dim; ++i) ainv_[i] = new double[dim];

    double **quadJacMat;
    quadJacMat = new double*[dim];
    for (int i = 0; i < dim; ++i) quadJacMat[i] = new double[dim];

    int index = 0;    
    for(double* it = sQuad.beginIso(); it != sQuad.endIso(); it++){
        
        if ((intPointCorrespElem_ISO[index] == ielem)){

          	//Fine mesh computation
          	//Defines the integration points adimentional coordinates
        	for (int k = 0; k < dim; k ++) xsi[k] = sQuad.PointListIso(index,k);
        	
        	//Returns the quadrature integration weight
        	weight_ = sQuad.WeightListIso(index);

        	//Computes the velocity shape functions
        	shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,(*iparameters),Npatch_);
        
        	//Computes the jacobian matrix
        	getJacobianMatrix_ISO(xsi,quadJacMat,ainv_);

        	//Computes spatial derivatives
        	getSpatialDerivatives_ISO(xsi,ainv_,dphi_dx);

            //Arlequin Stabilization parameter
            getParameterArlequin_ISO(phi_,dphi_dx);

        	
        	//Coarse mesh computatation
        	//Defines the equivalent integration point in coarse mesh
        	for (int k = 0; k < dim; k++) xsiC[k] = intPointCorrespXsi_ISO[index][k];

        	//Computes the velocity shape functions
        	shapeQuad.evaluateIso(xsiC,phiC_,wpcC,incC_,(*iparametersC),NpatchC_);

        	//Computes the jacobian matrix
        	getJacobianMatrix_COARSE_ISO(xsiC,quadJacMatC,ainvC_);

        	//Computes spatial derivatives
        	getSpatialDerivatives_COARSE_ISO(xsiC,ainvC_,dphiC_dx);

        	

            //Interpolates Lagrange multiplier and velocity
        	getInterpolatedVariablesDifferentMesh_ISO_ISO(phi_,dphi_dx,phiC_,dphiC_dx);

        	//Computes Matrix and vectors
        	getMatrixAndVectorsDifferentMesh_ISO_ISO(phi_,dphi_dx,phiC_,dphiC_dx,
        										    lagrMultMatrix,rhsVector1,rhsVector2,
                                                     arlequinStab,arlequinStabVector);


        };
        index++;        
    };  

    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] quadJacMat[i];
    delete [] quadJacMat;

    for (int i = 0; i < 2; ++i) delete [] dphiC_dx[i];
    delete [] dphiC_dx;
    for (int i = 0; i < 2; ++i) delete [] ainvC_[i];
    delete [] ainvC_;
    for (int i = 0; i < 2; ++i) delete [] quadJacMatC[i];
    delete [] quadJacMatC;
    

    return;

};
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//----------------MOUNTS EACH TYPE OF INCOMPRESSIBLE FLOW PROBEM----------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
//------------------------------------------------------------------------------
template<>
void Element<2>::setBoundaryConditions_ISO(double **jacobianNRMatrix,double *rhsVector){


    for (int i = 0; i < 9; i++){
        //direction x
        int constrain = (*nodes_)[connect_[i]] -> getConstrains(0);
        if ((constrain == 1) || (constrain == 3)){
            for (int j = 0; j < 27; j++){
                jacobianNRMatrix[2*i][j] = 0.0;
                jacobianNRMatrix[j][2*i] = 0.0;
            }
            jacobianNRMatrix[2*i][2*i] = 1.0;
            rhsVector[2*i] = 0.0;
        }   

        //direction y
        constrain = (*nodes_)[connect_[i]] -> getConstrains(1);
        if ((constrain == 1) || (constrain == 3)){
            for (int j = 0; j < 27; j++){
                jacobianNRMatrix[2*i+1][j] = 0.0;
                jacobianNRMatrix[j][2*i+1] = 0.0;
            }
            jacobianNRMatrix[2*i+1][2*i+1] = 1.0;
            rhsVector[2*i+1] = 0.0;
        }

        //if PRESSURE CONDITION ADD NODE!
    }
}
//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------

template<>
void Element<2>::getTransientNavierStokes_FEM(double **jacobianNRMatrix,double *rhsVector){

    
    //quadrature and functions local classes
    NormalQuad              nQuad = NormalQuad(); 
    QuadShapeFunction<2>    shapeQuad;

    int index = 0;

    for(double* it = nQuad.beginFem(); it != nQuad.endFem(); it++){

    	//variables
        double phi_[6],xsi[2];

        double **dphi_dx;
        dphi_dx = new double*[2];
        for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[6];

        double **ainv_;
        ainv_ = new double*[2];
        for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

        double **Jac;
        Jac = new double*[2];
        for (int i = 0; i < 2; ++i) Jac[i] = new double[2];

        //Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointListFem(index,0);
        xsi[1] = nQuad.PointListFem(index,1);
        //Returns the quadrature integration weight
        weight_ = nQuad.WeightListFem(index);

        //Computes the velocity shape functions
        shapeQuad.evaluateFem(xsi,phi_);

        //Computes the jacobian matrix
        getJacobianMatrix_FEM(xsi,Jac,ainv_);

        //Computes spatial derivatives
        getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

        //Interpolates variables and its derivatives values
        getVelAndDerivatives_FEM(phi_,dphi_dx);

        //Compute Stabilization Parameters
        getNewParameterSUPG_FEM(Jac,phi_,dphi_dx);

        //Computes the element matrix
        getElemMatrix_FEM(index,phi_,dphi_dx,jacobianNRMatrix);

        //Computes the RHS vector
        getResidualVector_FEM(index,phi_,dphi_dx,rhsVector); 

        
        for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
        delete [] dphi_dx;
        for (int i = 0; i < 2; ++i) delete [] ainv_[i];
        delete [] ainv_;
        for (int i = 0; i < 2; ++i) delete [] Jac[i];
        delete [] Jac;

        
        index++;       
    };

    
};


template<>
void Element<2>::getTransientNavierStokes_ISO(double **jacobianNRMatrix,double *rhsVector){
  
    
    //quadrature and functions local classes
    NormalQuad              nQuad = NormalQuad(); 
    QuadShapeFunction<2>    shapeQuad;
    
    //data for IGA elements
    double wpc[9],xsi[2], phi_[9];
    for (int i = 0; i < 9; i ++) wpc[i] = (*nodes_)[connect_[i]] -> getWeightPC();
    int *inc_ = (*nodes_)[connect_[8]] -> getINC();

    double **dphi_dx;
    dphi_dx = new double*[2];
    for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[9];

    double **ainv_;
    ainv_ = new double*[2];
    for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    double **quadJacMat;
    quadJacMat = new double*[2];
    for (int i = 0; i < 2; ++i) quadJacMat[i] = new double[2];  

    int index = 0;
    for(double* it = nQuad.beginIso(); it != nQuad.endIso(); it++){

        //Defines the integration points adimentional coordinates
        xsi[0] = nQuad.PointListIso(index,0);
        xsi[1] = nQuad.PointListIso(index,1);
        //Returns the quadrature integration weight
        weight_ = nQuad.WeightListIso(index);

        //Computes the velocity shape functions
        shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,(*iparameters),Npatch_);
        
        //Computes the jacobian matrix
        getJacobianMatrix_ISO(xsi,quadJacMat,ainv_);

        //Computes spatial derivatives
        getSpatialDerivatives_ISO(xsi,ainv_,dphi_dx);

        //Interpolates variables and pressure and its derivatives values
        getVelAndDerivatives_ISO(phi_,dphi_dx);

        //Compute Stabilization Parameters
        getNewParameterSUPG_ISO(quadJacMat,phi_,dphi_dx);

        //Computes the element matrix
        getElemMatrix_ISO(index,phi_,dphi_dx,jacobianNRMatrix);

        //Computes the RHS vector
        getResidualVector_ISO(index,phi_,dphi_dx,rhsVector); 
       
        index++; 
    };


    for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    delete [] dphi_dx;
    for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    delete [] ainv_;
    for (int i = 0; i < 2; ++i) delete [] quadJacMat[i];
    delete [] quadJacMat;
   
};






#endif

