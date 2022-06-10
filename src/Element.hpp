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
	FParameters                parameters;          // Fluid parameters
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
    double         weight_;							// Weight for numerical integration
    double         djac_;                           // Jacobian determinant
    double         x_, y_, z_;	           			// Interpolated coordinates
    double         u_, v_, w_, p_;         			// Interpolated velocity and pressure
    double         uprev_, vprev_, wprev_; 			// Interpolated previous velocity 
    double         ax_, ay_, az_;          			// Interpolated acceleration
    double         axprev_, ayprev_, azprev_;       // Interpolated previous acceleration
	double         umesh_, vmesh_, wmesh_; 			// Interpolated mesh velocity
	double         umeshprev_, vmeshprev_,
				   wmeshprev_; 			   			// Interpolated previous mesh velocity
    double         du_dx, du_dy, du_dz, 
    			   dv_dx, dv_dy, dv_dz, 
                   dw_dx, dw_dy, dw_dz;    			// Interpolated fluid spatial derivatives
    double         duprev_dx, duprev_dy, duprev_dz,
                   dvprev_dx, dvprev_dy, dvprev_dz,
                   dwprev_dx, dwprev_dy, dwprev_dz; // Interpolated previous fluid spatial derivatives
    double         dp_dx, dp_dy, dp_dz;				// Interpolated pressure spatial derivative   
    double         lamx_, lamy_;					// Lagrange multipliers
    double         lamx_dx, lamx_dy, 			    // Lagrange multipliers derivatives
    			   lamy_dx, lamy_dy;      

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
    // double intPointCorrespXsi[25][DIM];		       
    // double intPointCorrespElem[25];
    // double intPointCorrespXsi[64][DIM];         
    // double intPointCorrespElem[64];  
    double intPointCorrespXsi[9][DIM];            
    double intPointCorrespElem[9]; 
    
    // Integration point energy weight 
    // double intPointWeightFunction_FEM[8*DIM-9];
    double intPointWeightFunction_ISO[18*DIM-27];
    // double intPointWeightFunctionPrev_FEM[8*DIM-9];
    double intPointWeightFunctionPrev_ISO[18*DIM-27];
    
    // Integration point energy weight gluing zone
    double intPointWeightFunctionSpecial[9];
    double intPointWeightFunctionSpecialPrev[9];
    // double intPointWeightFunctionSpecial[25];
    // double intPointWeightFunctionSpecialPrev[25];
    // double intPointWeightFunctionSpecial[64];
    // double intPointWeightFunctionSpecialPrev[64];
    
public:
    // fluid element constructor
    Element(int index, int *connect, std::vector<Nodes *> &nodes, int ElemType, FParameters &param, std::vector<IParameters *> &iparam, int ipatch){
        index_ = index;
        connect_ = connect;
        nodes_ = &nodes;
        ElemType_ = ElemType;
        parameters = param;
        iparameters = &iparam;
        Npatch_ = ipatch;
        glueZone = false;

        // if (ElemType == 0){
        // 	for (int i = 0; i < 8*DIM-9; i++){
        // 		intPointWeightFunction_FEM[i] = 1.;
        // 		intPointWeightFunctionPrev_FEM[i] = 1.;
        // 	};
        // } else { //IGA elements
        	for (int i = 0; i < 18*DIM-27; i++){
        		intPointWeightFunction_ISO[i] = 1.;
        		intPointWeightFunctionPrev_ISO[i] = 1.;
        	};
        // };

        for (int i = 0; i < 9; i++){
        	intPointWeightFunctionSpecial[i] = 1.;
        	intPointWeightFunctionSpecialPrev[i] = 1.;
           
            for (int j = 0; j < DIM; j++) intPointCorrespXsi[i][j] = 0.0;
            intPointCorrespElem[i] = 0;

        };

        // for (int i = 0; i < 64; i++){
        // 	intPointWeightFunctionSpecial = 1.;
        // 	intPointWeightFunctionSpecialPrev = 1.;
        // }

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
    // int getNumberOfIntegrationPoints_FEM(){NormalQuad nQuad = NormalQuad(); return (nQuad.end() - nQuad.begin());}; //outside of the gluing zone for FEM elements
    int getNumberOfIntegrationPoints_ISO(){NormalQuad nQuad = NormalQuad(); return (nQuad.endIso() - nQuad.beginIso());}; //outside of the gluing zone for IGA elements
    int getNumberOfIntegrationPointsSpecial(){SpecialQuad sQuad = SpecialQuad(); return (sQuad.end() - sQuad.begin());}; //inside of gluing zone
   
    // Set the integration point and element correspondence in the coarse mesh
    void setIntegrationPointCorrespondence(int ipoint, double *x, int elem){
        intPointCorrespElem[ipoint] = elem;
        intPointCorrespXsi[ipoint][0] = x[0];
        intPointCorrespXsi[ipoint][1] = x[1]; 
    };
    
    // Returns correspondent element in the coarse mesh
    int getIntegPointCorrespondenceElement(int index){
    	return intPointCorrespElem[index];
    };

    // Returns the correspondent integration point coordinates in the coarse mesh
    double getIntegPointCoordinatesValue(int index, int dir){
        double x = intPointCorrespXsi[index][dir];
        return x;
    };

    // Compute the weighted integration points
    void setIntegPointWeightFunction_ISO();
    // void setIntegPointWeightFunction_FEM();

    //......................Jacobian Matrixes and Derivatives....................
    // Compute and store the spatial jacobian matrix
    // void getJacobianMatrix_FEM(double *xsi, double **Jac, double **ainv_);
    void getJacobianMatrix_ISO(double *xsi, double **quadJacMat, double **ainv_);
    // Compute and store the quadrature jacobian matrix
    void getQuadJacobianMatrix_ISO(double *xsi, double **quadJacMatInv);

 //    // gets Jacobian Matrix Values
 //    void getJacobianMatrixValues_FEM(double *xsi, double **ainv_){

 //        double **Jac;
 //        Jac= new double*[DIM];
 //        for (int i = 0; i < DIM; ++i) Jac[i] = new double[DIM];

 //        getJacobianMatrix_FEM(xsi,Jac,ainv_);

 //    	for (int i = 0; i < DIM; ++i) delete [] Jac[i];
 //    	delete [] Jac;
	// };

    void getJacobianMatrixValues_ISO(double *xsi, double **ainv_ ){
        
        double **quadJacMat;
        quadJacMat = new double*[DIM];
        for (int i = 0; i < DIM; ++i) quadJacMat[i] = new double[DIM];
        
        getJacobianMatrix_ISO(xsi,quadJacMat,ainv_);

    	for (int i = 0; i < DIM; ++i) delete [] quadJacMat[i];
    	delete [] quadJacMat;
	}

	// gets Jacobian Matrix determinant value
    // double getDetJac_FEM(double *xsi){
               
    //     double **ainv_;
    //     ainv_ = new double*[DIM];
    //     for (int i = 0; i < DIM; ++i) ainv_[i] = new double[DIM];

    //     double **Jac;
    //     Jac= new double*[DIM];
    //     for (int i = 0; i < DIM; ++i) Jac[i] = new double[DIM];
        
    //     getJacobianMatrix_FEM(xsi,Jac,ainv_);

    // 	for (int i = 0; i < DIM; ++i) delete [] ainv_[i];
    // 	delete [] ainv_;

    // 	for (int i = 0; i < DIM; ++i) delete [] Jac[i];
    // 	delete [] Jac;
        
    //     return djac_;
    // }

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
    // void getSpatialDerivatives_FEM(double *xsi, double **ainv_, double **dphi_dx);
    void getSpatialDerivatives_ISO(double *xsi, double **ainv_, double **phi_dx);

    // Compute and stores the interpolated variables
    // void getVelAndDerivatives_FEM(double *phi_, double **dphi_dx);
    void getVelAndDerivatives_ISO(double *phi_, double **dphi_dx);

    // Compute and stores the interpolated variables when Arlequin Problem
    void getInterpolatedVariablesSameMesh_ISO(double *phi_, double **dphi_dx);
    void getInterpolatedVariablesDifferentMesh_ISO(double *phi_, double **dphi_dx,double *phiC_, double **dphiC_dx);


    //......................Stabilization Parameters....................
    // Compute and store the SUPG, PSPG and LSIC stabilization parameters
    // void getNewParameterSUPG_FEM(double **Jac , double *phi, double **dphi_dx);
    void getNewParameterSUPG_ISO(double **quadJacMat, double *phi_, double **dphi_dx);
 

 	//......................Drag and Lift Parameters....................
    // Compute and store the drag and lift forces at the element boundary
    void computeDragAndLiftForces_ISO(int bound);
    // Gets the element pressure drag force
    double getPressureDragForce(){return pressureDragForce;};
    // Gets the element pressure lift force
    double getPressureLiftForce(){return pressureLiftForce;};
    // Gets the element friction drag force
    double getFrictionDragForce(){return frictionDragForce;};
    // Gets the element friction lift force
    double getFrictionLiftForce(){return frictionLiftForce;};
    // Gets the element drag force
    double getDragForce(){return dragForce;};
    // Gets the element lift force
    double getLiftForce(){return liftForce;};
    // Gets the element Pitching Moment
    double getPitchingMoment(){return pitchingMoment;};
    // Gets the element boundary perimeter
    double getPerimeter(){return perimeter;};
   
   
    //.......................Element vectors and matrices.......................
    // Compute and store the element matrix for the incompressible flow problem
    // void getElemMatrix_FEM(double *phi_, double **dphi_dx, double **jacobianNRMatrix);
    void getElemMatrix_ISO(int &index,double *phi_, double **dphi_dx, double **jacobianNRMatrix);   
    void getMatrixAndVectorsSameMesh_ISO(double *phi_, double **dphi_dx, double **lagrMultMatrix, 
        	                             double *rhsVector1, double *rhsVector2); 
    void getMatrixAndVectorsDifferentMesh_ISO(double *phi_, double **dphi_dx, double *phiC_, double **dphiC_dx,
        									  double **lagrMultMatrix, double *rhsVector1, double *rhsVector2);

    //Compute and store the residual vector for the incompressible flow problem
    // void getResidualVector_FEM(double *phi_, double **dphi_dx, double *rhsVector);
    void getResidualVector_ISO(int &index,double *phi_, double **dphi_dx, double *rhsVector);

    //.................Element vectors and matrices Arlequin.....................
    // Compute and store the Lagrange multiplier operator when integrating the same mesh portion
    void getLagrangeMultipliersSameMesh(double **lagrMultMatrix, double *lagrMultVector, double *rhsVector);

    //Compute and store the Lagrange multiplier operator when integrationg the different mesh portion
    void getLagrangeMultipliersDifferentMesh_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_,int *connecC,
												 std::vector<IParameters *> &iparamC, int &ielem, 
												 double **lagrMultMatrix, double *rhsVectorLM, double *rhsVector);

    //...............................Problem type...............................
    // Compute the Transient Navier-Stokes problem matrices and vectors
    // void getTransientNavierStokes_FEM(double **jacobianNRMatrix, double *rhsVector);
    void getTransientNavierStokes_ISO(double **jacobianNRMatrix, double *rhsVector);


    //.......................Bezier Element transformation.......................
    // Computes Bézier extract operator 
    void getMatrixC(double **MatrixC);
    //Computes inverse Bézier extractor operator    
    void getInvMatrixC(double **MatrixCuInv, double **MatrixCvInv);

    double getIntegPointWeightFunction_ISO(int ind) {return intPointWeightFunctionSpecial[ind];};

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------WEIGHTED INTEGRATION POINTS-------------------------
//------------------------------------------------------------------------------
// template<>
// void Element<2>::setIntegPointWeightFunction_FEM(){


    // NormalQuad              nQuad = NormalQuad(); 
    // QuadShapeFunction<2>    shapeQuad;

    // int index = 0;
    // for(double* it = nQuad.begin(); it != nQuad.end(); it++){
    // 	intPointWeightFunctionPrev_FEM[index] = intPointWeightFunction_FEM[index];
    // 	intPointWeightFunction_FEM[index] = 0.;
    // 	index++
    // };
   
    // index = 0;
    // for(double* it = nQuad.begin(); it != nQuad.end(); it++){

    // 	//variables
    //     double phi_[6],xsi[2];
    //     xsi[0] = nQuad.PointList(index,0);
    //     xsi[1] = nQuad.PointList(index,1);
    //     shapeQuad.evaluate(xsi,phi_);

    //     for (int i = 0; i < 6; i++){
    //     	intPointWeightFunction_FEM[index] += (*nodes_)[connect_[i]] -> getWeightFunction() * phi_[i];
    //     }
    	
    // 	index++;
    // };
    
// };

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
    SpecialQuad sQuad = SpecialQuad();
    index = 0;
    for(double* it = sQuad.begin(); it != sQuad.end(); it++){
    	intPointWeightFunctionSpecialPrev[index] = intPointWeightFunctionSpecial[index];
    	intPointWeightFunctionSpecial[index] = 0.;
    	index++;
    };

    index = 0;
    for(double* it = sQuad.begin(); it != sQuad.end(); it++){
    	
    	for (int i = 0; i < dim; i++) xsi[i] = sQuad.PointList(index,i);
        shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,(*iparameters),Npatch_);

        for (int i = 0; i < 9; i++){
        	intPointWeightFunctionSpecial[index] += (*nodes_)[connect_[i]] -> getWeightFunction() * phi_[i];
        };
    	
    	index++;
    };


};




//------------------------------------------------------------------------------
//-------------------------SPATIAL TRANSFORM - JACOBIAN-------------------------
//------------------------------------------------------------------------------
// template<>
// void Element<2>::getJacobianMatrix_FEM(double *xsi,  double **Jac,  double **ainv_) {

//     //Computes the spatial Jacobian matrix and its inverse
//     QuadShapeFunction<2>    shapeQuad;
    
//     double **dphi;
//     dphi = new double*[2];
//     for (int i = 0; i < 2; ++i) dphi[i] = new double[6];
    
//     shapeQuad.evaluateGradient(xsi,dphi);

//     double dx_dxsi1 = 0.;
//     double dx_dxsi2 = 0.;
//     double dy_dxsi1 = 0.;
//     double dy_dxsi2 = 0.;

//     double &alpha_f = parameters.getAlphaF();

//     for (int i=0; i<6; i++){
		
// 		double xna_ = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(0) + 
//                       (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(0);
//         double yna_ = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(1) + 
//                       (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(1);    

//         dx_dxsi1 += xna_ * dphi[0][i];
//         dx_dxsi2 += xna_ * dphi[1][i];
//         dy_dxsi1 += yna_ * dphi[0][i];
//         dy_dxsi2 += yna_ * dphi[1][i];   
//     };

//     // Defining the Jacobian matrix
//     Jac[0][0] = dx_dxsi1;
//     Jac[0][1] = dx_dxsi2;
//     Jac[1][0] = dy_dxsi1;
//     Jac[1][1] = dy_dxsi2;

//     //Computing the jacobian determinant
//     djac_ = Jac[0][0] * Jac[1][1] - Jac[0][1] * Jac[1][0];

//     //Computing Jacobian inverse (transposed)
//     ainv_[0][0] =  dy_dxsi2 / djac_;
//     ainv_[0][1] = -dy_dxsi1 / djac_;
//     ainv_[1][0] = -dx_dxsi2 / djac_;
//     ainv_[1][1] =  dx_dxsi1 / djac_;

//     for (int i = 0; i < 2; ++i) delete [] dphi[i];
//     delete [] dphi;

//     return;
    
// };


template<>
void Element<2>::getJacobianMatrix_ISO(double *xsi, double **quadJacMat, double **ainv_) {

    QuadShapeFunction<2>    shapeQuad;

    double &alpha_f = parameters.getAlphaF();
   
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
void Element<2>::getQuadJacobianMatrix_ISO(double *xsi, double **quadJacMatInv){

    QuadShapeFunction<2>    shapeQuad;

    double &alpha_f = parameters.getAlphaF();
    
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
// template<>
// void Element<2>::getSpatialDerivatives_FEM(double *xsi, double **ainv_,double **dphi_dx) {
    
//     QuadShapeFunction<2>    shapeQuad;

//     double **dphi;
//     dphi = new double*[2];
//     for (int i = 0; i < 2; ++i) dphi[i] = new double[6];
    
//     shapeQuad.evaluateGradient(xsi,dphi);
    
//     for (int i = 0; i < 2; ++i)
//         for (int j = 0; j < 6; ++j)
//             dphi_dx[i][j] = 0.;

//     //Quadratic shape functions spatial first derivatives
//     for (int i = 0; i < 2; i++)
//         for (int j = 0; j < 2; j++)
//             for (int k = 0; k < 6; k++)
//                 dphi_dx[i][k] += ainv_[i][j]*dphi[j][k];

//     for (int i = 0; i < 2; ++i) delete [] dphi[i];
//     delete [] dphi;

//     return;
// };


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


//------------------------------------------------------------------------------
//-------------INTERPOLATES VELOCITY, PRESSURE AND ITS DERIVATIVES--------------
//------------------------------------------------------------------------------
// template<>
// void Element<2>::getVelAndDerivatives_FEM(double *phi_, double **dphi_dx) {

// 	x_ = 0.;          y_ = 0.;
//     u_ = 0.;          v_ = 0.;
//     du_dx = 0.;       du_dy = 0.;    
//     dv_dx = 0.;       dv_dy = 0.;
//     uprev_ = 0.;      vprev_ = 0.;
//     duprev_dx = 0.;   duprev_dy = 0.;
//     dvprev_dx = 0.;   dvprev_dy = 0.;
//     ax_ = 0.;         ay_ = 0.;
//     axprev_ = 0.;     ayprev_ = 0.;
//     umesh_ = 0.;      vmesh_ = 0.;
//     umeshprev_ = 0.;  vmeshprev_ = 0.;
//     p_ = 0.;
//     dp_dx = 0.;       dp_dy = 0.;
      
//     double &alpha_f = parameters.getAlphaF();

//     //Interpolates the velocity components and its spatial derivatives
//     for (int i = 0; i < 6; i++){
        
//         //coordinates
//         double xna_ = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(0) + 
//                       (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(0);
//         double yna_ = alpha_f * (*nodes_)[connect_[i]] -> getCoordinateValue(1) + 
//                       (1. - alpha_f) * (*nodes_)[connect_[i]] -> getPreviousCoordinateValue(1);
        
//         x_ += xna_ * phi_[i];
//         y_ += yna_ * phi_[i];

//         //velocity
//         u_ += (*nodes_)[connect_[i]] -> getVelocity(0) * phi_[i];
//         v_ += (*nodes_)[connect_[i]] -> getVelocity(1) * phi_[i];

//         //previous velocity
//         uprev_ += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * phi_[i];
//         vprev_ += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * phi_[i];

//         //velocity derivatives
//         du_dx += (*nodes_)[connect_[i]] -> getVelocity(0) * dphi_dx[0][i];
//         du_dy += (*nodes_)[connect_[i]] -> getVelocity(0) * dphi_dx[1][i];
//         dv_dx += (*nodes_)[connect_[i]] -> getVelocity(1) * dphi_dx[0][i];
//         dv_dy += (*nodes_)[connect_[i]] -> getVelocity(1) * dphi_dx[1][i];

//         //previous velocity derivatives
//         duprev_dx += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * dphi_dx[0][i];
//         duprev_dy += (*nodes_)[connect_[i]] -> getPreviousVelocity(0) * dphi_dx[1][i];
//         dvprev_dx += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * dphi_dx[0][i];
//         dvprev_dy += (*nodes_)[connect_[i]] -> getPreviousVelocity(1) * dphi_dx[1][i];

//         //acceleration
//         ax_ += (*nodes_)[connect_[i]] -> getAcceleration(0) * phi_[i];
//         ay_ += (*nodes_)[connect_[i]] -> getAcceleration(1) * phi_[i];
        
//         //previous acceleration
//         axprev_ += (*nodes_)[connect_[i]] -> getPreviousAcceleration(0) * phi_[i];
//         ayprev_ += (*nodes_)[connect_[i]] -> getPreviousAcceleration(1) * phi_[i];

//         //mesh velocity
//         umesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * phi_[i];
//         vmesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * phi_[i];
        
//         //previous mesh velocity
//         umeshprev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * phi_[i];
//         vmeshprev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * phi_[i];

//         //pressure
//         p_ += (*nodes_)[connect_[i]] -> getPressure() * phi_[i];

//         //pressure derivatives
//         dp_dx += (*nodes_)[connect_[i]] -> getPressure() * dphi_dx[0][i];
//         dp_dy += (*nodes_)[connect_[i]] -> getPressure() * dphi_dx[1][i];
//     };  

//     return;
// };

template<>
void Element<2>::getVelAndDerivatives_ISO(double *phi_, double **dphi_dx) {

	x_ = 0.;          y_ = 0.;
    u_ = 0.;          v_ = 0.;
    du_dx = 0.;       du_dy = 0.;    
    dv_dx = 0.;       dv_dy = 0.;
    uprev_ = 0.;      vprev_ = 0.;
    duprev_dx = 0.;   duprev_dy = 0.;
    dvprev_dx = 0.;   dvprev_dy = 0.;
    ax_ = 0.;         ay_ = 0.;
    axprev_ = 0.;     ayprev_ = 0.;
    umesh_ = 0.;      vmesh_ = 0.;
    umeshprev_ = 0.;  vmeshprev_ = 0.;
    p_ = 0.;
    dp_dx = 0.;       dp_dy = 0.;
    
    double &alpha_f = parameters.getAlphaF();

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

        //mesh velocity
        umesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(0) * phi_[i];
        vmesh_ += (*nodes_)[connect_[i]] -> getMeshVelocity(1) * phi_[i];

        //previous mesh velocity
        umeshprev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(0) * phi_[i];
        vmeshprev_ += (*nodes_)[connect_[i]] -> getPreviousMeshVelocity(1) * phi_[i];

        //pressure
        p_ += (*nodes_)[connect_[i]] -> getPressure() * phi_[i];

        //Pressure derivatives
        dp_dx += (*nodes_)[connect_[i]] -> getPressure() * dphi_dx[0][i];
        dp_dy += (*nodes_)[connect_[i]] -> getPressure() * dphi_dx[1][i];

    };  

    return;
};

template<>
void Element<2>::getInterpolatedVariablesSameMesh_ISO(double *phi_, double **dphi_dx) {

    u_ = 0.;          v_ = 0.;
    du_dx = 0.;       du_dy = 0.;    
    dv_dx = 0.;       dv_dy = 0.;
    uprev_ = 0.;      vprev_ = 0.;
    duprev_dx = 0.;   duprev_dy = 0.;
    dvprev_dx = 0.;   dvprev_dy = 0.;
    lamx_ = 0.;       lamy_ = 0.;	
    lamx_dx = 0.;     lamx_dy = 0.;
    lamy_dx = 0.;     lamy_dy = 0.;      
        
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

    };  

    return;
};


template<>
void Element<2>::getInterpolatedVariablesDifferentMesh_ISO(double *phi_, double **dphi_dx,double *phiC_, double **dphiC_dx) {

	u_ = 0.;          v_ = 0.;
    du_dx = 0.;       du_dy = 0.;    
    dv_dx = 0.;       dv_dy = 0.;
    uprev_ = 0.;      vprev_ = 0.;
    duprev_dx = 0.;   duprev_dy = 0.;
    dvprev_dx = 0.;   dvprev_dy = 0.;
    lamx_ = 0.;       lamy_ = 0.;	
    lamx_dx = 0.;     lamx_dy = 0.;
    lamy_dx = 0.;     lamy_dy = 0.;  
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
void Element<2>::computeDragAndLiftForces_ISO(int bound) {

    // double    &visc_ = parameters.getViscosity();
    // double    localNodesBoundaryIso_[3][2];
    // double    pcWeightl[3];
    // double    pcWeight[9];
    // int       *incL_;
   
    // for (int i = 0; i<9; i++) pcWeight[i] = (*nodes_)[connect_[i]] -> getWeightPC();

    // int size = elemBound_.size();

    // for (int iS=0; iS< size;iS++){

    //     if (elemBound_[iS] == bound){

    //         if(elemSide_[iS] == 0){

    //             pcWeightl[0] = (*nodes_)[connect_[0]] -> getWeightPC();
    //             pcWeightl[1] = (*nodes_)[connect_[1]] -> getWeightPC();
    //             pcWeightl[2] = (*nodes_)[connect_[2]] -> getWeightPC();
                
    //             for (int i=0; i<2; i++){
    //                 localNodesBoundaryIso_[0][i] = ((*nodes_)[connect_[0]] -> getCoordinateValue(i))/pcWeightl[0];
    //                 localNodesBoundaryIso_[1][i] = ((*nodes_)[connect_[1]] -> getCoordinateValue(i))/pcWeightl[1];
    //                 localNodesBoundaryIso_[2][i] = ((*nodes_)[connect_[2]] -> getCoordinateValue(i))/pcWeightl[2];           
    //             };

    //             incL_ = (*nodes_)[connect_[2]] -> getINC(); 

    //         }else{

    //             if(elemSide_[iS] == 1){

    //                 pcWeightl[0] = (*nodes_)[connect_[2]] -> getWeightPC();
    //                 pcWeightl[1] = (*nodes_)[connect_[5]] -> getWeightPC();
    //                 pcWeightl[2] = (*nodes_)[connect_[8]] -> getWeightPC();

    //                 for (int i=0; i<2; i++){
    //                     localNodesBoundaryIso_[0][i] = ((*nodes_)[connect_[2]] -> getCoordinateValue(i))/pcWeightl[0];
    //                     localNodesBoundaryIso_[1][i] = ((*nodes_)[connect_[5]] -> getCoordinateValue(i))/pcWeightl[1];
    //                     localNodesBoundaryIso_[2][i] = ((*nodes_)[connect_[8]] -> getCoordinateValue(i))/pcWeightl[2];
    //                 };

    //                 incL_ = (*nodes_)[connect_[8]] -> getINC(); 

    //             if(elemSide_[iS] == 2){

    //                 pcWeightl[0] = (*nodes_)[connect_[8]] -> getWeightPC();
    //                 pcWeightl[1] = (*nodes_)[connect_[7]] -> getWeightPC();
    //                 pcWeightl[2] = (*nodes_)[connect_[6]] -> getWeightPC();
                        
    //                     for (int i=0; i<2; i++){
    //                     localNodesBoundaryIso_[0][i] = ((*nodes_)[connect_[8]] -> getCoordinateValue(i))/pcWeightl[0];
    //                     localNodesBoundaryIso_[1][i] = ((*nodes_)[connect_[7]] -> getCoordinateValue(i))/pcWeightl[1];
    //                     localNodesBoundaryIso_[2][i] = ((*nodes_)[connect_[6]] -> getCoordinateValue(i))/pcWeightl[2];
    //                 };

    //                 incL_ = (*nodes_)[connect_[8]] -> getINC(); 
    //             };

    //             }else{

    //                 pcWeightl[0] = (*nodes_)[connect_[6]] -> getWeightPC();
    //                 pcWeightl[1] = (*nodes_)[connect_[3]] -> getWeightPC();
    //                 pcWeightl[2] = (*nodes_)[connect_[0]] -> getWeightPC();
                    
    //                 for (int i=0; i<2; i++){
    //                     localNodesBoundaryIso_[0][i] = ((*nodes_)[connect_[6]] -> getCoordinateValue(i))/pcWeightl[0];
    //                     localNodesBoundaryIso_[1][i] = ((*nodes_)[connect_[3]] -> getCoordinateValue(i))/pcWeightl[1];
    //                     localNodesBoundaryIso_[2][i] = ((*nodes_)[connect_[0]] -> getCoordinateValue(i))/pcWeightl[2];
    //                 };

    //                 incL_ = (*nodes_)[connect_[6]] -> getINC(); 
    //             };        
    //         };



    //         BoundaryIntegQuadratureIso<2>                      bQuad;    
    //         QuadShapeFunction<2>                           shapeQuad; 
    //         BoundShapeFunction<2>                          shapeBound;

    //         double n_vector[2] = {};
    //         double shearStress[2][2] = {};
    //         double ident[2][2];
    //         double load_friction[2] = {};
    //         double load_pressure[2] = {};

    //         ident[0][0] = 1.;
    //         ident[0][1] = 0.;
    //         ident[1][0] = 0.;
    //         ident[1][1] = 1.;

    //         double *phi_;
    //         phi_ = new double[9];

    //         double *phiCB_;
    //         phiCB_ = new double[3];

    //         double **dphiCB_;
    //         dphiCB_ = new double*[1];
    //         for (int i = 0; i < 1; ++i) dphiCB_[i] = new double[3];

    //         double **dphi_dx;
    //         dphi_dx = new double*[2];
    //         for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[9];

    //         double **ainv_;
    //         ainv_ = new double*[2];
    //         for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

    //         double **quadJacMat;
    //         quadJacMat = new double*[2];
    //         for (int i = 0; i < 2; ++i) quadJacMat[i] = new double[2];

    //         double xsi[2] = {};
    //         double xsiB[1] = {};
    //         double weightB[1] = {};

    //         double moment = 0.;
    //         double per = 0.;
    //         int index = 0;
            
    //         // indexes the base fuction (index space) that starts in the left lower corner of the element
    //         int *inc_;
    //         inc_ = (*nodes_)[connect_[8]] -> getINC();  

    //         for(double* it = bQuad.begin(); it != bQuad.end(); it++){
                
    //             xsiB[0] = bQuad.PointList(index,0);
    //             weightB[0] = bQuad.WeightList(index);

    //             if(elemSide_[iS] == 0){
    //                 xsi[0] = xsiB[0];
    //                 xsi[1] = -1;
    //             };
    //             if(elemSide_[iS] == 1){
    //                 xsi[0] = -1;
    //                 xsi[1] = xsiB[0];
    //             };
    //             if(elemSide_[iS] == 2){
    //                 xsi[0] = xsiB[0];
    //                 xsi[1] = 1;
    //             };
    //             if(elemSide_[iS] == 3){
    //                 xsi[0] = -1;
    //                 xsi[1] = xsiB[0];
    //             };

    //             //Computes the shape functions
    //             shapeQuad.evaluateIso(xsi,phi_,pcWeight,inc_,(*iparameters),Npatch_);
                
    //             //Computes the jacobian matrix
    //             getJacobianMatrixIso(xsi,quadJacMatInv,quadJacMat,ainv_);

    //             //Computes spatial derivatives
    //             getSpatialDerivativesIso(xsi,ainv_,dphi_dx);

    //             //Interpolates velocity and its derivatives values
    //             getVelAndDerivativesIso(phi_,dphi_dx);   

    //             //Get Functions and Derivatives on the boundary
    //             shapeBound.evaluateBoundaryIso(xsiB,phiCB_,pcWeightl,elemSide_[iS],incL_,(*iparameters),Npatch_);
    //             shapeBound.evaluateGradientBoundaryIso(xsiB,dphiCB_,pcWeightl,elemSide_[iS],incL_,(*iparameters),Npatch_);
               
    //             double Tx=0.; double Ty = 0.;
                
    //             // (Jacobian parametric space to physical one)
    //             for (int i=0; i<3; i++){
    //                 Tx += localNodesBoundaryIso_[i][0] * dphiCB_[0][i];
    //                 Ty += localNodesBoundaryIso_[i][1] * dphiCB_[0][i];
    //             };

    //             // (Jacobian parental space to the physical one)
    //             if ((elemSide_[iS] == 0) || (elemSide_[iS] == 2)){
                    
    //                 int degm = (*iparameters)[Npatch_] -> getDegree(0);
    //                 int npcm = (*iparameters)[Npatch_] -> getNcp(0);  
    //                 double *uknot_;
    //                 uknot_ = (*iparameters)[Npatch_] -> getuKnot();

    //                 int uind = incL_[0];

    //                 double u1 = uknot_[uind];
    //                 double u2 = uknot_[uind+1];

    //                 double dxsi1_dqxsi1 = (u2-u1)*0.5;
                    
    //                 Tx *= dxsi1_dqxsi1;
    //                 Ty *= dxsi1_dqxsi1;
         
    //             } else {

    //                 int degn = (*iparameters)[Npatch_] -> getDegree(1);
    //                 int npcn = (*iparameters)[Npatch_] -> getNcp(1);   
    //                 double *vknot_;
    //                 vknot_ = (*iparameters)[Npatch_] -> getvKnot();  

    //                 int vind = incL_[1]; 

    //                 double v1 = vknot_[vind];
    //                 double v2 = vknot_[vind +1];

    //                 double dxsi2_dqxsi2 = (v2-v1)*0.5;

    //                 Tx *= dxsi2_dqxsi2;
    //                 Ty *= dxsi2_dqxsi2;
    //             }
                
    //             double jacb_ = sqrt(Tx*Tx + Ty*Ty);
                
    //             n_vector[0] =  Ty / jacb_;
    //             n_vector[1] = -Tx / jacb_;

    //             shearStress[0][0] = 2. * visc_ * du_dx;
    //             shearStress[0][1] = visc_ * (du_dy + dv_dx);
    //             shearStress[1][0] = visc_ * (du_dy + dv_dx);
    //             shearStress[1][1] = 2. * visc_ * dv_dy;

    //             for (int i = 0; i < 2; i++){
    //                 for (int j = 0; j < 2; j++){
    //                     load_pressure[i] += -p_*ident[i][j]*n_vector[j] * jacb_ * weightB[0];
    //                     load_friction[i] += shearStress[i][j]*n_vector[j] * jacb_ * weightB[0];
    //                 }
    //             }

    //             moment +=  (-(load_pressure[0]+load_friction[0]) * (y_- 4.) 
    //                        + (load_pressure[0]+load_friction[0]) * (x_ - 4.)) * jacb_ * weightB[0];
    //             per += jacb_ * weightB[0];
                
    //             index++;
    //         };

    //         perimeter = per;
    //         pitchingMoment = moment;

    //         pressureDragForce = -load_pressure[0];
    //         pressureLiftForce = -load_pressure[1];
            
    //         frictionDragForce = -load_friction[0];
    //         frictionLiftForce = -load_friction[1]; 

    //         dragForce = pressureDragForce + frictionDragForce;
    //         liftForce = pressureLiftForce + frictionLiftForce;


    //         delete [] phi_;
    //         delete [] phiCB_;
    //         for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
    //         delete [] dphi_dx;
    //         for (int i = 0; i < 1; ++i) delete [] dphiCB_[i];
    //         delete [] dphiCB_;   
    //         for (int i = 0; i < 2; ++i) delete [] ainv_[i];
    //         delete [] ainv_;
    //         for (int i = 0; i < 2; ++i) delete [] quadJacMat[i];
    //         delete [] quadJacMat;
    //     };
    // };


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
// template<>
// void Element<2>::getNewParameterSUPG_FEM(double **Jac, double *phi_, double **dphi_dx) {
    
//     double MatrixD[2][2],MatrixInvD[2][2],MatrixQh[2][2],MatrixInvQh[2][2],MatrixG[2][2],rr[2][2];   
//     double r[2] = {}; double ua_[2] = {}; 
//     double hrqd,tSUGN1_,tSUGN2_,tSUGN3_;

//     double &dTime_ = parameters.getTimeStep();
//     double &visc_ = parameters.getViscosity();
//     double &dens_ = parameters.getDensity();

//     //Matrix D = polynomial order of the base functions
//     MatrixD[0][0] = 2.;
//     MatrixD[1][1] = 2.;
//     MatrixD[1][0] = 0.;
//     MatrixD[0][1] = 0.;

//     //Inverse Matrix D
//     double detMatrixD = MatrixD[0][0] * MatrixD[1][1] - MatrixD[0][1] * MatrixD[1][0];
//     MatrixInvD[0][0] = (1./detMatrixD) * MatrixD[1][1];
//     MatrixInvD[1][1] = (1./detMatrixD) * MatrixD[0][0];
//     MatrixInvD[0][1] = -(1./detMatrixD) * MatrixD[0][1];
//     MatrixInvD[1][0] = -(1./detMatrixD) * MatrixD[1][0];

//     // Matrix Q "hat"
//     MatrixQh[0][0] = Jac[0][0] * MatrixInvD[0][0] + Jac[0][1] * MatrixInvD[1][0];
//     MatrixQh[0][1] = Jac[0][0] * MatrixInvD[0][1] + Jac[0][1] * MatrixInvD[1][1];
//     MatrixQh[1][0] = Jac[1][0] * MatrixInvD[0][0] + Jac[1][1] * MatrixInvD[1][0];
//     MatrixQh[1][1] = Jac[1][0] * MatrixInvD[0][1] + Jac[1][1] * MatrixInvD[1][1];

//     //Matrix Inverse Q "hat"
//     double detQh = MatrixQh[0][0] * MatrixQh[1][1] - MatrixQh[0][1] * MatrixQh[1][0];

//     MatrixInvQh[0][0] = (1./detQh) *  MatrixQh[1][1];
//     MatrixInvQh[1][1] = (1./detQh) *  MatrixQh[0][0];
//     MatrixInvQh[0][1] = -(1./detQh) * MatrixQh[0][1];
//     MatrixInvQh[1][0] = -(1./detQh) * MatrixQh[1][0];

//     // Matrix G
//     MatrixG[0][0] = MatrixInvQh[0][0] * MatrixInvQh[0][0] + MatrixInvQh[1][0] * MatrixInvQh[1][0];
//     MatrixG[0][1] = MatrixInvQh[0][0] * MatrixInvQh[0][1] + MatrixInvQh[1][0] * MatrixInvQh[1][1];
//     MatrixG[1][0] = MatrixInvQh[0][1] * MatrixInvQh[0][0] + MatrixInvQh[1][1] * MatrixInvQh[1][0];
//     MatrixG[1][1] = MatrixInvQh[0][1] * MatrixInvQh[0][1] + MatrixInvQh[1][1] * MatrixInvQh[1][1];

//     // Calculation hrqd
//     for (int i = 0; i < 6; i++){

//         double ua = (*nodes_)[connect_[i]] -> getVelocity(0);
//         double va = (*nodes_)[connect_[i]] -> getVelocity(1);

//         double uma = (*nodes_)[connect_[i]] -> getMeshVelocity(0);
//         double vma = (*nodes_)[connect_[i]] -> getMeshVelocity(1);

//         ua -= uma;
//         va -= vma;

//         ua_[0] += ua * phi_[i];
//         ua_[1] += va * phi_[i];
        
//         r[0] += sqrt(ua * ua + va * va) * dphi_dx[0][i];
//         r[1] += sqrt(ua * ua + va * va) * dphi_dx[1][i];
//    };

//     double rNorm = sqrt(r[0]*r[0] + r[1]*r[1]);;
//     double referpar = 0.000000000000000001;

//     //physical unitary direction from gradient velocity 
//     r[0] /= (rNorm + referpar);
//     r[1] /= (rNorm + referpar);
   
//     rr[0][0] = r[0] * r[0];
//     rr[0][1] = r[0] * r[1];
//     rr[1][0] = r[1] * r[0];
//     rr[1][1] = r[1] * r[1];

//     hrqd = 2./((sqrt(rr[0][0] * MatrixG[0][0] + rr[0][1] * MatrixG[0][1] + rr[1][0] * MatrixG[1][0] + rr[1][1] * MatrixG[1][1])) + referpar);

//     //hmin e hmax
//     //calculating the eigen-values from G
//     double a_ = 1.;
//     double b_ = -(MatrixG[0][0] + MatrixG[1][1]);
//     double c_ = (MatrixG[0][0] * MatrixG[1][1] - MatrixG[0][1] * MatrixG[1][0]);

//     double lambdaMax = (- b_ + sqrt(b_ * b_ - 4. * a_ * c_))/ (2. * a_);
//     double lambdaMin = (- b_ - sqrt(b_ * b_ - 4. * a_ * c_))/ (2. * a_);

//     double hmin = 2. / sqrt(lambdaMax);
//     double hmax = 2. / sqrt(lambdaMin);

//     if (hrqd < hmin) hrqd = hmin;
//     if (hrqd > hmax) hrqd = hmax;

//     //tSUGN1_ (-2)
//     tSUGN1_ = (ua_[0] * ua_[0]) * MatrixG[0][0] +
//               (ua_[0] * ua_[1]) * MatrixG[0][1] + 
//               (ua_[1] * ua_[0]) * MatrixG[1][0] + 
//               (ua_[1] * ua_[1]) * MatrixG[1][1];
    
//     tSUGN2_ = dTime_/2.;

//     //tSUGN3_ (-1)
//     tSUGN3_ = (visc_/dens_)*(4./(hrqd*hrqd));

//     //Computing tSUPG parameter
//     tSUPG_ = 1./(sqrt(tSUGN1_ + (1./(tSUGN2_ * tSUGN2_)) + (tSUGN3_ * tSUGN3_)));
   
//     tLSIC_ = hrqd * hrqd / tSUPG_;
    
//     return;

// };

template<>
void Element<2>::getNewParameterSUPG_ISO(double **quadJacMat, double *phi_, double **dphi_dx) {

    double MatrixD[2][2] = {};
    double MatrixInvD[2][2],MatrixQh[2][2],MatrixInvQh[2][2],MatrixG[2][2],rr[2][2];   
    double r[2] = {}; double ua_[2] = {}; 
    double hrqd,tSUGN1_,tSUGN2_,tSUGN3_;
    double &dTime_ = parameters.getTimeStep();
    double &visc_ = parameters.getViscosity();
    double &dens_ = parameters.getDensity(); 

    
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
   
    tLSIC_ = hrqd * hrqd / tSUPG_;

    for (int i = 0; i < 3; ++i) delete [] MatrixCuInv[i];
    delete [] MatrixCuInv;
    for (int i = 0; i < 3; ++i) delete [] MatrixCvInv[i];
    delete [] MatrixCvInv;
    
    return;
};


//------------------------------------------------------------------------------
//-----------------------------ELEMENT LOCAL MATRIX-----------------------------
//------------------------------------------------------------------------------

// template<>
// void Element<2>::getElemMatrix_FEM(double *phi_, double **dphi_dx,double **jacobianNRMatrix){
    
//     double &dTime_ = parameters.getTimeStep();
//     double &visc_ = parameters.getViscosity();
//     double &dens_ = parameters.getDensity();
//     double &alpha_f = parameters.getAlphaF();
//     double &alpha_m = parameters.getAlphaM();
//     double &gamma = parameters.getGamma();

//     double una_ = alpha_f * u_ + (1. - alpha_f) * uprev_;
//     double vna_ = alpha_f * v_ + (1. - alpha_f) * vprev_;

//     double axm_ = alpha_m * ax_ + (1. - alpha_m) * axprev_;
//     double aym_ = alpha_m * ay_ + (1. - alpha_m) * ayprev_;

//     double duna_dx = alpha_f * du_dx + (1. - alpha_f) * duprev_dx;
//     double duna_dy = alpha_f * du_dy + (1. - alpha_f) * duprev_dy;
//     double dvna_dx = alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx;
//     double dvna_dy = alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy;

//     double umeshna_ = alpha_f * umesh_ + (1. - alpha_f) * umeshprev_;
//     double vmeshna_ = alpha_f * vmesh_ + (1. - alpha_f) * vmeshprev_;


//     for (int i = 0; i < 6; i++){
//         for (int j = 0; j < 6; j++){

//             double wSUPGi = (una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i];
//             double wSUPGj = (una_ - umeshna_) * dphi_dx[0][j] + (vna_ - vmeshna_) * dphi_dx[1][j];
            
//             //Mass matrix (use for both directions) 
//             double mM =  phi_[i] * phi_[j] * dens_ * alpha_m + 
//                          wSUPGi * phi_[j] * tSUPG_ * dens_ * alpha_m;
                                 
//             //Difusion matrix (viscosity)
//             double Kxx = (2. * dphi_dx[0][i] * dphi_dx[0][j] + dphi_dx[1][i] * dphi_dx[1][j]) 
//                          * visc_* alpha_f * gamma * dTime_;
//             double Kxy = dphi_dx[1][i] * dphi_dx[0][j] * visc_* alpha_f * gamma * dTime_;
//             double Kyx = dphi_dx[0][i] * dphi_dx[1][j] * visc_ * alpha_f * gamma * dTime_;
//             double Kyy = (2. * dphi_dx[1][i] * dphi_dx[1][j] + dphi_dx[0][i] * dphi_dx[0][j])
//                          * visc_* alpha_f * gamma * dTime_;
            
//             //Convection matrixes
//             double Cxx = (dphi_dx[0][j] * (una_ - umeshna_) + dphi_dx[1][j] * (vna_ - vmeshna_)) 
//                         * phi_[i] * dens_ * alpha_f * gamma * dTime_
//                         + wSUPGi * wSUPGj * tSUPG_ * dens_* alpha_f * gamma * dTime_;
//             double Cyy  = Cxx;
            
//             double Cuu = phi_[i] * duna_dx * phi_[j] * dens_* alpha_f * gamma * dTime_;
//             double Cuv = phi_[i] * duna_dy * phi_[j] * dens_ * alpha_f * gamma * dTime_;
//             double Cvu = phi_[i] * dvna_dx * phi_[j] * dens_ * alpha_f * gamma * dTime_;
//             double Cvv = phi_[i] * dvna_dy * phi_[j] * dens_ * alpha_f * gamma * dTime_;

//             //Stabilization LSIC matrix          
//             double KLSxx = dphi_dx[0][i] * dphi_dx[0][j] * tLSIC_ * dens_ * alpha_f * gamma * dTime_;
//             double KLSxy = dphi_dx[0][i] * dphi_dx[1][j] * tLSIC_ * dens_ * alpha_f * gamma * dTime_;
//             double KLSyx = dphi_dx[1][i] * dphi_dx[0][j] * tLSIC_ * dens_ * alpha_f * gamma * dTime_;
//             double KLSyy = dphi_dx[1][i] * dphi_dx[1][j] * tLSIC_ * dens_ * alpha_f * gamma * dTime_;

//             jacobianNRMatrix[2*i  ][2*j  ] += (mM + Kxx + KLSxx + Cxx + Cuu) * weight_ * djac_;
//             jacobianNRMatrix[2*i+1][2*j+1] += (mM + Kyy + KLSyy + Cyy + Cvv) * weight_ * djac_;
//             jacobianNRMatrix[2*i  ][2*j+1] += (Kxy + Cuv + KLSxy) * weight_ * djac_;
//             jacobianNRMatrix[2*i+1][2*j  ] += (Kyx + Cvu + KLSyx) * weight_ * djac_; 
            
//             //multipy pressure direction x and y
//             double QSUPGx = - (dphi_dx[0][i] * phi_[j]) + 
//                               wSUPGi * dphi_dx[0][j] * tSUPG_;
//             double QSUPGy = - (dphi_dx[1][i] * phi_[j]) +
//                                wSUPGi *dphi_dx[1][j] * tSUPG_;
            
//             //multiply velocity direction x and y
//             double Qx = dphi_dx[0][j] * phi_[i] * alpha_f * gamma * dTime_;
//             double Qy = dphi_dx[1][j] * phi_[i] * alpha_f * gamma * dTime_;                

//             jacobianNRMatrix[12+i][2*j  ] += Qx * weight_ * djac_;
//             jacobianNRMatrix[12+i][2*j+1] += Qy * weight_ * djac_;
//             jacobianNRMatrix[2*i  ][12+j] += QSUPGx * weight_ * djac_;
//             jacobianNRMatrix[2*i+1][12+j] += QSUPGy * weight_ * djac_;

//             //PSPG stabilization matrixes
//             double Hx = dphi_dx[0][i] * phi_[j] * tPSPG_ * alpha_m;
//             double Hy = dphi_dx[1][i] * phi_[j] * tPSPG_* alpha_m;
            
//             double Gx = dphi_dx[0][i] * wSUPGj * tPSPG_ * alpha_f * gamma * dTime_;
//             double Gy = dphi_dx[1][i] * wSUPGj * tPSPG_ * alpha_f * gamma * dTime_;

//             double Q = (dphi_dx[0][i] * dphi_dx[0][j] + 
//                         dphi_dx[1][i] * dphi_dx[1][j]) * tPSPG_ / (dens_);         

//             jacobianNRMatrix[12+i][2*j  ] += (Hx + Gx) * weight_ * djac_;
//             jacobianNRMatrix[12+i][2*j+1] += (Hy + Gy) * weight_ * djac_;
//             jacobianNRMatrix[12+i][12+j]  +=  Q * weight_ * djac_;              
//         };
//     };

//     return;
// };

template<>
void Element<2>::getElemMatrix_ISO(int &index, double *phi_, double **dphi_dx, double **jacobianNRMatrix){
    
    double &dTime_ = parameters.getTimeStep();
    double &visc_ = parameters.getViscosity();
    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &alpha_m = parameters.getAlphaM();
    double &gamma = parameters.getGamma();

    double una_ = alpha_f * u_ + (1. - alpha_f) * uprev_;
    double vna_ = alpha_f * v_ + (1. - alpha_f) * vprev_;

    double duna_dx = alpha_f * du_dx + (1. - alpha_f) * duprev_dx;
    double duna_dy = alpha_f * du_dy + (1. - alpha_f) * duprev_dy;
    double dvna_dx = alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx;
    double dvna_dy = alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy;

    double axm_ = alpha_m * ax_ + (1. - alpha_m) * axprev_;
    double aym_ = alpha_m * ay_ + (1. - alpha_m) * ayprev_;

    double umeshna_ = alpha_f * umesh_ + (1. - alpha_f) * umeshprev_;
    double vmeshna_ = alpha_f * vmesh_ + (1. - alpha_f) * umeshprev_;

    double wna_ = alpha_f * intPointWeightFunction_ISO[index] + (1. - alpha_f) * intPointWeightFunctionPrev_ISO[index];

    //Stokes Problem   
    for (int i = 0; i < 9; i++){
        for (int j = 0; j < 9; j++){
                   
            //Difusion matrix (viscosity)
            double Kxx = (2. * dphi_dx[0][i] * dphi_dx[0][j] + 
                         dphi_dx[1][i] * dphi_dx[1][j]) * visc_* alpha_f * gamma * dTime_;
            double Kxy = dphi_dx[1][i] * dphi_dx[0][j] * visc_ * alpha_f * gamma * dTime_;
            double Kyx = dphi_dx[0][i] * dphi_dx[1][j] * visc_ * alpha_f * gamma * dTime_;
            double Kyy = (2. * dphi_dx[1][i] * dphi_dx[1][j] + 
                         dphi_dx[0][i] * dphi_dx[0][j]) * visc_* alpha_f * gamma * dTime_;
            
            jacobianNRMatrix[2*i  ][2*j  ] += Kxx * weight_ * djac_ * wna_;
            jacobianNRMatrix[2*i+1][2*j+1] += Kyy * weight_ * djac_* wna_;

			//multiply pressure 
            double QSUPGx = - dphi_dx[0][i] * phi_[j];
            double QSUPGy = - dphi_dx[1][i] * phi_[j];

            //multiply velocity 
            double Qx = dphi_dx[0][j] * phi_[i] * alpha_f * gamma * dTime_;
            double Qy = dphi_dx[1][j] * phi_[i] * alpha_f * gamma * dTime_;

            jacobianNRMatrix[18+i][2*j  ] += Qx * weight_ * djac_* wna_;
            jacobianNRMatrix[18+i][2*j+1] += Qy * weight_ * djac_* wna_;
            jacobianNRMatrix[2*i  ][18+j] += QSUPGx * weight_ * djac_* wna_;
            jacobianNRMatrix[2*i+1][18+j] += QSUPGy * weight_ * djac_* wna_;

            double Q = (dphi_dx[0][i] * dphi_dx[0][j] + 
                        dphi_dx[1][i] * dphi_dx[1][j]) * tPSPG_ / (dens_);
           
            jacobianNRMatrix[18+i][18+j] +=  Q * weight_ * djac_* wna_;

        };
    };

    return;
};

template<>
void Element<2>::getMatrixAndVectorsSameMesh_ISO(double *phi_, double **dphi_dx, double **lagrMultMatrix, 
        	                                     double *rhsVector1, double *rhsVector2){

	//Fluid Data
    double &visc_ = parameters.getViscosity();
    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &alpha_m = parameters.getAlphaM();
    double &k1 = parameters.getArlequinK1();
    double &k2 = parameters.getArlequinK2();

    double una_ = alpha_f * u_ + (1. - alpha_f) * uprev_;
	double vna_ = alpha_f * v_ + (1. - alpha_f) * vprev_;

	double duna_dx = alpha_f * du_dx + (1. - alpha_f) * duprev_dx;
	double duna_dy = alpha_f * du_dy + (1. - alpha_f) * duprev_dy;
	double dvna_dx = alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx;
	double dvna_dy = alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy;


	for (int i = 0; i < 9; i++){
		for (int j = 0; j < 9; j++){

			//L2 operator
			double L2 = phi_[i] * phi_[j] * k1;

			lagrMultMatrix[2*i][2*j] -= L2 * djac_ * weight_;
			lagrMultMatrix[2*i+1][2*j+1] -= L2 * djac_ * weight_;

			// //H1 operator
			// double H1xx = (2. * dphi_dx[0][i] * dphi_dx[0][j] + 
   //                         dphi_dx[1][i] * dphi_dx[1][j]) * k2;
			// double H1xy = dphi_dx[1][i] * dphi_dx[0][j]* k2;
			// double H1yx = dphi_dx[0][i] * dphi_dx[1][j]* k2;
			// double H1yy = (2. * dphi_dx[1][i] * dphi_dx[1][j] + 
   //                         dphi_dx[0][i] * dphi_dx[0][j])* k2;

			// lagrMultMatrix[2*i][2*j] -= H1xx * djac_ * weight_;
			// lagrMultMatrix[2*i+1][2*j] -= H1yx * djac_ * weight_;
			// lagrMultMatrix[2*i][2*j+1] -= H1xy * djac_ * weight_;
			// lagrMultMatrix[2*i+1][2*j+1] -= H1yy * djac_ * weight_;


		};

		//L2 operator - Lagrange
		double l2x_ = phi_[i] * lamx_ * k1;
		double l2y_ = phi_[i] * lamy_ * k1;

		//H1 operator - Lagrange
		double h1x_ = 0.;//(2. * dphi_dx[0][i] * lamx_dx + 
                      // dphi_dx[1][i] * lamx_dy + 
                      // dphi_dx[1][i] * lamy_dx) * k2 ;
        
        double h1y_ = 0.;//(dphi_dx[0][i] * lamx_dy + 
                      // 2. * dphi_dx[1][i] * lamy_dy + 
                      // dphi_dx[0][i] * lamy_dx) * k2;  

		rhsVector1[2*i] += (l2x_+ h1x_)  * djac_ * weight_;
		rhsVector1[2*i+1] += (l2y_+ h1y_)  * djac_ * weight_;

		//L2 operator - Velocity
		double l2ux_ = phi_[i] * una_ * k1;
		double l2uy_ = phi_[i] * vna_ * k1;

		//H1 operator - Velocity
		double h1ux_ = 0.;//(2. * dphi_dx[0][i] * duna_dx+ 
                       // dphi_dx[1][i] * duna_dy + 
                       // dphi_dx[1][i] * dvna_dx) * k2;
        
        double h1uy_= 0.;//(dphi_dx[0][i] * duna_dy + 
                      // 2. * dphi_dx[1][i] * dvna_dy + 
                      // dphi_dx[0][i] * dvna_dx) * k2; 

		rhsVector2[2*i] += (l2ux_ + h1ux_) * djac_ * weight_;
		rhsVector2[2*i+1] += (l2uy_ + h1uy_) * djac_ * weight_;

	};

};

template<>
void Element<2>::getMatrixAndVectorsDifferentMesh_ISO(double *phi_, double **dphi_dx, double *phiC_, double **dphiC_dx,
        										      double **lagrMultMatrix, double *rhsVector1, double *rhsVector2){

	double &visc_ = parameters.getViscosity();
    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &alpha_m = parameters.getAlphaM();
    double &k1 = parameters.getArlequinK1();
    double &k2 = parameters.getArlequinK2();

    double una_ = alpha_f * u_ + (1. - alpha_f) * uprev_;
	double vna_ = alpha_f * v_ + (1. - alpha_f) * vprev_;

	double duna_dx = alpha_f * du_dx + (1. - alpha_f) * duprev_dx;
	double duna_dy = alpha_f * du_dy + (1. - alpha_f) * duprev_dy;
	double dvna_dx = alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx;
	double dvna_dy = alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy;
	

	for (int i = 0; i < 9; i++){
		for (int j = 0; j < 9; j++){

			//L2 operator
			double L2 = phi_[i] * phiC_[j] * k1;

			lagrMultMatrix[2*i][2*j] += L2 * djac_ * weight_;
			lagrMultMatrix[2*i+1][2*j+1] += L2 * djac_ * weight_;

			// //H1 operator
			// double H1xx = (2. * dphi_dx[0][i] * dphiC_dx[0][j] + 
   //                         dphi_dx[1][i] * dphiC_dx[1][j]) * k2;
			// double H1xy = dphi_dx[1][i] * dphiC_dx[0][j]* k2;
			// double H1yx = dphi_dx[0][i] * dphiC_dx[1][j]* k2;
			// double H1yy = (2. * dphi_dx[1][i] * dphiC_dx[1][j] + 
   //                         dphi_dx[0][i] * dphiC_dx[0][j])* k2;

			// lagrMultMatrix[2*i][2*j] += H1xx * djac_ * weight_;
			// lagrMultMatrix[2*i+1][2*j] += H1yx * djac_ * weight_;
			// lagrMultMatrix[2*i][2*j+1] += H1xy * djac_ * weight_;
			// lagrMultMatrix[2*i+1][2*j+1] += H1yy * djac_ * weight_;

		};

		//L2 operator - Lagrange
		double l2x_ = phiC_[i] * lamx_ * k1;
		double l2y_ = phiC_[i] * lamy_ * k1;

		//H1 operator - Lagrange
		double h1x_ = 0.;//(2. * dphiC_dx[0][i] * lamx_dx + 
                      // dphiC_dx[1][i] * lamx_dy + 
                      // dphiC_dx[1][i] * lamy_dx) * k2 ;
        
        double h1y_ = 0.;//(dphiC_dx[0][i] * lamx_dy + 
                      // 2. * dphiC_dx[1][i] * lamy_dy + 
                      // dphiC_dx[0][i] * lamy_dx) * k2;  

		rhsVector1[2*i] -= (l2x_+ h1x_)  * djac_ * weight_;
		rhsVector1[2*i+1] -= (l2y_+ h1y_)  * djac_ * weight_;

		//L2 operator - Velocity
		double l2ux_ = phi_[i] * una_ * k1;
		double l2uy_ = phi_[i] * vna_ * k1;

		//H1 operator - Velocity
		double h1ux_ = 0.;//(2. * dphi_dx[0][i] * duna_dx+ 
                       // dphi_dx[1][i] * duna_dy + 
                       // dphi_dx[1][i] * dvna_dx) * k2;
        
        double h1uy_= 0.;//(dphi_dx[0][i] * duna_dy + 
                      // 2. * dphi_dx[1][i] * dvna_dy + 
                      // dphi_dx[0][i] * dvna_dx) * k2; 

		rhsVector2[2*i] -= (l2ux_ + h1ux_) * djac_ * weight_;
		rhsVector2[2*i+1] -= (l2uy_ + h1uy_) * djac_ * weight_;

	};
};



//------------------------------------------------------------------------------
//-----------------------------RESIDUAL - RHS VECTOR----------------------------
//------------------------------------------------------------------------------
// template<>
// void Element<2>::getResidualVector_FEM(double *phi_, double **dphi_dx, double *rhsVector){
    
//     double &visc_ = parameters.getViscosity();
//     double &dens_ = parameters.getDensity();
//     double &alpha_f = parameters.getAlphaF();
//     double &alpha_m = parameters.getAlphaM();

//     double ff[2];
//     ff[0] = parameters.getFieldForce(0);
//     ff[1] = parameters.getFieldForce(1);

//     double una_ = alpha_f * u_ + (1. - alpha_f) * uprev_;
//     double vna_ = alpha_f * v_ + (1. - alpha_f) * vprev_;

//     double duna_dx = alpha_f * du_dx + (1. - alpha_f) * duprev_dx;
//     double duna_dy = alpha_f * du_dy + (1. - alpha_f) * duprev_dy;
//     double dvna_dx = alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx;
//     double dvna_dy = alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy;

//     double axm_ = alpha_m * ax_ + (1. - alpha_m) * axprev_;
//     double aym_ = alpha_m * ay_ + (1. - alpha_m) * ayprev_;

//     double umeshna_ = alpha_f * umesh_ + (1. - alpha_f) * umeshprev_;
//     double vmeshna_ = alpha_f * vmesh_ + (1. - alpha_f) * vmeshprev_;

    
//     for (int i = 0; i < 6; i++){
       
//         double mx = phi_[i] * (axm_) * dens_ + 
//                     ((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i]) * 
//                     axm_ * tSUPG_ * dens_;
//         double my = phi_[i] * (aym_) * dens_ + 
//                     ((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i]) * 
//                     aym_ * tSUPG_ * dens_;

//         double Kx = 2. * dphi_dx[0][i] * duna_dx * visc_ + 
//                     dphi_dx[1][i] * duna_dy * visc_ + 
//                     dphi_dx[1][i] * dvna_dx * visc_ ;
        
//         double Ky= dphi_dx[0][i] * duna_dy * visc_ + 
//                    2. * dphi_dx[1][i] * dvna_dy * visc_ + 
//                    dphi_dx[0][i] * dvna_dx * visc_;            

//         double KLSx = dphi_dx[0][i] * (duna_dx + dvna_dy) * tLSIC_ * dens_;
//         double KLSy = dphi_dx[1][i] * (duna_dx + dvna_dy) * tLSIC_ * dens_;

//         double Cx = (duna_dx * (una_ - umeshna_) + 
//                     duna_dy * (vna_ - vmeshna_)) * phi_[i] * dens_ +
//                     ((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i]) *
//                     ((una_ - umeshna_) * duna_dx + (vna_ - vmeshna_) * duna_dy) * tSUPG_ *dens_;
//         double Cy = (dvna_dx * (una_ - umeshna_) + 
//                     dvna_dy * (vna_ - vmeshna_)) * phi_[i] * dens_ +
//                     ((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i]) *
//                     ((una_ - umeshna_) * dvna_dx + (vna_ - vmeshna_) * dvna_dy) * tSUPG_ *dens_;
 
//         double Px = - (dphi_dx[0][i] * p_) + 
//                       ((dphi_dx[0][i] * (una_ - umeshna_) + dphi_dx[1][i] * (vna_ - vmeshna_)) * dp_dx * tSUPG_);
//         double Py = - (dphi_dx[1][i] * p_) + 
//                       ((dphi_dx[0][i] * (una_ - umeshna_) + dphi_dx[1][i] * (vna_ - vmeshna_)) * dp_dy * tSUPG_);
           
//         double Q = ((duna_dx + dvna_dy) * phi_[i]) +
//                     (dphi_dx[0][i] * dp_dx + dphi_dx[1][i] * dp_dy) * tPSPG_ / dens_ +
//                      dphi_dx[0][i] * ((una_ - umeshna_) * duna_dx +
//                                     (vna_ - vmeshna_) * duna_dy) * tPSPG_ +
//                     dphi_dx[1][i] * ((una_ - umeshna_) * dvna_dx +
//                                     (vna_ - vmeshna_) * dvna_dy) * tPSPG_ +
//                     dphi_dx[0][i] * axm_* tPSPG_ +
//                     dphi_dx[1][i] * aym_ * tPSPG_;
                   
//         double Ffvx = phi_[i]* dens_ * ff[0] + 
//                       tSUPG_* dens_ * ((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i]) * ff[0];
//         double Ffvy = phi_[i]* dens_ * ff[1] + 
//                       tSUPG_* dens_ * ((una_ - umeshna_) * dphi_dx[0][i] + (vna_ - vmeshna_) * dphi_dx[1][i]) * ff[1];

//         double Ffp = tPSPG_ * (dphi_dx[0][i] * ff[0] + dphi_dx[1][i] * ff[1]);
        

//         rhsVector[2*i  ] += (-mx + Ffvx -Kx - Px - Cx - KLSx) * weight_ * djac_;
//         rhsVector[2*i+1] += (-my + Ffvy -Ky - Py - Cy - KLSy) * weight_ * djac_;
//         rhsVector[12+i] +=  (Ffp -Q) * weight_ * djac_;

//     };

//     return;
// };

template<>
void Element<2>::getResidualVector_ISO(int &index,double *phi_, double **dphi_dx, double *rhsVector){
    
    double &visc_ = parameters.getViscosity();
    double &dens_ = parameters.getDensity();
    double &alpha_f = parameters.getAlphaF();
    double &alpha_m = parameters.getAlphaM();

    double ff[2];
    ff[0] = parameters.getFieldForce(0);
    ff[1] = parameters.getFieldForce(1);

    double duna_dx = alpha_f * du_dx + (1. - alpha_f) * duprev_dx;
    double duna_dy = alpha_f * du_dy + (1. - alpha_f) * duprev_dy;
    double dvna_dx = alpha_f * dv_dx + (1. - alpha_f) * dvprev_dx;
    double dvna_dy = alpha_f * dv_dy + (1. - alpha_f) * dvprev_dy;

    double wna_ = alpha_f * intPointWeightFunction_ISO[index] + (1. - alpha_f) * intPointWeightFunctionPrev_ISO[index];

    //Stokes problem
    for (int i = 0; i < 9; i++){

        double Kx = 2. * dphi_dx[0][i] * duna_dx * visc_ + 
                    dphi_dx[1][i] * duna_dy * visc_ + 
                    dphi_dx[1][i] * dvna_dx * visc_;
        
        double Ky = dphi_dx[0][i] * duna_dy * visc_ + 
                   2. * dphi_dx[1][i] * dvna_dy * visc_ + 
                   dphi_dx[0][i] * dvna_dx * visc_;                
  
        double Px = - (dphi_dx[0][i] * p_);
        double Py = - (dphi_dx[1][i] * p_);
           
        double Q = ((duna_dx + dvna_dy) * phi_[i]) +
                    (dphi_dx[0][i] * dp_dx + dphi_dx[1][i] * dp_dy) * tPSPG_ / dens_;
             
        double Ffvx = phi_[i]* dens_ * ff[0];
        double Ffvy = phi_[i]* dens_ * ff[1];

        double Ffp = tPSPG_ * (dphi_dx[0][i] * ff[0] + dphi_dx[1][i] * ff[1]);


        rhsVector[2*i  ] += (Ffvx - Kx - Px) * weight_ * djac_ * wna_;
        rhsVector[2*i+1] += (Ffvy - Ky - Py) * weight_ * djac_ * wna_;
       
        rhsVector[18+i] +=  (-Q + Ffp) * weight_ * djac_ * wna_;

    };
   
    return;
};

//------------------------------------------------------------------------------
//----------------------ARLEQUIN ELEMENT LOCAL MATRIX/VECTOR--------------------
//------------------------------------------------------------------------------

template<>
void Element<2>::getLagrangeMultipliersSameMesh(double **lagrMultMatrix, double *rhsVector1, double *rhsVector2){

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
void Element<2>::getLagrangeMultipliersDifferentMesh_ISO(int &ipatchC, std::vector<Nodes *> &nodesCoarse_,int *connecC,
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
    for(double* it = sQuad.begin(); it != sQuad.end(); it++){
        
        if ((intPointCorrespElem[index] == ielem)){

          	//Fine mesh computation
          	//Defines the integration points adimentional coordinates
        	for (int k = 0; k < dim; k ++) xsi[k] = sQuad.PointList(index,k);
        	
        	//Returns the quadrature integration weight
        	weight_ = sQuad.WeightList(index);

        	//Computes the velocity shape functions
        	shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,(*iparameters),Npatch_);
        
        	//Computes the jacobian matrix
        	getJacobianMatrix_ISO(xsi,quadJacMat,ainv_);

        	//Computes spatial derivatives
        	getSpatialDerivatives_ISO(xsi,ainv_,dphi_dx);

        	
        	//Coarse mesh computatation
        	//Defines the equivalent integration point in coarse mesh
        	for (int k = 0; k < dim; k++) xsiC[k] = intPointCorrespXsi[index][k];

        	//Computes the velocity shape functions
        	shapeQuad.evaluateIso(xsiC,phiC_,wpcC,incC_,(*iparametersC),NpatchC_);

        	//Computes the jacobian matrix
        	getJacobianMatrix_ISO(xsiC,quadJacMatC,ainvC_);

        	//Computes spatial derivatives
        	getSpatialDerivatives_ISO(xsiC,ainvC_,dphiC_dx);


        	//Interpolates Lagrange multiplier and velocity
        	getInterpolatedVariablesDifferentMesh_ISO(phi_,dphi_dx,phiC_,dphiC_dx);

        	//Computes Matrix and vectors
        	getMatrixAndVectorsDifferentMesh_ISO(phi_,dphi_dx,phiC_,dphiC_dx,
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

//------------------------------------------------------------------------------
//-----------------------TRANSIENT NAVIER-STOKES PROBEM-------------------------
//------------------------------------------------------------------------------
// template<>
// void Element<2>::getTransientNavierStokes_FEM(double **jacobianNRMatrix,double *rhsVector){

//     //quadrature and functions local classes
//     NormalQuad              nQuad = NormalQuad(); 
//     QuadShapeFunction<2>    shapeQuad;

//     int index = 0;

//     for(double* it = nQuad.begin(); it != nQuad.end(); it++){

//     	//variables
//         double phi_[6],xsi[2];

//         double **dphi_dx;
//         dphi_dx = new double*[2];
//         for (int i = 0; i < 2; ++i) dphi_dx[i] = new double[6];

//         double **ainv_;
//         ainv_ = new double*[2];
//         for (int i = 0; i < 2; ++i) ainv_[i] = new double[2];

//         double **Jac;
//         Jac = new double*[2];
//         for (int i = 0; i < 2; ++i) Jac[i] = new double[2];

//         //Defines the integration points adimentional coordinates
//         xsi[0] = nQuad.PointList(index,0);
//         xsi[1] = nQuad.PointList(index,1);
//         //Returns the quadrature integration weight
//         weight_ = nQuad.WeightList(index);

//         //Computes the velocity shape functions
//         shapeQuad.evaluate(xsi,phi_);

//         //Computes the jacobian matrix
//         getJacobianMatrix_FEM(xsi,Jac,ainv_);

//         //Computes spatial derivatives
//         getSpatialDerivatives_FEM(xsi,ainv_,dphi_dx);

//         //Interpolates variables and its derivatives values
//         getVelAndDerivatives_FEM(phi_,dphi_dx);

//         //Compute Stabilization Parameters
//         getNewParameterSUPG_FEM(Jac,phi_,dphi_dx);
//         tPSPG_ = tSUPG_;

//         //Computes the element matrix
//         getElemMatrix_FEM(phi_,dphi_dx,jacobianNRMatrix);

//         //Computes the RHS vector
//         getResidualVector_FEM(phi_,dphi_dx,rhsVector); 

        
//         for (int i = 0; i < 2; ++i) delete [] dphi_dx[i];
//         delete [] dphi_dx;
//         for (int i = 0; i < 2; ++i) delete [] ainv_[i];
//         delete [] ainv_;
//         for (int i = 0; i < 2; ++i) delete [] Jac[i];
//         delete [] Jac;

        
//         index++;       
//     };

    
// };


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
        tPSPG_ = tSUPG_;

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

