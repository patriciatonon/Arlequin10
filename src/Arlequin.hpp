//------------------------------------------------------------------------------
// 
//           Jeferson W D Fernandes, Rodolfo A K Sanches and Patricia Tonon
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//-----------------------------------------------------------------------------

//------------------------------------------------------------------------------
//----------------------------------ARLEQUIN-------------------------------------
//------------------------------------------------------------------------------

#ifndef ARLEQUIN_H
#define ARLEQUIN_H
#include "FluidData.hpp"
#include "Glue.hpp"

template<int DIM>
class Arlequin{

public:
	//Defines locally the classes:
    // Class of fluid data
    typedef FluidData<DIM>                      FluidMesh;
    // Class of Element
    typedef typename FluidMesh::Elements        Element;
    // Class of Node 
    typedef typename FluidMesh::Node            Nodes;
    // Class of bondary
    typedef typename FluidMesh::Boundaries      Boundary;
    // Class of quadrature
    typedef typename Element::NormalQuad        Quadrature;
    // Class of special quadrature
    typedef typename Element::SpecialQuad       SpecialQuadrature;
    // Class of isogeometric parameters
    typedef typename FluidMesh::IsoParameters   IsoParameters;
    // Class of fluid parameters
    typedef FluidParameters<DIM>                Parameters;
    // Class glue
    typedef Glue<DIM> 					        GlueZone;


    //Defines the objects from the classes:
    FluidMesh coarseModel, fineModel;

    Parameters *parametersCoarse,*parametersFine;

    std::vector<Nodes *>     nodesCoarse_;
    std::vector<Nodes *>     nodesFine_;
    std::vector<Nodes *>     nodesGlueFine_;
    std::vector<Nodes *>     nodesGlueCoarse_;
    std::vector<Nodes *>     nodesLagrangeFine_;

    std::vector<Element *>   elementsCoarse_;
    std::vector<Element *>   elementsFine_;

    std::vector<IsoParameters* > IsoParCoarse;
    std::vector<IsoParameters* > IsoParFine;

    std::vector<Boundary *>  boundaryCoarse_;
    std::vector<Boundary *>  boundaryFine_;

    std::vector<GlueZone *>  glueZoneFine_;
	
private:
    
    // Meshes data
    int elemType;                           // 0 - FEM coarse mesh; 1 - IGA coarse mesh 
    int numElemCoarse;                      // Number of elements in the coarse mesh
    int numElemFine;                        // Number of elements in the fine mesh
    int numBoundElemCoarse;                 // Number of bondary elements in the coarse mesh
    int numBoundElemFine;                   // Number of bondary elements in the fine mesh
    int numElemGlueZoneFine;                // Number of fine elements in the gluing zone
    int numElemGlueZoneCoarse;              // Number of coarse elements in the gluing zone
    int numNodesCoarse;                     // Number of nodes in the coarse mesh
    int numNodesFine;		                // Number of nodes in the fine mesh
    int numNodesGlueZoneFine;               // Number of fine nodes in the gluing zone
    int numNodesGlueZoneCoarse;             // Number of coarse nodes in the gluing zone
    std::vector<int>elementsGlueZoneFine_;  // Vector with fine elements in the gluing zone
    std::vector<int>nodesGlueZoneFine_;     // Vector with fine nodes in the gluing zone
    std::vector<int>elementsGlueZoneCoarse_;// Vector with coarse elements in the gluing zone
    std::vector<int>nodesGlueZoneCoarse_;   // Vector with coarse nodes in the gluing zone
    
    // Integration parameters
    int numTimeSteps;           // Number of times steps
    double dTime;               // Time step
    int iTimeStep;              // Time counter   
    int integScheme;            // Integration scheme (0 - max. dissipation; 1 - without dissipation)
    
    //Arlequin parameters
    double glueZoneThickness;	// Thickness from gluing zone
    double arlequinEpsilon;		// Constant > 0
         
    //MPI parameters
    int rank;                   				// Name of the process
    std::pair<idx_t*,idx_t*> domDecompCoarse;   // Coarse Model Domain Decomposition
    std::pair<idx_t*,idx_t*> domDecompFine;     // Fine Model Domain Decomposition
    double pi = M_PI;
    
    int NumBezierNodesLagrange;   

public:
    // Sets the coarse and fine mesh models
    void setFluidModels(FluidMesh& coarse, FluidMesh& fine);
    
    // Computes and stores the element boxes for improving the correspondence searching process
    void setElementBoxes();
    
    // Computes and store the signaled distance function from Nodes or CP to a defined boundary
    void setSignaledDistance();

    // Defines the Gluing zone
    void setGluingZone();

    // Computes the energy weight function
    void setWeightFunction();

    //Sets the nodal and integration points correspondence in the coarse mesh
    void setCorrespondenceFine();

   	//Sets the Dirichelet Constrains in the domain
   	void setDirichletConstrain(std::vector<int> &dofTemp);
   
    // Searchs point correspondence in the coarse mesh
    void searchPointCorrespondence(double *x,std::vector<Nodes *> nodes,
                                  std::vector<Element *> elements,int numElem,
                                  double *xsiC, int &elemC, int elSearch);
    
    //Solves the Arlequin Problem
    int solveArlequinProblem(int iterNumber,double tolerance);

    //Prints the results for Paraview post-processing
    void printResults(int step);

    //Print the results of the gluing in the integration points for Paraview post-processing
    void printResultsIP(int step);

};

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------



//------------------------------------------------------------------------------
//-----------------------------SET ELEMENT BOXES--------------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::setElementBoxes() {
   
	int dim = 2;

    //Compute element boxes for coarse model (FEM mesh or IGA mesh)
    for (int jel = 0; jel < numElemCoarse; jel++){
        
        int *connec = elementsCoarse_[jel] -> getConnectivity();
        double xk[dim], Xk[dim];

        // if (elemType == 0) { //a box that envelop the triangle

        //     double *x1 = nodesCoarse_[connec[0]] -> getCoordinates();
        //     double *x2 = nodesCoarse_[connec[1]] -> getCoordinates();
        //     double *x3 = nodesCoarse_[connec[2]] -> getCoordinates();
                       
        //     for (int i = 0; i < dim; i++) {
        //     	xk[i] = std::min(x1[i],std::min(x2[i], x3[i]));
        //     	Xk[i] = std::max(x1[i],std::max(x2[i], x3[i]));
        //     } 
       
        //     elementsCoarse_[jel] -> setIntersectionParameters(xk, Xk);
        
        // } else { 

            //Bezier transformation matrix
        	double **MatrixC_;
        	MatrixC_ = new double*[9];
        	for (int i = 0; i < 9; i++) MatrixC_[i] = new double[9];
        	
        	elementsCoarse_[jel] -> getMatrixC(MatrixC_);
            
            //Transposing MatrixC
            double transMatrixC_[9][9];
            for (int i = 0; i<9; i++){
                for (int j = 0; j<9; j++){
                    transMatrixC_[i][j] = MatrixC_[j][i];
                };
            };
            
            //control points coordinates
            double coord_[9][2];
            for (int i = 0; i < 9; i++){
                double *xx = nodesCoarse_[connec[i]]->getCoordinates();
                for (int j = 0; j < dim; j++) coord_[i][j] = xx[j];
            };
            
            //Bezier control points coordinates
            double Bcoord_[9][2] = {};
            for (int i = 0; i<9; i++){
                for (int j = 0; j<9; j++){
                	for (int k = 0; k < dim; k++) Bcoord_[i][k] += transMatrixC_[i][j] * coord_[j][k];
                };
            };
            
            for (int i = 0; i < dim; i++){
            	xk[i] = Bcoord_[0][i];
            	Xk[i] = Bcoord_[8][i];
            }
                        
            elementsCoarse_[jel] -> setIntersectionParameters(xk, Xk);
        // }; //IGA elements

    };

    return;
};

//------------------------------------------------------------------------------
//--------------------------SEARCH CORRESPONDENCE-------------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::searchPointCorrespondence(double *x,std::vector<Nodes *> nodes, 
                                            std::vector<Element *> elements, 
                                            int numElem, double *xsiC, int &elemC, int elSearch){

    int dim = 2;
    QuadShapeFunction<2> shapeQuad;
    double x_[dim],deltaX[dim],deltaXsi[dim],xsi[dim],xsiCC[dim+1];

    double **ainv;
    ainv = new double*[dim];
    for (int i = 0; i < dim; ++i) ainv[i] = new double[dim];

    for (int i = 0; i<dim; i++){
       xsi[i] = 0.0; //central element cooordinates
       x_[i] = 0.0;
       xsiC[i] = 1.e50;
    };
    for (int i = 0; i<dim+1; i++) xsiCC[i] = 1.e10;
  
    int *connec = elements[elSearch] -> getConnectivity();
            
    //Computing basis functions
    double   phi_[9],wpc[9];
    for (int i = 0; i < 9; i++) wpc[i] = nodes[connec[i]] -> getWeightPC();
    int *INC_ = nodes[connec[8]] -> getINC();
    int patch = elements[elSearch] -> getPatch();
    shapeQuad.evaluateIso(xsi,phi_,wpc,INC_,IsoParCoarse,patch);

    for (int i = 0; i < 9; i++){
        double *xint = nodes[connec[i]] -> getCoordinates();
        for (int j = 0; j < dim; j++) x_[j] += xint[j] * phi_[i];                 
    };

    double error = 1.e6;
    int iterations = 0;

    while ((error > 1.e-8) && (iterations < 4)) {
        
        iterations++;

        for (int i = 0; i < dim; i++){
            deltaX[i] = x[i] - x_[i]; 
            deltaXsi[i] = 0.0;
        };

        elements[elSearch] -> getQuadJacobianMatrix_ISO(xsi,ainv);
        
        for (int i = 0; i <dim; i++){
            for (int j = 0; j<dim; j++){
                deltaXsi[i] += ainv[i][j] * deltaX[j];
            };
        };
        
        for (int i = 0; i < dim; i++){
           xsi[i] += deltaXsi[i];
           x_[i] = 0.0;
        };
                           
        shapeQuad.evaluateIso(xsi,phi_,wpc,INC_,IsoParCoarse,patch);
        
        for (int i = 0; i < 9; i++){
            double *xint = nodes[connec[i]] -> getCoordinates();
            for (int j = 0; j < dim; j++) x_[j] += xint[j]* phi_[i];                  
        };
                          
        error = sqrt(deltaXsi[0]*deltaXsi[0] + deltaXsi[1]*deltaXsi[1]);
    };
    
    double t1 = -1 - 1.e-2;
    double t2 =  1. + 1.e-2;   

    if ((xsi[0] >= t1) && (xsi[1] >= t1) &&
        (xsi[0] <= t2) && (xsi[1] <= t2)){

        xsiC[0] = xsi[0];
        xsiC[1] = xsi[1];
        elemC = elSearch;

    } else {

        for (int jel = 0; jel < numElem; jel++){

            int *connec = elements[jel] -> getConnectivity();

            // if (elemType == 0){ // FEM coarse mesh

            //     //get boxes information        
            //     std::pair<double*,double*> XK;
            //     XK = elements[jel] -> getXIntersectionParameter();

            //     //Chech if the node is inside the element box
            //     if ((x[0] < XK.first[0]) || (x[0] > XK.second[0]) ||
            //         (x[1] < XK.first[1]) || (x[1] > XK.second[1])) continue;
                
            //     //Compute nodal correspondence
            //     xsiCC[0] = 1.e10;
            //     xsiCC[1] = 1.e10;
            //     xsiCC[2] = 1.e10;
                
            //     xsi[0] = 1. / 3.;
            //     xsi[1] = 1. / 3.;

            //     double phi_[6];
            //     shapeQuad.evaluate(xsi,phi_);

            //     for (int i = 0; i <2 ; i++) x_[i] = 0.;
            //     for (int i = 0; i < 6; i++){
            //         double *xint = nodes[connec[i]] -> getCoordinates();
            //         double *xintp = nodes[connec[i]] -> getPreviousCoordinates();
            //         x_[0] += (alpha_f * xint[0] + (1. - alpha_f) * xintp[0]) * phi_[i];
            //         x_[1] += (alpha_f * xint[1] + (1. - alpha_f) * xintp[1]) * phi_[i];                    
            //     };

            //     double error = 1.e6;
                
            //     int iterations = 0;

            //     while ((error > 1.e-8) && (iterations < 4)) {
                    
            //         iterations++;
                    
            //         deltaX[0] = x[0] - x_[0];
            //         deltaX[1] = x[1] - x_[1];
                    
            //         elements[jel] -> getJacobianMatrixValues(xsi,ainv);

            //         double tempAinv[2][2];
            //         for (int i = 0; i <2; i++){
            //             for (int j =0; j<2; j++){
            //                 tempAinv[i][j] = ainv[j][i];
            //             }
            //         }

            //         for (int i = 0; i <2; i++){
            //             for (int j =0; j<2; j++){
            //                 ainv[i][j] = tempAinv[i][j];
            //             }
            //         }

            //         for (int i = 0; i <2; i++) deltaXsi[i] = 0.0;

            //         for (int i = 0; i <2; i++){
            //             for (int j = 0; j<2; j++){
            //                 deltaXsi[i] += ainv[i][j] * deltaX[j];
            //             }
            //         }
                    
            //         xsi[0] += deltaXsi[0];
            //         xsi[1] += deltaXsi[1];
                              
            //         shapeQuad.evaluate(xsi,phi_);
                    
            //         for (int i = 0; i <2 ; i++) x_[i] = 0.;
                    
            //         for (int i=0; i<6; i++){
            //             double *xint = nodes[connec[i]] -> getCoordinates();
            //             double *xintp = nodes[connec[i]] -> getPreviousCoordinates();
            //             x_[0] += (alpha_f * xint[0] + (1. - alpha_f) * xintp[0]) * phi_[i];
            //             x_[1] += (alpha_f * xint[1] + (1. - alpha_f) * xintp[1]) * phi_[i];                    
            //         };                   
                    

            //         error = sqrt(deltaXsi[0]*deltaXsi[0] + deltaXsi[1]*deltaXsi[1]);
                
            //     };
                
            //     double t1 = -1.e-2;
            //     double t2 =  1. - t1;
                
            //     xsiCC[0] = xsi[0];
            //     xsiCC[1] = xsi[1];       
            //     xsiCC[2] = 1. - xsiCC[0] - xsiCC[1];

            //     if ((xsiCC[0] >= t1) && (xsiCC[1] >= t1) && (xsiCC[2] >= t1) &&
            //         (xsiCC[0] <= t2) && (xsiCC[1] <= t2) && (xsiCC[2] <= t2)){

            //         xsiC[0] = xsi[0];
            //         xsiC[1] = xsi[1];
            //         elemC = jel;
            //     };           
       
            // } else { //IGA coarse mesh
                          
                //get boxes information  
                std::pair<double*,double*> XK;      
                XK = elements[jel] -> getXIntersectionParameter();

                //Chech if the node is inside the element box
                if ((x[0] < XK.first[0]) || (x[0] > XK.second[0]) ||
                    (x[1] < XK.first[1]) || (x[1] > XK.second[1])) continue;
                     
                //Compute the basis functions
                double phi_[9],wpc[9];
                for (int i = 0; i < 9; i++) wpc[i] = nodes[connec[i]] -> getWeightPC();
                int *INC_ = nodes[connec[8]] -> getINC();
                int patch = elements[jel] -> getPatch();
                
               
                for (int i = 0; i < dim; i++) {
                    xsi[i] = 0.0; //central element cooordinates
                    x_[i] = 0.0;
                };

                shapeQuad.evaluateIso(xsi,phi_,wpc,INC_,IsoParCoarse,patch);
               
                for (int i = 0; i < 9; i++){
                    double *xint = nodes[connec[i]] -> getCoordinates();
                    for (int j = 0; j < dim; j++) x_[j] += xint[j] * phi_[i];                  
                };

                double error = 1.e6;
                
                int iterations = 0;

                while ((error > 1.e-8) && (iterations < 4)) {
                    
                    iterations++;
                    
                    for (int i = 0; i < dim; i++){
                        deltaX[i] = x[i] - x_[i]; 
                        deltaXsi[i] = 0.0;
                    };

                    elements[jel] -> getQuadJacobianMatrix_ISO(xsi,ainv);
                    
                    for (int i = 0; i <dim; i++){
                        for (int j = 0; j<dim; j++){
                            deltaXsi[i] += ainv[i][j] * deltaX[j];
                        };
                    };
                    
                    for (int i = 0; i < dim; i++){
                        xsi[i] += deltaXsi[i];
                        x_[i] = 0.0;
                    };
                               
                    shapeQuad.evaluateIso(xsi,phi_,wpc,INC_,IsoParCoarse,patch);
                    
                    for (int i = 0; i < 9; i++){
                        double *xint = nodes[connec[i]] -> getCoordinates();
                        for (int j = 0; j < dim; j++) x_[j] += xint[j] * phi_[i];                    
                    };
                                      
                    error = sqrt(deltaXsi[0]*deltaXsi[0] + deltaXsi[1]*deltaXsi[1]);
                };
                
                double t1 = -1 - 1.e-2;
                double t2 =  1. + 1.e-2;
                    
                if ((xsi[0] >= t1) && (xsi[1] >= t1) &&
                    (xsi[0] <= t2) && (xsi[1] <= t2)){
                    xsiC[0] = xsi[0];
                    xsiC[1] = xsi[1];
                    elemC = jel;
                };           
            // };   
        }; 

    };     
    
    if (fabs(xsi[0]) > 2.) std::cout << "PROBEM SEARCHING NODE CORRESPONDENCE " 
                                         << std::endl;  


    for (int i = 0; i < 2; ++i) delete [] ainv[i];
    delete [] ainv;

};

template<>
void Arlequin<2>::setCorrespondenceFine() {

    int dim = 2;
    double &alpha_f = parametersFine -> getAlphaF();

    //Node correspondence
    for (int inode = 0; inode < numNodesGlueZoneFine; inode++) {


    	//std::cout << nodesGlueZoneFine_[inode] << std::endl;
        
        double* x = nodesFine_[nodesGlueZoneFine_[inode]] -> getCoordinates();

        int elemC = 0;
        double xsiC[dim] = {};

        searchPointCorrespondence(x, nodesCoarse_, elementsCoarse_,elementsCoarse_.size(),xsiC,elemC,
                                  nodesFine_[nodesGlueZoneFine_[inode]] -> getNodalElemCorrespondence());
      
        nodesFine_[nodesGlueZoneFine_[inode]] -> setNodalCorrespondence(elemC,xsiC);             

    };
    
    //integration points correspondence
    for (int i = 0; i< numElemGlueZoneFine; i++){

        SpecialQuadrature squad;
        
        int *connec = elementsFine_[elementsGlueZoneFine_[i]] -> getConnectivity();
        int patch = elementsFine_[elementsGlueZoneFine_[i]] -> getPatch();
        int *inc = nodesFine_[connec[8]] -> getINC();

        double x1[9], x2[9], wpc[9];
        for (int j = 0; j < 9; j++){
            double *x = nodesFine_[connec[j]] -> getCoordinates();
            double *xp = nodesFine_[connec[j]] -> getPreviousCoordinates();
            wpc[j] = nodesFine_[connec[j]] -> getWeightPC();
            x1[j] = alpha_f * x[0] + (1. - alpha_f) * xp[0];
            x2[j] = alpha_f * x[1] + (1. - alpha_f) * xp[1];
        };   

            
        int numberIntPoints = elementsFine_[elementsGlueZoneFine_[i]] -> 
                              getNumberOfIntegrationPointsSpecial();

        for (int ip = 0; ip < numberIntPoints; ip++){

            int elemC = 0;
            double xsiC[dim] = {};
            double x_[dim] = {};

            x_[0] = squad.interpolateQuadraticVariable(x1,ip,wpc,inc,IsoParFine,patch);
            x_[1] = squad.interpolateQuadraticVariable(x2,ip,wpc,inc,IsoParFine,patch);

            searchPointCorrespondence(x_,nodesCoarse_,elementsCoarse_,
                                      elementsCoarse_.size(),xsiC,elemC, 
                                      elementsFine_[elementsGlueZoneFine_[i]] -> getIntegPointCorrespondenceElement(ip));

            
            elementsFine_[elementsGlueZoneFine_[i]] -> setIntegrationPointCorrespondence(ip,xsiC,elemC);
        }            
    }
}

//------------------------------------------------------------------------------
//-------------------------FIND SIGNALADED DISTANCE-----------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::setSignaledDistance(){

	int dim = 2;

    double  x[dim], x1[dim], x1B[dim], x2[dim], x2B[dim], x3B[dim];
    double  n[dim], test[dim];
    int bconnec[3];
    double dist;

    for (int inode = 0 ; inode < numNodesFine; inode++){ 
        nodesFine_[inode] -> clearInnerNormal();
        nodesFine_[inode] -> setDistFunction(0.0);
    };
   
    for (int inode = 0 ; inode < numNodesCoarse; inode++){ 
        nodesCoarse_[inode] -> setDistFunction(0.0);
    };
  
    // Normal vector from a defined boundary in fine mesh
     for (int ibound = 0; ibound < numBoundElemFine; ibound++){

        if (boundaryFine_[ibound]-> getConstrain(0) == 2){

        	int *connec = elementsFine_[boundaryFine_[ibound] -> getElement()] -> getConnectivity();

        	//Recognizing the element side
            std::pair<std::vector<int>, std::vector<int> > elemBound;
            elemBound = elementsFine_[boundaryFine_[ibound] -> getElement()] -> getElemSideInBoundary();
            int numofBoundaries = elemBound.first.size(); 
            int side;
            for (int j = 0; j < numofBoundaries; j++) {
                if (elemBound.second[j] == ibound ) {
                    side = elemBound.first[j];
                };
            };

            //Finding the coorrespondent coordinates in bezier element
            double **MatrixC_;
            MatrixC_ = new double*[9];
            for (int i = 0; i < 9; i++) MatrixC_[i] = new double[9];
            elementsFine_[boundaryFine_[ibound] -> getElement()] -> getMatrixC(MatrixC_);

            //Transposing MatrixC
            double transMatrixC_[9][9];
            for (int i = 0; i<9; i++){
                for (int j = 0; j<9; j++){
                    transMatrixC_[i][j] = MatrixC_[j][i];
                };
            };
        	
        	double coord_[9][2];
        	for (int i = 0; i < 9; i++){
                double *xx = nodesFine_[connec[i]]->getCoordinates();
                for (int j = 0; j < dim; j++) coord_[i][j] = xx[j];
            };
            
            //Bezier Coordinates
            double Bcoord_[9][2] = {};
            for (int i = 0; i<9; i++){
                for (int j = 0; j<9; j++){
                    for (int k = 0; k < dim; k++) Bcoord_[i][k] += transMatrixC_[i][j] * coord_[j][k];
                };
            };
   
            if (side == 0){                 
                bconnec[0] = connec[0];
                bconnec[1] = connec[2];
                bconnec[2] = connec[1];
                for (int i = 0; i < dim; i++){
                	x1B[i]= Bcoord_[0][i];
                	x2B[i]= Bcoord_[2][i];
                	x3B[i]= Bcoord_[1][i];
                };
            };
            if (side == 1){
                bconnec[0] = connec[2];
                bconnec[1] = connec[8];
                bconnec[2] = connec[5];
                for (int i = 0; i < dim; i++){
	                x1B[i]= Bcoord_[2][i];
	                x2B[i]= Bcoord_[8][i];
	                x3B[i]= Bcoord_[5][i];
                };
			};
            if (side == 2){
                bconnec[0] = connec[8];
                bconnec[1] = connec[6];
                bconnec[2] = connec[7];
                for (int i = 0; i < dim; i++){
                	x1B[i]= Bcoord_[8][i];
	                x2B[i]= Bcoord_[6][i];
	                x3B[i]= Bcoord_[7][i];
                };	                
            };
            if (side == 3){
                bconnec[0] = connec[6];
                bconnec[1] = connec[0];
                bconnec[2] = connec[3];
                for (int i = 0; i < dim; i++){
		            x1B[i]= Bcoord_[6][i];
	                x2B[i]= Bcoord_[0][i];
	                x3B[i]= Bcoord_[3][i];
	            }; 
            };

            //first segment
            int no1 = bconnec[0];
            int no2 = bconnec[2];
            for (int i = 0; i < dim; i++){
				x1[i] = x1B[i];
	            x2[i] = x3B[i];
            }
            
            double sLength = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                                  (x1[0] - x2[0]) * (x1[0] - x2[0]));
            
            n[0] = (x2[1] - x1[1]) / sLength;
            n[1] = (x1[0] - x2[0]) / sLength;

            nodesFine_[no1] -> setInnerNormal(n);
            nodesFine_[no2] -> setInnerNormal(n);

            //second segment
            no1 = bconnec[2];
            no2 = bconnec[1];
            for (int i = 0; i <dim; i++){
            	x1[i] = x3B[i];
	            x2[i] = x2B[i];
            };
            
            sLength = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                           (x1[0] - x2[0]) * (x1[0] - x2[0]));

            n[0] = (x2[1] - x1[1]) / sLength;
            n[1] = (x1[0] - x2[0]) / sLength;

            nodesFine_[no1] -> setInnerNormal(n);
            nodesFine_[no2] -> setInnerNormal(n);   

        };//constrain == 2
    }; //numBoundFine



    // lembrar que a normal do ultimo elemento com aquele ponto de controle (bezier control point) 
    // vai me dar a normal final do ponto de controle

    //Gambiarra para os cantos para o problema da cavidade com todas as bordas com malhas refinidas
    // n[0] = 1.;
    // n[1] = 1.;

    // nodesFine_[334] -> setInnerNormal(n);
    // nodesFine_[469] -> setInnerNormal(n);

    // n[0] = -1.;
    // n[1] = 1.;

    // nodesFine_[359] -> setInnerNormal(n);
    // nodesFine_[720] -> setInnerNormal(n);


    // n[0] = -1.;
    // n[1] = -1.;

    // nodesFine_[1105] -> setInnerNormal(n);
    // nodesFine_[970] -> setInnerNormal(n);

    
    // n[0] = 1.;
    // n[1] = -1.;

    // nodesFine_[1080] -> setInnerNormal(n);
    // nodesFine_[719] -> setInnerNormal(n);


    //Coarse mesh - closer distance for nodes or control points from defined fine boundary 
    for (int ino = 0; ino < numNodesCoarse; ino++){
        
        double *xx = nodesCoarse_[ino]->getCoordinates();
        
        for (int i = 0; i < dim; i++) x[i] = xx[i];
        
        dist=10000000000000000000000000000.;
        
        for (int ibound = 0; ibound < numBoundElemFine; ibound++){
            
            if (boundaryFine_[ibound] -> getConstrain(0) == 2){

                int *connec = elementsFine_[boundaryFine_[ibound] -> getElement()] -> getConnectivity();
                
                //side in the boundary
                std::pair<std::vector<int>, std::vector<int> > elemBound;
                elemBound = elementsFine_[boundaryFine_[ibound] -> getElement()] -> getElemSideInBoundary();
                int numofBoundaries = elemBound.first.size();
                int side;
                for (int j = 0; j < numofBoundaries; j++) {
                    if (elemBound.second[j] == ibound ) {
                        side = elemBound.first[j];
                    };
                };
                
                //gets the Bezier element coordinates        
                double **MatrixC_;
                MatrixC_ = new double*[9];
                for (int i = 0; i < 9; i++) MatrixC_[i] = new double[9];
                elementsFine_[boundaryFine_[ibound] -> getElement()] -> getMatrixC(MatrixC_);

                //Transposing MatrixC
                double transMatrixC_[9][9];
                for (int i = 0; i<9; i++){
                    for (int j = 0; j<9; j++){
                        transMatrixC_[i][j] = MatrixC_[j][i];
                    };
                };
            	
            	double coord_[9][2];
            	for (int i = 0; i < 9; i++){
                    double *xx = nodesFine_[connec[i]]->getCoordinates();
                    for (int j = 0; j < dim; j++) coord_[i][j] = xx[j];
                };

	            //Bezier Coordinates
	            double Bcoord_[9][2] = {};
	            for (int i = 0; i<9; i++){
	                for (int j = 0; j<9; j++){
	                	for (int k = 0; k < dim; k++) Bcoord_[i][k] += transMatrixC_[i][j] * coord_[j][k];
	                };
	            };
	       
	            if (side == 0){                 
	                bconnec[0] = connec[0];
	                bconnec[1] = connec[2];
	                bconnec[2] = connec[1];
	                for (int i = 0; i < dim; i++){
	                	x1B[i]= Bcoord_[0][i];
	                	x2B[i]= Bcoord_[2][i];
	                	x3B[i]= Bcoord_[1][i];
	                };
	            };
	            if (side == 1){
	                bconnec[0] = connec[2];
	                bconnec[1] = connec[8];
	                bconnec[2] = connec[5];
	                for (int i = 0; i < dim; i++){
		                x1B[i]= Bcoord_[2][i];
		                x2B[i]= Bcoord_[8][i];
		                x3B[i]= Bcoord_[5][i];
	                };
				};
	            if (side == 2){
	                bconnec[0] = connec[8];
	                bconnec[1] = connec[6];
	                bconnec[2] = connec[7];
	                for (int i = 0; i < dim; i++){
	                	x1B[i]= Bcoord_[8][i];
		                x2B[i]= Bcoord_[6][i];
		                x3B[i]= Bcoord_[7][i];
	                };	                
	            };
	            if (side == 3){
	                bconnec[0] = connec[6];
	                bconnec[1] = connec[0];
	                bconnec[2] = connec[3];
	                for (int i = 0; i < dim; i++){
			            x1B[i]= Bcoord_[6][i];
		                x2B[i]= Bcoord_[0][i];
		                x3B[i]= Bcoord_[3][i];
		            }; 
	            };

                //first segment
                int no1,no2;
                no1 = bconnec[0];
                no2 = bconnec[2];
                for (int i = 0; i <dim; i++){
                	x1[i]= x1B[i];
                	x2[i]= x3B[i];
                }
                
                double aux0 =  sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                                    (x2[0] - x1[0]) * (x2[0] - x1[0]));
                double aux1 = ((x[0] - x1[0]) * (x2[0] - x1[0])+
                               (x[1] - x1[1]) * (x2[1] - x1[1])) / aux0;
                double dist2 =-((x2[1] - x1[1]) * x[0] - 
                                (x2[0] - x1[0]) * x[1] +
                                x2[0] * x1[1] - x2[1] * x1[0]) / aux0;
                
                if (aux1 > aux0){
                    
                    dist2 = sqrt((x2[1] - x[1]) * (x2[1] - x[1]) +
                                 (x2[0] - x[0]) * (x2[0] - x[0]));
                    //find signal
                    //side normal vector
                    double *normal = nodesFine_[no2] -> getInnerNormal();
                    
                    test[0] = x[0] - x2[0];
                    test[1] = x[1] - x2[1];

                    double signaltest = 0.0;
                    for (int i=0; i<2; i++){
                    	signaltest += normal[i] * test[i];
                    }
                   
                    double signal = -1.;
                    if (signaltest <= -0.001)signal = 1.;
                    
                    dist2 *= signal;
                };

                if (aux1 < 0.){
                    
                    dist2 = sqrt((x[1] - x1[1]) * (x[1] - x1[1]) +
                                 (x[0] - x1[0]) * (x[0] - x1[0]));
                    //find signal
                    //side normal vector
                    double *normal = nodesFine_[no1] -> getInnerNormal();
                    
                    test[0] = x[0] - x1[0];
                    test[1] = x[1] - x1[1];
                    
                    double signaltest = 0.0;
                    for (int i=0; i<2; i++){
                    	signaltest += normal[i] * test[i];
                    }
                    
                    double signal = -1.;
                    if (signaltest <= -0.001) signal = 1.;
                    
                    dist2 *= signal;
                };
                
                if (fabs(dist2) < fabs(dist)) dist = dist2;
                
                //second segment
                no1 = bconnec[2];
                no2 = bconnec[1];
                for (int i = 0; i < dim; i++){
                	x1[i]= x3B[i];
                	x2[i]= x2B[i];
                }

                aux0 = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                            (x2[0] - x1[0]) * (x2[0] - x1[0]));
                aux1 = ((x[0] - x1[0]) * (x2[0] - x1[0]) + 
                        (x[1] - x1[1]) * (x2[1] - x1[1])) / aux0;
                dist2 = -((x2[1] - x1[1]) * x[0] - (x2[0] - x1[0]) * x[1] +
                          x2[0] * x1[1] - x2[1] * x1[0]) / aux0;

                if (aux1 > aux0){
                    dist2 = sqrt((x2[1] - x[1]) * (x2[1] - x[1]) +
                                 (x2[0] - x[0]) * (x2[0] - x[0]));
                    
                    double *normal = nodesFine_[no2] -> getInnerNormal();
                    
                    test[0] = x[0] - x2[0];
                    test[1] = x[1] - x2[1];
                    
                    double signaltest = 0.0;
                    for (int i=0; i<2; i++){
                    	signaltest += normal[i] * test[i];
                    }

                    double signal = -1.;
                    if (signaltest <= -0.001)signal = 1.;
                    
                    dist2 *= signal;
                                       
                };

                if (aux1 < 0.){
                    dist2 = sqrt((x[1] - x1[1]) * (x[1] - x1[1]) +  
                                 (x[0] - x1[0]) * (x[0] - x1[0]));
                    
                    //find signal
                    //side normal vector
                    double *normal;
                    normal = nodesFine_[no1] -> getInnerNormal();
                    
                    test[0] = x[0] - x1[0];
                    test[1] = x[1] - x1[1];
                    
                    double signaltest = 0.0;
                    for (int i=0; i<2; i++){
                    	signaltest += normal[i] * test[i];
                    }
                    
                    double signal = -1.;
                    if (signaltest <= -0.001)signal = 1.;
                    
                    dist2 *= signal;
                    
                };
                if (fabs(dist2) < fabs(dist)) dist = dist2;
            }; //if bf is the blend boundary
        }; // numboundaryfine
    
        if (fabs(nodesCoarse_[ino] -> getDistFunction()) < 1.e-2){
            nodesCoarse_[ino] -> setDistFunction(dist); 
        };
    
    }; //numnodescoarse


    //Fine mesh - closer distance for nodes or control points from defined fine boundary 
    for (int ino = 0; ino < numNodesFine; ino++){
        
        double &alpha_f = parametersFine -> getAlphaF();  

        double *xx=nodesFine_[ino]->getCoordinates();

        for (int i = 0; i < dim; i++) x[i] = xx[i];

        
        dist=10000000000000000000000000000.;
               
        for (int ibound = 0; ibound < numBoundElemFine; ibound++){
            
            if (boundaryFine_[ibound] -> getConstrain(0) == 2){
                
                int *connec = elementsFine_[boundaryFine_[ibound] -> getElement()] -> getConnectivity();

                std::pair<std::vector<int>, std::vector<int> > elemBound;
                elemBound = elementsFine_[boundaryFine_[ibound] -> getElement()] -> getElemSideInBoundary();
                int numofBoundaries = elemBound.first.size();
                int side;
                for (int j = 0; j < numofBoundaries; j++) {
                    if (elemBound.second[j] == ibound ) {
                        side = elemBound.first[j];
                    };
                }  ;
                    
                //Finding the coorrespondent coordinates in bezier element
                double **MatrixC_;
                MatrixC_ = new double*[9];
                for (int i = 0; i < 9; i++) MatrixC_[i] = new double[9];
                elementsFine_[boundaryFine_[ibound] -> getElement()] -> getMatrixC(MatrixC_);

            	//Transposing MatrixC
                double transMatrixC_[9][9];
                for (int i = 0; i<9; i++){
                    for (int j = 0; j<9; j++){
                        transMatrixC_[i][j] = MatrixC_[j][i];
                    };
                };
                
                double coord_[9][2];
                for (int i = 0; i < 9; i++){
                    double *xx = nodesFine_[connec[i]]->getCoordinates();
                    for (int j = 0; j < dim; j++) coord_[i][j] = xx[j];
                };

                //Bezier Coordinates
                double Bcoord_[9][2] = {};
                for (int i = 0; i<9; i++){
                    for (int j = 0; j<9; j++){
                    	for (int k = 0; k < dim; k++) Bcoord_[i][k] += transMatrixC_[i][j] * coord_[j][k];
                    };
                };
                
                if (side == 0){                 
	                bconnec[0] = connec[0];
	                bconnec[1] = connec[2];
	                bconnec[2] = connec[1];
	                for (int i = 0; i < dim; i++){
	                	x1B[i]= Bcoord_[0][i];
	                	x2B[i]= Bcoord_[2][i];
	                	x3B[i]= Bcoord_[1][i];
	                };
	            };
	            if (side == 1){
	                bconnec[0] = connec[2];
	                bconnec[1] = connec[8];
	                bconnec[2] = connec[5];
	                for (int i = 0; i < dim; i++){
		                x1B[i]= Bcoord_[2][i];
		                x2B[i]= Bcoord_[8][i];
		                x3B[i]= Bcoord_[5][i];
	                };
				};
	            if (side == 2){
	                bconnec[0] = connec[8];
	                bconnec[1] = connec[6];
	                bconnec[2] = connec[7];
	                for (int i = 0; i < dim; i++){
	                	x1B[i]= Bcoord_[8][i];
		                x2B[i]= Bcoord_[6][i];
		                x3B[i]= Bcoord_[7][i];
	                };	                
	            };
	            if (side == 3){
	                bconnec[0] = connec[6];
	                bconnec[1] = connec[0];
	                bconnec[2] = connec[3];
	                for (int i = 0; i < dim; i++){
			            x1B[i]= Bcoord_[6][i];
		                x2B[i]= Bcoord_[0][i];
		                x3B[i]= Bcoord_[3][i];
		            }; 
	            };

                //first segment
                int no1 = bconnec[0];
                int no2 = bconnec[2];
                for (int i = 0; i < dim; i++){
                	x1[i]= x1B[i];
                	x2[i]= x3B[i];
                }

                double aux0 =  sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                                    (x2[0] - x1[0]) * (x2[0] - x1[0]));
                double aux1 = ((x[0] - x1[0]) * (x2[0] - x1[0])+
                               (x[1] - x1[1]) * (x2[1] - x1[1])) / aux0;
                double dist2 =-((x2[1] - x1[1]) * x[0] - 
                                (x2[0] - x1[0]) * x[1] +
                                x2[0] * x1[1] - x2[1] * x1[0]) / aux0;
                
                
                if (aux1 > aux0){
                    dist2 = sqrt((x2[1] - x[1]) * (x2[1] - x[1]) +
                                 (x2[0] - x[0]) * (x2[0] - x[0]));
                    //find signal
                    //side normal vector
                    
                    double *normal = nodesFine_[no2] -> getInnerNormal();
                    
                    test[0] = x[0] - x2[0];
                    test[1] = x[1] - x2[1];

                    double signaltest = 0.0;
                    for (int i = 0; i< 2; i++){
                        signaltest += test[i] * normal[i];
                    }
                    
                    double signal = -1.;
                    if (signaltest <= -0.001)signal = 1.;
                    
                    dist2 *= signal;
                };

                if (aux1 < 0.){
                    dist2 = sqrt((x[1] - x1[1]) * (x[1] - x1[1]) +
                                 (x[0] - x1[0]) * (x[0] - x1[0]));
                    //find signal
                    //side normal vector
                    double *normal =  nodesFine_[no1] -> getInnerNormal();
                    
                    test[0] = x[0] - x1[0];
                    test[1] = x[1] - x1[1];

                    double signaltest = 0.0;
                    for (int i = 0; i< 2; i++){
                        signaltest += test[i] * normal[i];
                    }
                    
                    double signal = -1.;
                    if (signaltest <= -0.001) signal = 1.;
                    
                    dist2 *= signal;
                };
                
                if (fabs(dist2) < fabs(dist)) {
                    dist = dist2;
                };
                
                //second segment
                no1 = bconnec[2];
                no2 = bconnec[1];
                for (int i = 0; i < dim; i++){
					x1[i]= x3B[i];
	                x2[i]= x2B[i];
                }

                aux0 = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
                            (x2[0] - x1[0]) * (x2[0] - x1[0]));
                aux1 = ((x[0] - x1[0]) * (x2[0] - x1[0]) + 
                        (x[1] - x1[1]) * (x2[1] - x1[1])) / aux0;
                dist2 = -((x2[1] - x1[1]) * x[0] - (x2[0] - x1[0]) * x[1] +
                          x2[0] * x1[1] - x2[1] * x1[0]) / aux0;

                if (aux1 > aux0){
                    dist2 = sqrt((x2[1] - x[1]) * (x2[1] - x[1]) +
                                 (x2[0] - x[0]) * (x2[0] - x[0]));
                    
                    double *normal = nodesFine_[no2] -> getInnerNormal();
                    
                    test[0] = x[0] - x2[0];
                    test[1] = x[1] - x2[1];

                    double signaltest = 0.0;
                    for (int i = 0; i< 2; i++){
                        signaltest += test[i] * normal[i];
                    };

                    double signal = -1.;
                    if (signaltest <= -0.001)signal = 1.;
                    
                    dist2 *= signal;
                                       
                };

                if (aux1 < 0.){
                    dist2 = sqrt((x[1] - x1[1]) * (x[1] - x1[1]) +
                                 (x[0] - x1[0]) * (x[0] - x1[0]));
                    //find signal
                    //side normal vector
                    double *normal = nodesFine_[no1] -> getInnerNormal();
                    
                    test[0] = x[0] - x1[0];
                    test[1] = x[1] - x1[1];
                    
                    double signaltest = 0.0;
                    for (int i = 0; i< 2; i++){
                        signaltest += test[i] * normal[i];
                    };
                    
                    double signal = -1.;
                    if (signaltest <= -0.001)signal = 1.;
                    
                    dist2 *= signal;
                    
                };

                if (fabs(dist2) < fabs(dist)) {
                    dist = dist2;
                };

            }; //if bf is the blend boundary  
            
        }; // numboundfine

        if(dist < 0) dist = 0;
        nodesFine_[ino] -> setDistFunction(dist);

    };//numNodes

}; //função


//------------------------------------------------------------------------------
//------------------------------SETS GLUING ZONE--------------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::setGluingZone(){

	int dim = 2;
    int flag;
    int nodesCZ[numNodesFine];
    int nodesCZ2[numNodesCoarse];
    QuadShapeFunction<2> shapeQuad;

    for (int i = 0; i < numNodesFine; i++) nodesCZ[i] = 0;    
    elementsGlueZoneFine_.reserve(numElemFine / 3);
    glueZoneFine_.reserve(numElemFine/3);
    nodesGlueZoneFine_.reserve(numNodesFine / 3);
    
    for (int i = 0; i < numNodesCoarse; i++) nodesCZ2[i] = 0;    
    elementsGlueZoneCoarse_.reserve(numElemCoarse / 3);
    nodesGlueZoneCoarse_.reserve(numNodesCoarse / 3);

    //parametric coordinates from Bézier control points
    double xsiCP[9][2];
    xsiCP[0][0] = -1.; xsiCP[0][1] = -1.;
    xsiCP[1][0] = 0.; xsiCP[1][1] = -1.;
    xsiCP[2][0] = 1.; xsiCP[2][1] = -1.;
    xsiCP[3][0] = -1.; xsiCP[3][1] = 0.;
    xsiCP[4][0] = 0.; xsiCP[4][1] = 0.;
    xsiCP[5][0] = 1.; xsiCP[5][1] = 0.;
    xsiCP[6][0] = -1.; xsiCP[6][1] = 1.;
    xsiCP[7][0] = 0.; xsiCP[7][1] = 1.;
    xsiCP[8][0] = 1.; xsiCP[8][1] = 1.;   

    // Defines a criterion to select the fine elements that are in the gluing zone
    int index = 0;
    for (int jel = 0; jel < numElemFine; jel++){
        
        int *connec = elementsFine_[jel] -> getConnectivity();
       	               
        double distance_[9];
        for (int i = 0; i < 9; i++){
            distance_[i] = nodesFine_[connec[i]] -> getDistFunction(); 
        };

        double Bdistance_[9] = {};
        for (int icp = 0; icp< 9; icp++){
            
            double wpc[9],phi_[9],xsi[dim];
            for (int i = 0; i <9; i++) wpc[i] = nodesFine_[connec[i]] -> getWeightPC();
            int *inc_ = nodesFine_[connec[8]] -> getINC();
            int patch_ = elementsFine_[jel] -> getPatch();
            for (int i = 0; i < dim; i++) xsi[i] = xsiCP[icp][i];
            shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,IsoParFine,patch_);
            
            for (int i = 0; i < 9; i++){
                Bdistance_[icp] += phi_[i] * distance_[i];
            };
        };

        flag = 0;
        for (int ino = 0; ino < 9; ino++){
            if (Bdistance_[ino] <= glueZoneThickness + 0.00001){
                flag += 1;
            };
        };

        if (flag == 9) {
            elementsGlueZoneFine_.push_back(jel);
            elementsFine_[jel] -> setGlueZone();

            GlueZone *el = new GlueZone(index++,jel);
            glueZoneFine_.push_back(el);

    	};
    };

    //Defines which nodes are in the gluing zone
    numElemGlueZoneFine = elementsGlueZoneFine_.size();

    for (int i = 0; i < numElemGlueZoneFine; i++){

        int *connec = elementsFine_[elementsGlueZoneFine_[i]] -> getConnectivity();

        for (int ino = 0; ino < 9; ino++){
            nodesCZ[connec[ino]] += 1;
            nodesFine_[connec[ino]] -> setGlueZone();
        };
    };

    //Compute and save number of nodes in the gluing zone
    numNodesGlueZoneFine = 0;
    for (int i = 0; i < numNodesFine; i++){
        if(nodesCZ[i] > 0) {
            numNodesGlueZoneFine += 1;
            nodesGlueZoneFine_.push_back(i);
        };
    };

    //Defining the Lagrange multipliers in the fine mesh
    for (int i = 0; i < numNodesGlueZoneFine; i++){
        
        double* x = nodesFine_[nodesGlueZoneFine_[i]] -> getCoordinates();
        
        Nodes *no = new Nodes(x,i,1.);
        nodesLagrangeFine_.push_back(no);
    };

    // Define the Lagrange Multipliers connectivity in the fine mesh
    for (int i = 0; i < numElemGlueZoneFine; i++){
        
        int *connecAux;
        connecAux = new int[9];
        
        int *connec = elementsFine_[elementsGlueZoneFine_[i]] -> getConnectivity();
       
        for (int ino = 0; ino < numNodesGlueZoneFine; ino++)
            for (int k = 0; k < 9; k++)
                if (nodesGlueZoneFine_[ino] == connec[k]) connecAux[k] = ino;
        
        glueZoneFine_[i] -> setConnectivity(connecAux);

    };

    
    // Defines a criterion to select the coarse elements that are in the gluing zone
    for (int jel = 0; jel < numElemCoarse; jel++){
        
        int *connec = elementsCoarse_[jel] -> getConnectivity();
        
        flag = 0;
        // if (elemType == 0){          
            
        //     for (int ino = 0; ino < 6; ino++){
        //         double dist = nodesCoarse_[connec[ino]] -> getDistFunction();
        //         if ((dist <= glueZoneThickness + 0.00001) && (dist >= 0.00001)){
        //             flag += 1;
        //         };
        //     };

        //     if (flag != 0) {
        //         elementsGlueZoneCoarse_.push_back(jel);
        //         elementsCoarse_[jel] -> setGlueZone();
        //     };  
        
        // } else {
            
            double distance_[9];
            for (int i = 0; i < 9; i++){
                distance_[i] = nodesCoarse_[connec[i]]->getDistFunction(); 
            };

            double Bdistance_[9] = {};
            for (int icp = 0; icp< 9; icp++){
                
                double wpc[9],phi_[9],xsi[dim];
                for (int i = 0; i <9; i++) wpc[i] = nodesCoarse_[connec[i]] -> getWeightPC();
                int *inc_ = nodesCoarse_[connec[8]] -> getINC();
                int patch_ = elementsCoarse_[jel] -> getPatch();
                for (int i = 0; i < dim; i++) xsi[i] = xsiCP[icp][i];
                shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,IsoParCoarse,patch_);

                for (int i = 0; i < 9; i++){
                    Bdistance_[icp] += phi_[i] * distance_[i];
                };
            };

            for (int ino = 0; ino < 9; ino++){
                if ((Bdistance_[ino] <= glueZoneThickness + 0.00001) && (Bdistance_[ino]  >= 0.00001)){
                     flag += 1;
                };
            };

            if (flag != 0) {
                elementsGlueZoneCoarse_.push_back(jel);
                elementsCoarse_[jel] -> setGlueZone();
        	};
        // }; //else IGA element
      
    };


    //Defines which coarse nodes are in the gluing zone
    numElemGlueZoneCoarse = elementsGlueZoneCoarse_.size();
    for (int i = 0; i < numElemGlueZoneCoarse; i++){
        int *connec = elementsCoarse_[elementsGlueZoneCoarse_[i]] -> getConnectivity();
        // if (elemType == 0){
        //     for (int ino = 0; ino < 6; ino++){
        //         nodesCZ2[connec[ino]] += 1;
        //         nodesCoarse_[connec[ino]] -> setGlueZone();
        //     };
        // } else{
            for (int ino = 0; ino < 9; ino++){
                nodesCZ2[connec[ino]] += 1;
                nodesCoarse_[connec[ino]] -> setGlueZone();
            };
        // };
    };

    //Compute number of coarse nodes in the gluing zone
    numNodesGlueZoneCoarse = 0;
    for (int i = 0; i < numNodesCoarse; i++){
        if(nodesCZ2[i] > 0) {
            numNodesGlueZoneCoarse += 1;
            nodesGlueZoneCoarse_.push_back(i);
        };
    };
   
};

template<>
void Arlequin<2>::setWeightFunction(){

	double wFuncValue;

	glueZoneThickness *= 1.01;	// Thickness from gluing zone

    //IGA COARSE MESH or FEM mesh
    for (int iNode = 0; iNode< numNodesCoarse; iNode++){

    	double r = nodesCoarse_[iNode] ->getDistFunction();
        
        if (r < 0){
            wFuncValue = 1.;
        } else {
            if (r >= glueZoneThickness){
                wFuncValue = arlequinEpsilon;
            } else {
                wFuncValue = 1. - (1. - arlequinEpsilon) / glueZoneThickness * r;                
                if (wFuncValue < arlequinEpsilon) wFuncValue = arlequinEpsilon;
            };
        };  
              
        nodesCoarse_[iNode] -> setWeightFunction(wFuncValue);

    };  //ielem    


    for (int jel = 0; jel < numElemCoarse; jel++){
        elementsCoarse_[jel] -> setIntegPointWeightFunction_ISO();        
    };


    //IGA FINE MESH
    for (int iNode = 0; iNode< numNodesFine; iNode++){

    	double r = nodesFine_[iNode] -> getDistFunction();
        
		if (r >= glueZoneThickness){
        	wFuncValue = 1. - arlequinEpsilon;
    	} else {
        	wFuncValue = (1. - arlequinEpsilon) / glueZoneThickness * r;
        	if (wFuncValue > (1. - arlequinEpsilon)) wFuncValue = 1. - arlequinEpsilon;

    	}; 
              
        nodesFine_[iNode] -> setWeightFunction(wFuncValue);

    };  //ielem


    for (int jel = 0; jel < numElemFine; jel++){
        elementsFine_[jel] -> setIntegPointWeightFunction_ISO();        
    };  
     
    return;

};


//------------------------------------------------------------------------------
//------------------------PRINT RESULTS IN PARAVIEW-----------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::printResults(int step) {

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (rank == 0){

        //PRINT COARSE MODEL RESULTS 
        std::string result;
        std::ostringstream convert;
        convert << step+100000;
        result = convert.str();

        int dim = 2;
        
        std::string s = "COARSEoutput"+result+".vtu";
        std::fstream output_v(s.c_str(), std::ios_base::out);

        // if (elemType == 0){ //FEM mesh

        // 	output_v << "<?xml version=\"1.0\"?>" << std::endl
        //          	<< "<VTKFile type=\"UnstructuredGrid\">" << std::endl
        //          	<< "  <UnstructuredGrid>" << std::endl
        //          	<< "  <Piece NumberOfPoints=\"" << numNodesCoarse
        //          	<< "\"  NumberOfCells=\"" << numElemCoarse
        //          	<< "\">" << std::endl;

	       //  //WRITE NODAL COORDINATES
	       //  output_v << "    <Points>" << std::endl
	       //           << "      <DataArray type=\"Float64\" "
	       //           << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

	       //  for (int i = 0; i < numNodesCoarse; i++){
	       //      double *x = nodesCoarse_[i] -> getCoordinates();
	       //      double x_[dim];
	       //      for (int i = 0; i < dim; i++) x_[i] = x[i];
	       //      output_v << x_[0] << " " << x_[1] << " " << -0.01 << std::endl;
	       //  };

	       //  output_v << "      </DataArray>" << std::endl
	       //           << "    </Points>" << std::endl;
        
	       //  //WRITE ELEMENT CONNECTIVITY
	       //  output_v << "    <Cells>" << std::endl
	       //           << "      <DataArray type=\"Int32\" "
	       //           << "Name=\"connectivity\" format=\"ascii\">" << std::endl;
	        
	       //  for (int i = 0; i < numElemCoarse; i++){
	       //      int *connec = elementsCoarse_[i] -> getConnectivity();
	       //      int con[6];
	       //      for (int i = 0; i < 6; i++) con[i] = connec[i];
	       //      output_v << con[0] << " " << con[1] << " " << con[2] << " " 
	       //               << con[3] << " " << con[4] << " " << con[5] << std::endl;
	       //  };
	       //  output_v << "      </DataArray>" << std::endl;
	      
	       //  //WRITE OFFSETS IN DATA ARRAY
	       //  output_v << "      <DataArray type=\"Int32\""
	       //           << " Name=\"offsets\" format=\"ascii\">" << std::endl;
	        
	       //  int aux = 0;
	       //  for (int i = 0; i < numElemCoarse; i++){
	       //      output_v << aux + 6 << std::endl;
	       //      aux += 6;
	       //  };
	       //  output_v << "      </DataArray>" << std::endl;
	      
	       //  //WRITE ELEMENT TYPES
	       //  output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
	       //           << "format=\"ascii\">" << std::endl;
        
	       //  for (int i = 0; i < numElemCoarse; i++){
	       //      output_v << 22 << std::endl;
	       //  };

	       //  output_v << "      </DataArray>" << std::endl
	       //           << "    </Cells>" << std::endl;

	       //  //WRITE NODAL RESULTS
	       //  output_v << "    <PointData>" << std::endl;

	       //  if (coarseModel.printVelocity){
	       //      output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
	       //               << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
	       //      for (int i=0; i<numNodesCoarse; i++){ 
	       //      	     output_v << nodesCoarse_[i] -> getVelocity(0) << " "              
	       //               << nodesCoarse_[i] -> getVelocity(1)  << " " 
	       //               << 0. << std::endl;
	       //      };
	       //      output_v << "      </DataArray> " << std::endl;
	       //  };
       
	       //  if (coarseModel.printDistFunction){
	       //      output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
	       //               << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;
	       //      for (int i=0; i<numNodesCoarse; i++){        
	       //          output_v << nodesCoarse_[i] -> getDistFunction() << std::endl;
	       //      };
	       //      output_v << "      </DataArray> " << std::endl;
	       //  };
        
	       //  if (coarseModel.printPressure){
	       //      output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
	       //               << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
	       //      for (int i=0; i<numNodesCoarse; i++){
	       //      	output_v << 0. << " " << 0. << " " 
	       //          << nodesCoarse_[i] -> getPressure()<< std::endl;
	       //      };
	       //      output_v << "      </DataArray> " << std::endl;
	       //  };

	       //  output_v << "    </PointData>" << std::endl; 

	       //  //WRITE ELEMENT RESULTS
	       //  output_v << "    <CellData>" << std::endl;
	        
	       //  if (coarseModel.printProcess){
	       //      output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
	       //               << "Name=\"Process\" format=\"ascii\">" << std::endl;
	       //      for (int i=0; i<numElemCoarse; i++){
	       //          output_v << domDecompCoarse.first[i] << std::endl;
	       //      };
	       //      output_v << "      </DataArray> " << std::endl;
	       //  };

	       //  int cont=0;
	       //  if (coarseModel.printGlueZone){
	       //      output_v <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
	       //               << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
	       //      cont = 0;
	       //      for (int i=0; i<numElemCoarse; i++){
	       //          if (elementsGlueZoneCoarse_[cont] == i){
	       //              output_v << 1.0 << std::endl;
	       //              cont++;
	       //          }else{
	       //              output_v << 0.0 << std::endl;
	       //          };
	       //      };
	       //      output_v << "      </DataArray> " << std::endl;
	       //  };

	       //  output_v << "    </CellData>" << std::endl; 

	       //  //FINALIZE OUTPUT FILE
	        
	       //  output_v << "  </Piece>" << std::endl
	       //     << "  </UnstructuredGrid>" << std::endl
	       //     << "</VTKFile>" << std::endl;
        

        // } else { //IGA mesh

	        int numBezierNodes = coarseModel.NumBezierNodes;

	        output_v << "<?xml version=\"1.0\"?>" << std::endl
	                 << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
	                 << "  <UnstructuredGrid>" << std::endl
	                 << "  <Piece NumberOfPoints=\"" << numBezierNodes
	                 << "\"  NumberOfCells=\"" << numElemCoarse
	                 << "\">" << std::endl;

	        //WRITE NODAL COORDINATES
	        output_v << "    <Points>" << std::endl
	                 << "      <DataArray type=\"Float64\" "
	                 << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

	        // Bezier Extraction
	        double Vel[numBezierNodes][dim], realVel[numBezierNodes][dim], Coord[numBezierNodes][dim];
	        double Press[numBezierNodes], realPress[numBezierNodes], Distance[numBezierNodes], EnergyW[numBezierNodes];

            //Interpolated variables 
            double xsiCP[9][dim];
            xsiCP[0][0] = -1.; xsiCP[0][1] = -1.;
            xsiCP[1][0] = 0.; xsiCP[1][1] = -1.;
            xsiCP[2][0] = 1.; xsiCP[2][1] = -1.;
            xsiCP[3][0] = -1.; xsiCP[3][1] = 0.;
            xsiCP[4][0] = 0.; xsiCP[4][1] = 0.;
            xsiCP[5][0] = 1.; xsiCP[5][1] = 0.;
            xsiCP[6][0] = -1.; xsiCP[6][1] = 1.;
            xsiCP[7][0] = 0.; xsiCP[7][1] = 1.;
            xsiCP[8][0] = 1.; xsiCP[8][1] = 1.;

	        for (int iElem = 0; iElem < numElemCoarse; iElem++){
	                            
	            int *Beconnec = elementsCoarse_[iElem] -> getBezierConnectivity();
	            int *connec = elementsCoarse_[iElem] -> getConnectivity();
		
				double coord_[9][dim],vel_[9][dim], realvel_[9][dim];
	            double press_[9], realpress_[9], distance_[9], energyW_[9];
	        	
                //Data in the NURBS control points
	            for (int i = 0; i < 9; i++){
					double *x = nodesCoarse_[connec[i]] -> getCoordinates();
	                for (int j = 0; j < dim; j++){
	                	coord_[i][j] = x[j];
	                	vel_[i][j] = nodesCoarse_[connec[i]] -> getVelocity(j);
	                	realvel_[i][j] = nodesCoarse_[connec[i]] -> getVelocityArlequin(j);
	                };

                    press_[i] = nodesCoarse_[connec[i]] -> getPressure();
                    realpress_[i] = nodesCoarse_[connec[i]] -> getPressureArlequin();
	                distance_[i] = nodesCoarse_[connec[i]] -> getDistFunction(); 
	                energyW_[i] = nodesCoarse_[connec[i]] -> getWeightFunction(); 
	         	};
    
                //interpolated values (Bézier variables)
                double Bcoord_[9][2] = {};
                double Bvel_[9][2]= {};
                double Brealvel_[9][2]= {};
                double Bpress_[9]= {};
                double Brealpress_[9]= {};
                double Bdistance_[9]= {};
                double BenergyW_[9] = {};

                for (int i = 0; i < 9; i++){

                    QuadShapeFunction<2> shapeQuad;
                    double phi_[9],wpc[9],xsi[dim];
                    for (int k = 0; k < 9; k ++) wpc[k] = nodesCoarse_[connec[k]] -> getWeightPC();  
                    int *inc_ = nodesCoarse_[connec[8]] -> getINC(); 
                    int patch = elementsCoarse_[iElem] -> getPatch();
                    for (int j = 0; j < dim; j++) xsi[j] = xsiCP[i][j]; 
                    shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,IsoParCoarse,patch);

                    for (int j = 0; j < 9; j++){
                        for (int k = 0; k < dim; k++){
                        	Bcoord_[i][k] += phi_[j] * coord_[j][k];
                        	Bvel_[i][k] += phi_[j] * vel_[j][k];
                        	Brealvel_[i][k] += phi_[j] * realvel_[j][k];
                        };
                        Bpress_[i] += phi_[j] * press_[j];
                        Brealpress_[i] += phi_[j] * realpress_[j];
                        Bdistance_[i] += phi_[j] * distance_[j];
                        BenergyW_[i] += phi_[j] * energyW_[j];
                    };

                };

	            for (int i = 0; i< 9; i++){
	                for (int j = 0; j < dim; j++){
	                	Coord[Beconnec[i]][j] = Bcoord_[i][j];
	                	Vel[Beconnec[i]][j] = Bvel_[i][j];
	                	realVel[Beconnec[i]][j] = Brealvel_[i][j];
	                };
	                Press[Beconnec[i]] = Bpress_[i];
	                realPress[Beconnec[i]] = Brealpress_[i];
	                Distance[Beconnec[i]] = Bdistance_[i];
	                EnergyW[Beconnec[i]] = BenergyW_[i];
	            };
	        
	        };//iElem

	                
	        for (int i = 0; i< numBezierNodes;i++){
	        	output_v << Coord[i][0] << " " << Coord[i][1] << " " << 0. << std::endl;
	        }

	        output_v << "      </DataArray>" << std::endl
	                 << "    </Points>" << std::endl;
	        
	        //WRITE ELEMENT CONNECTIVITY
	        output_v << "    <Cells>" << std::endl
	                 << "      <DataArray type=\"Int32\" "
	                 << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

	        for (int iElem = 0; iElem < numElemCoarse; ++iElem){
	           
	            int *Bconnec_ = elementsCoarse_[iElem]->getBezierConnectivity();  
	            int Bconnec[9];
	            for (int i = 0; i < 9; i++) Bconnec[i] = Bconnec_[i];            
	            
	        	output_v << Bconnec[0] << " " << Bconnec[2] << " " << Bconnec[8] << " "
	            << Bconnec[6] << " " << Bconnec[1] << " " << Bconnec[5] << " " 
	            << Bconnec[7] << " " << Bconnec[3] << " " << Bconnec[4]<<  std::endl;

	        };
	        output_v << "      </DataArray>" << std::endl;
	      
	        //WRITE OFFSETS IN DATA ARRAY
	        output_v << "      <DataArray type=\"Int32\""
	                 << " Name=\"offsets\" format=\"ascii\">" << std::endl;
	        int aux = 0;
	        for (int i = 0; i < numElemCoarse; i++){
	            output_v << aux + 9 << std::endl;
	            aux += 9;
	        };
	        output_v << "      </DataArray>" << std::endl;
	      
	        //WRITE ELEMENT TYPES
	        output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
	                 << "format=\"ascii\">" << std::endl;
	        
	        for (int i = 0; i < numElemCoarse; i++){
	            output_v << 70 << std::endl;
	        };
	        output_v << "      </DataArray>" << std::endl
	                 << "    </Cells>" << std::endl;

	        // WRITE NODAL RESULTS
	        output_v << "    <PointData>" << std::endl;

	        if (coarseModel.printVelocity){
	            output_v<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
	                      << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
	            for (int i=0; i<numBezierNodes; i++){
	                output_v << Vel[i][0] << " "              
	                          << Vel[i][1] << " " 
	                          << 0. << std::endl;
	            };
	            output_v << "      </DataArray> " << std::endl;
	        };

	        if (coarseModel.printRealVelocity){
	            output_v<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
	                      << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
	            for (int i=0; i<numBezierNodes; i++){
	                output_v << realVel[i][0] << " "              
	                          << realVel[i][1] << " " 
	                          << 0. << std::endl;
	            };
	            output_v << "      </DataArray> " << std::endl;
	        };


	        if (coarseModel.printPressure){
	            output_v<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
	                     << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
	            for (int i=0; i<numBezierNodes; i++){
	                output_v << 0. << " " << 0. << " " 
	                          << Press[i] << std::endl;
	            };
	            output_v << "      </DataArray> " << std::endl;
	        };

	        if (coarseModel.printRealPressure){
	            output_v<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
	                     << "Name=\"Real Pressure\" format=\"ascii\">" << std::endl;
	            for (int i=0; i<numBezierNodes; i++){
	                output_v << 0. << " " << 0. << " " 
	                          << realPress[i] << std::endl;
	            };
	            output_v << "      </DataArray> " << std::endl;
	        };


	        if (coarseModel.printDistFunction){
	            output_v<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
	                     << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;
	            for (int i=0; i<numBezierNodes; i++){
	                output_v << Distance[i] << std::endl;
	            };
	            output_v << "      </DataArray> " << std::endl;
	        };

	        if (coarseModel.printEnergyWeightFunction){
	            output_v<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
	                     << "Name=\"Energy Weight Function\" format=\"ascii\">" << std::endl;
	            for (int i=0; i<numBezierNodes; i++){
	                output_v << EnergyW[i] << std::endl;
	            };
	            output_v << "      </DataArray> " << std::endl;
	        };



	        output_v << "    </PointData>" << std::endl; 

	        //WRITE ELEMENT RESULTS
	        output_v << "    <CellData>" << std::endl;
	        
	        if (coarseModel.printProcess){
	            output_v<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
	                     << "Name=\"Process\" format=\"ascii\">" << std::endl;
	            for (int i=0; i<numElemCoarse; i++){
	                output_v << domDecompFine.first[i] << std::endl;
	            };
	            output_v << "      </DataArray> " << std::endl;
	        };


	        if (fineModel.printGlueZone){
	            output_v<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
	                     << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
	            int cont=0;
	            for (int i=0; i<numElemCoarse; i++){
	                if (elementsGlueZoneCoarse_[cont] == i){
	                    output_v << 1.0 << std::endl;
	                    cont += 1; 
	                }else{
	                    output_v << 0.0 << std::endl;
	                };
	            };
	            output_v << "      </DataArray> " << std::endl;
	        };

	        output_v << "    </CellData>" << std::endl; 

	        //FINALIZE OUTPUT FILE
	        output_v << "  </Piece>" << std::endl
	               << "  </UnstructuredGrid>" << std::endl
	               << "</VTKFile>" << std::endl;

        // };


    
        //PRINT FINE MODEL RESULTS - IGA MESH
        numBezierNodes = fineModel.NumBezierNodes;

        std::string f = "FINEoutput"+result+".vtu";
        
        std::fstream output_vf(f.c_str(), std::ios_base::out);

        output_vf << "<?xml version=\"1.0\"?>" << std::endl
                 << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                 << "  <UnstructuredGrid>" << std::endl
                 << "  <Piece NumberOfPoints=\"" << numBezierNodes
                 << "\"  NumberOfCells=\"" << numElemFine
                 << "\">" << std::endl;

        //WRITE NODAL COORDINATES
        output_vf << "    <Points>" << std::endl
                 << "      <DataArray type=\"Float64\" "
                 << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;


        // Bezier Extraction
        double VelF[numBezierNodes][dim], realVelF[numBezierNodes][dim], CoordF[numBezierNodes][dim], LagMult[numBezierNodes][dim],
        PressF[numBezierNodes], realPressF[numBezierNodes], DistanceF[numBezierNodes],EnergyWF[numBezierNodes];

        for (int iElem = 0; iElem < numElemFine; ++iElem){
                         
            int *Beconnec = elementsFine_[iElem] -> getBezierConnectivity();
            int *connec = elementsFine_[iElem] -> getConnectivity();

        	// Data in the NURBS control points
            double coord_[9][dim],vel_[9][dim],realvel_[9][dim],lagmult_[9][dim],
            	   press_[9],realpress_[9],distance_[9],energyW_[9];
            
            for (int i = 0; i < 9; i++){
				double *x = nodesFine_[connec[i]] -> getCoordinates();
                for (int j = 0; j < dim; j++){
                	coord_[i][j] = x[j];
                	vel_[i][j] = nodesFine_[connec[i]] -> getVelocity(j);
                	realvel_[i][j] = nodesFine_[connec[i]] -> getVelocityArlequin(j);
                	lagmult_[i][j] = nodesFine_[connec[i]] -> getLagrangeMultiplier(j);
                };
                press_[i] = nodesFine_[connec[i]] -> getPressure();
                realpress_[i] = nodesFine_[connec[i]] -> getPressureArlequin();
                distance_[i] = nodesFine_[connec[i]] -> getDistFunction();   
                energyW_[i] = nodesFine_[connec[i]] -> getWeightFunction();              
         	};  

            //interpolated values (Bezier variables)
            double Bcoord_[9][2] = {};
            double Bvel_[9][2]= {};
            double Brealvel_[9][2]= {};
            double Bpress_[9]= {};
            double Brealpress_[9]= {};
            double Bdistance_[9]= {};
            double BenergyW_[9] = {};
            double Blagmult_[9][2] = {};

            for (int i = 0; i < 9; i++){

                QuadShapeFunction<2> shapeQuad;
                double phi_[9],wpc[9],xsi[dim];
                for (int k = 0; k < 9; k ++) wpc[k] = nodesFine_[connec[k]] -> getWeightPC();  
                int *inc_ = nodesFine_[connec[8]] -> getINC(); 
                int patch = elementsFine_[iElem] -> getPatch();
                for (int j = 0; j < dim; j++) xsi[j] = xsiCP[i][j]; 
                shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,IsoParFine,patch);

                for (int j = 0; j < 9; j++){
                	for (int k = 0; k < dim; k++){
                		Bcoord_[i][k] += phi_[j] * coord_[j][k];
                    	Bvel_[i][k] += phi_[j] * vel_[j][k];
                    	Brealvel_[i][k] += phi_[j] * realvel_[j][k];
                    	Blagmult_[i][k] += phi_[j] * lagmult_[j][k];
                	};
                    Bpress_[i] += phi_[j] * press_[j];
                    Brealpress_[i] += phi_[j] * realpress_[j];
                    Bdistance_[i] += phi_[j] * distance_[j];
                    BenergyW_[i] += phi_[j] * energyW_[j];
                };
            };

            for (int i = 0; i< 9; i++){
                for (int j = 0; j < dim; j++){
                	CoordF[Beconnec[i]][j] = Bcoord_[i][j];
                	VelF[Beconnec[i]][j] = Bvel_[i][j];
                	realVelF[Beconnec[i]][j] = Brealvel_[i][j];
                	LagMult[Beconnec[i]][j] = Blagmult_[i][j];
                };
                PressF[Beconnec[i]] = Bpress_[i];
                realPressF[Beconnec[i]] = Brealpress_[i];
                DistanceF[Beconnec[i]] = Bdistance_[i];
                EnergyWF[Beconnec[i]] = BenergyW_[i];
	       };

        };

              
        for (int i = 0; i< numBezierNodes;i++){
        	output_vf << CoordF[i][0] << " " << CoordF[i][1] << " " << 0.1 << std::endl;
        };

        output_vf << "      </DataArray>" << std::endl
                 << "    </Points>" << std::endl;
        
        //WRITE ELEMENT CONNECTIVITY
        output_vf << "    <Cells>" << std::endl
                 << "      <DataArray type=\"Int32\" "
                 << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

        for (int iElem = 0; iElem < numElemFine; ++iElem){
           
            int *Bconnec_= elementsFine_[iElem]->getBezierConnectivity();  
            int Bconnec[9];
            for (int i = 0; i < 9; i++) Bconnec[i] = Bconnec_[i];            
            
        	output_vf << Bconnec[0] << " " << Bconnec[2] << " " << Bconnec[8] << " "
            << Bconnec[6] << " " << Bconnec[1] << " " << Bconnec[5] << " " 
            << Bconnec[7] << " " << Bconnec[3] << " " << Bconnec[4]<<  std::endl;

        };

        output_vf << "      </DataArray>" << std::endl;
      
        //WRITE OFFSETS IN DATA ARRAY
        output_vf << "      <DataArray type=\"Int32\""
                 << " Name=\"offsets\" format=\"ascii\">" << std::endl;
        
        aux = 0;
        for (int i = 0; i < numElemFine; i++){
            output_vf << aux + 9 << std::endl;
            aux += 9;
        };
        output_vf << "      </DataArray>" << std::endl;
      
        //WRITE ELEMENT TYPES
        output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                 << "format=\"ascii\">" << std::endl;
        for (int i = 0; i < numElemFine; i++){
            output_vf << 70 << std::endl;
        };
        output_vf << "      </DataArray>" << std::endl
                 << "    </Cells>" << std::endl;

        // WRITE NODAL RESULTS
        output_vf << "    <PointData>" << std::endl;

        if (fineModel.printVelocity){
            output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                      << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
            for (int i=0; i<numBezierNodes; i++){
                output_vf << VelF[i][0] << " "              
                          << VelF[i][1] << " " 
                          << 0. << std::endl;
            };
            output_vf << "      </DataArray> " << std::endl;
        };


        if (fineModel. printRealVelocity){
            output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                      << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
            for (int i=0; i<numBezierNodes; i++){
                output_vf << realVelF[i][0] << " "              
                          << realVelF[i][1] << " " 
                          << 0. << std::endl;
            };
            output_vf << "      </DataArray> " << std::endl;
        };

        if (fineModel.printPressure){
            output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                     << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
            for (int i=0; i<numBezierNodes; i++){
                output_vf << 0. << " " << 0. << " " 
                          << PressF[i] << std::endl;
            };
            output_vf << "      </DataArray> " << std::endl;
        };

         if (fineModel.printRealPressure){
            output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                     << "Name=\"Real Pressure\" format=\"ascii\">" << std::endl;
            for (int i=0; i<numBezierNodes; i++){
                output_vf << 0. << " " << 0. << " " 
                          << realPressF[i] << std::endl;
            };
            output_vf << "      </DataArray> " << std::endl;
        };

        if (fineModel.printLagrangeMultipliers){
            output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                     << "Name=\"Lagrange Multipliers\" format=\"ascii\">" << std::endl;
            for (int i=0; i<numBezierNodes; i++){
                output_vf << LagMult[i][0] <<  " "
                          << LagMult[i][1] << " " 
                          << 0.0 << " " << std::endl;
            };
            output_vf << "      </DataArray> " << std::endl;
        };


        if (fineModel.printDistFunction){
            output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                     << "Name=\"printDistFunction\" format=\"ascii\">" << std::endl;
            for (int i=0; i<numBezierNodes; i++){
                output_vf << DistanceF[i] << std::endl;
            };
            output_vf << "      </DataArray> " << std::endl;
        };


        if (fineModel.printEnergyWeightFunction){
            output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                     << "Name=\"EnergyWeightFunction\" format=\"ascii\">" << std::endl;
            for (int i=0; i<numBezierNodes; i++){
                output_vf << EnergyWF[i] << std::endl;
            };
            output_vf << "      </DataArray> " << std::endl;
        };


        if (fineModel.printNodalCorrespondence){
            output_vf<<"      <DataArray type=\"Int32\" NumberOfComponents=\"1\" "
                     << "Name=\"printNodalCorrespondece\" format=\"ascii\">" << std::endl;
            for (int i=0; i<numBezierNodes; i++){
                output_vf << 0 << std::endl;
            };
            output_vf << "      </DataArray> " << std::endl;
        };

        output_vf << "    </PointData>" << std::endl; 

        //WRITE ELEMENT RESULTS
        output_vf << "    <CellData>" << std::endl;
        
        if (fineModel.printProcess){
            output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                     << "Name=\"Process\" format=\"ascii\">" << std::endl;
            for (int i=0; i<numElemFine; i++){
                output_vf << domDecompFine.first[i] << std::endl;
            };
            output_vf << "      </DataArray> " << std::endl;
        };


        if (fineModel.printGlueZone){
            output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                     << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
            int cont=0;
            for (int i=0; i<numElemFine; i++){
                if (elementsGlueZoneFine_[cont] == i){
                    output_vf << 1.0 << std::endl;
                    cont += 1; 
                }else{
                    output_vf << 0.0 << std::endl;
                };
            };
            output_vf << "      </DataArray> " << std::endl;
        };

        output_vf << "    </CellData>" << std::endl; 

        //FINALIZE OUTPUT FILE
        output_vf << "  </Piece>" << std::endl
               << "  </UnstructuredGrid>" << std::endl
               << "</VTKFile>" << std::endl;














    //           //LAGRANGE MULTIPLIERS
    //     std::string f = "LAGRANGEoutput"+result+".vtu";
        
    //     std::fstream output_vf(f.c_str(), std::ios_base::out);

    //     output_vf << "<?xml version=\"1.0\"?>" << std::endl
    //              << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
    //              << "  <UnstructuredGrid>" << std::endl
    //              << "  <Piece NumberOfPoints=\"" << numNodesGlueZoneFine
    //              << "\"  NumberOfCells=\"" << numElemGlueZoneFine
    //              << "\">" << std::endl;

    //     //WRITE NODAL COORDINATES
    //     output_vf << "    <Points>" << std::endl
    //              << "      <DataArray type=\"Float64\" "
    //              << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;


    //     // Bezier Extraction
    //     double CoordLA[numBezierNodes][dim];
    //     double PressF[numBezierNodes], realPressF[numBezierNodes], DistanceF[numBezierNodes],EnergyWF[numBezierNodes];

    //     for (int iElem = 0; iElem < numElemFine; ++iElem){
                         
    //         int *Beconnec = elementsFine_[iElem] -> getBezierConnectivity();
    //         int *connec = elementsFine_[iElem] -> getConnectivity();

    //     	// Data in the NURBS control points
    //         double coord_[9][dim],vel_[9][dim],realvel_[9][dim],
    //         	   press_[9],realpress_[9],distance_[9],energyW_[9];
            
    //         for (int i = 0; i < 9; i++){
				// double *x = nodesFine_[connec[i]] -> getCoordinates();
    //             for (int j = 0; j < dim; j++){
    //             	coord_[i][j] = x[j];
    //             	vel_[i][j] = nodesFine_[connec[i]] -> getVelocity(j);
    //             	realvel_[i][j] = nodesFine_[connec[i]] -> getVelocityArlequin(j);
    //             };
    //             press_[i] = nodesFine_[connec[i]] -> getPressure();
    //             realpress_[i] = nodesFine_[connec[i]] -> getPressureArlequin();
    //             distance_[i] = nodesFine_[connec[i]] -> getDistFunction();   
    //             energyW_[i] = nodesFine_[connec[i]] -> getWeightFunction();              
    //      	};  

    //         //interpolated values (Bezier variables)
    //         double Bcoord_[9][2] = {};
    //         double Bvel_[9][2]= {};
    //         double Brealvel_[9][2]= {};
    //         double Bpress_[9]= {};
    //         double Brealpress_[9]= {};
    //         double Bdistance_[9]= {};
    //         double BenergyW_[9] = {};

    //         for (int i = 0; i < 9; i++){

    //             QuadShapeFunction<2> shapeQuad;
    //             double phi_[9],wpc[9],xsi[dim];
    //             for (int k = 0; k < 9; k ++) wpc[k] = nodesFine_[connec[k]] -> getWeightPC();  
    //             int *inc_ = nodesFine_[connec[8]] -> getINC(); 
    //             int patch = elementsFine_[iElem] -> getPatch();
    //             for (int j = 0; j < dim; j++) xsi[j] = xsiCP[i][j]; 
    //             shapeQuad.evaluateIso(xsi,phi_,wpc,inc_,IsoParFine,patch);

    //             for (int j = 0; j < 9; j++){
    //             	for (int k = 0; k < dim; k++){
    //             		Bcoord_[i][k] += phi_[j] * coord_[j][k];
    //                 	Bvel_[i][k] += phi_[j] * vel_[j][k];
    //                 	Brealvel_[i][k] += phi_[j] * realvel_[j][k];
    //             	};
    //                 Bpress_[i] += phi_[j] * press_[j];
    //                 Brealpress_[i] += phi_[j] * realpress_[j];
    //                 Bdistance_[i] += phi_[j] * distance_[j];
    //                 BenergyW_[i] += phi_[j] * energyW_[j];
    //             };
    //         };

    //         for (int i = 0; i< 9; i++){
    //             for (int j = 0; j < dim; j++){
    //             	CoordF[Beconnec[i]][j] = Bcoord_[i][j];
    //             	VelF[Beconnec[i]][j] = Bvel_[i][j];
    //             	realVelF[Beconnec[i]][j] = Brealvel_[i][j];
    //             };
    //             PressF[Beconnec[i]] = Bpress_[i];
    //             realPressF[Beconnec[i]] = Brealpress_[i];
    //             DistanceF[Beconnec[i]] = Bdistance_[i];
    //             EnergyWF[Beconnec[i]] = BenergyW_[i];
	   //     };

    //     };

              
    //     for (int i = 0; i< numBezierNodes;i++){
    //     	output_vf << CoordF[i][0] << " " << CoordF[i][1] << " " << 0. << std::endl;
    //     };

    //     output_vf << "      </DataArray>" << std::endl
    //              << "    </Points>" << std::endl;
        
    //     //WRITE ELEMENT CONNECTIVITY
    //     output_vf << "    <Cells>" << std::endl
    //              << "      <DataArray type=\"Int32\" "
    //              << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

    //     for (int iElem = 0; iElem < numElemFine; ++iElem){
           
    //         int *Bconnec_= elementsFine_[iElem]->getBezierConnectivity();  
    //         int Bconnec[9];
    //         for (int i = 0; i < 9; i++) Bconnec[i] = Bconnec_[i];            
            
    //     	output_vf << Bconnec[0] << " " << Bconnec[2] << " " << Bconnec[8] << " "
    //         << Bconnec[6] << " " << Bconnec[1] << " " << Bconnec[5] << " " 
    //         << Bconnec[7] << " " << Bconnec[3] << " " << Bconnec[4]<<  std::endl;

    //     };

    //     output_vf << "      </DataArray>" << std::endl;
      
    //     //WRITE OFFSETS IN DATA ARRAY
    //     output_vf << "      <DataArray type=\"Int32\""
    //              << " Name=\"offsets\" format=\"ascii\">" << std::endl;
        
    //     aux = 0;
    //     for (int i = 0; i < numElemFine; i++){
    //         output_vf << aux + 9 << std::endl;
    //         aux += 9;
    //     };
    //     output_vf << "      </DataArray>" << std::endl;
      
    //     //WRITE ELEMENT TYPES
    //     output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
    //              << "format=\"ascii\">" << std::endl;
    //     for (int i = 0; i < numElemFine; i++){
    //         output_vf << 70 << std::endl;
    //     };
    //     output_vf << "      </DataArray>" << std::endl
    //              << "    </Cells>" << std::endl;

    //     // WRITE NODAL RESULTS
    //     output_vf << "    <PointData>" << std::endl;

    //     if (fineModel.printVelocity){
    //         output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
    //                   << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
    //         for (int i=0; i<numBezierNodes; i++){
    //             output_vf << VelF[i][0] << " "              
    //                       << VelF[i][1] << " " 
    //                       << 0. << std::endl;
    //         };
    //         output_vf << "      </DataArray> " << std::endl;
    //     };


    //     if (fineModel. printRealVelocity){
    //         output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
    //                   << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
    //         for (int i=0; i<numBezierNodes; i++){
    //             output_vf << realVelF[i][0] << " "              
    //                       << realVelF[i][1] << " " 
    //                       << 0. << std::endl;
    //         };
    //         output_vf << "      </DataArray> " << std::endl;
    //     };

    //     if (fineModel.printPressure){
    //         output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
    //                  << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
    //         for (int i=0; i<numBezierNodes; i++){
    //             output_vf << 0. << " " << 0. << " " 
    //                       << PressF[i] << std::endl;
    //         };
    //         output_vf << "      </DataArray> " << std::endl;
    //     };

    //      if (fineModel.printRealPressure){
    //         output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
    //                  << "Name=\"Real Pressure\" format=\"ascii\">" << std::endl;
    //         for (int i=0; i<numBezierNodes; i++){
    //             output_vf << 0. << " " << 0. << " " 
    //                       << realPressF[i] << std::endl;
    //         };
    //         output_vf << "      </DataArray> " << std::endl;
    //     };


    //     if (fineModel.printDistFunction){
    //         output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
    //                  << "Name=\"printDistFunction\" format=\"ascii\">" << std::endl;
    //         for (int i=0; i<numBezierNodes; i++){
    //             output_vf << DistanceF[i] << std::endl;
    //         };
    //         output_vf << "      </DataArray> " << std::endl;
    //     };


    //     if (fineModel.printEnergyWeightFunction){
    //         output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
    //                  << "Name=\"EnergyWeightFunction\" format=\"ascii\">" << std::endl;
    //         for (int i=0; i<numBezierNodes; i++){
    //             output_vf << EnergyWF[i] << std::endl;
    //         };
    //         output_vf << "      </DataArray> " << std::endl;
    //     };

    //     if (fineModel.printNodalCorrespondence){
    //         output_vf<<"      <DataArray type=\"Int32\" NumberOfComponents=\"1\" "
    //                  << "Name=\"printNodalCorrespondece\" format=\"ascii\">" << std::endl;
    //         for (int i=0; i<numBezierNodes; i++){
    //             output_vf << 0 << std::endl;
    //         };
    //         output_vf << "      </DataArray> " << std::endl;
    //     };

    //     output_vf << "    </PointData>" << std::endl; 

    //     //WRITE ELEMENT RESULTS
    //     output_vf << "    <CellData>" << std::endl;
        
    //     if (fineModel.printProcess){
    //         output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
    //                  << "Name=\"Process\" format=\"ascii\">" << std::endl;
    //         for (int i=0; i<numElemFine; i++){
    //             output_vf << domDecompFine.first[i] << std::endl;
    //         };
    //         output_vf << "      </DataArray> " << std::endl;
    //     };


    //     if (fineModel.printGlueZone){
    //         output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
    //                  << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
    //         int cont=0;
    //         for (int i=0; i<numElemFine; i++){
    //             if (elementsGlueZoneFine_[cont] == i){
    //                 output_vf << 1.0 << std::endl;
    //                 cont += 1; 
    //             }else{
    //                 output_vf << 0.0 << std::endl;
    //             };
    //         };
    //         output_vf << "      </DataArray> " << std::endl;
    //     };

    //     output_vf << "    </CellData>" << std::endl; 

    //     //FINALIZE OUTPUT FILE
    //     output_vf << "  </Piece>" << std::endl
    //            << "  </UnstructuredGrid>" << std::endl
    //            << "</VTKFile>" << std::endl;

    };
};



template<>
void Arlequin<2>::printResultsIP(int step) {

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    //Coarse Mesh (IGA MESH)
    if (rank == 0){
        std::string result;
        std::ostringstream convert;

        convert << step+100000;
        result = convert.str();
        std::string s = "saidaCoarseIP"+result+".vtu";
        
        std::fstream output_v(s.c_str(), std::ios_base::out);

        int numberIntPoints = elementsFine_[0] -> getNumberOfIntegrationPointsSpecial();

    	output_v << "<?xml version=\"1.0\"?>" << std::endl
                 << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                 << "  <UnstructuredGrid>" << std::endl
                 << "  <Piece NumberOfPoints=\"" << numElemGlueZoneFine*numberIntPoints
                 << "\"  NumberOfCells=\"" << numElemGlueZoneFine*numberIntPoints
                 << "\">" << std::endl;

        output_v << "    <Points>" << std::endl
                 << "      <DataArray type=\"Float64\" "
                 << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

        for (int i = 0; i < numElemGlueZoneFine; i++){

        	

        	for (int ip = 0; ip < numberIntPoints; ip++){
            	
            	QuadShapeFunction<2>  shapeQuad;
            	
            	//correspondent integration point in the coarse mesh
            	double qxsiC[2];
            	qxsiC[0] = elementsFine_[elementsGlueZoneFine_[i]] -> getIntegPointCoordinatesValue(ip,0);
            	qxsiC[1] = elementsFine_[elementsGlueZoneFine_[i]] -> getIntegPointCoordinatesValue(ip,1);
            	
            	//coarse element index
            	int indCoarseElem = elementsFine_[elementsGlueZoneFine_[i]] -> getIntegPointCorrespondenceElement(ip);
            	
                int *connec = elementsCoarse_[indCoarseElem] -> getConnectivity();


                //Computes nurbs basis functions
                int patch = elementsCoarse_[indCoarseElem] -> getPatch();
        		int *inc = nodesCoarse_[connec[8]] -> getINC();
        		double wpc[9],phi_[9];
        		for (int k = 0; k<9; k++) wpc[k] = nodesCoarse_[connec[k]] -> getWeightPC();
                shapeQuad.evaluateIso(qxsiC,phi_,wpc,inc,IsoParCoarse,patch);
        			            		
        		double coord[2] = {};
        		for (int j = 0; j < 9; j++){
        			double *x = nodesCoarse_[connec[j]] -> getCoordinates();
        			coord[0] += x[0]*phi_[j];
        			coord[1] += x[1]*phi_[j];
        		}

        		output_v << coord[0] << " " << coord[1] << " " << 0.0 << std::endl;
            	
        	} //loop integration points
        } // loop glue fine mesh elements   

        output_v << "      </DataArray>" << std::endl
                 << "    </Points>" << std::endl;

        

        
      
        //WRITE ELEMENT CONNECTIVITY
        output_v << "    <Cells>" << std::endl
                << "      <DataArray type=\"Int32\" "
                << "Name=\"connectivity\" format=\"ascii\">" << std::endl;
        
        for (int numN = 0; numN < numElemGlueZoneFine*numberIntPoints; ++numN){
                output_v << numN << std::endl;
        }

        output_v << "      </DataArray>" << std::endl;
  
        //WRITE OFFSETS IN DATA ARRAY
        output_v << "      <DataArray type=\"Int32\""
                << " Name=\"offsets\" format=\"ascii\">" << std::endl;
    
        int aux = 0;
        for (int i=0; i<numElemGlueZoneFine*numberIntPoints; i++){
            output_v << aux + 1 << std::endl;
            aux += 1;
        };

        output_v << "      </DataArray>" << std::endl;
  
        //WRITE ELEMENT TYPES
        output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                 << "format=\"ascii\">" << std::endl;
    
        for (int i=0; i<numElemGlueZoneFine*numberIntPoints; i++){
            output_v << 1 << std::endl;
        };

        output_v << "      </DataArray>" << std::endl
                 << "    </Cells>" << std::endl;


        //FINALIZE OUTPUT FILE
            output_v << "  </Piece>" << std::endl;
            output_v << "  </UnstructuredGrid>" << std::endl
                     << "</VTKFile>" << std::endl;
        


    	//PRINT FINE MODEL RESULTS
        std::string f = "saidaFineIP"+result+".vtu";
        std::fstream output_vf(f.c_str(), std::ios_base::out);
    
        	
    	output_vf << "<?xml version=\"1.0\"?>" << std::endl
                 << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                 << "  <UnstructuredGrid>" << std::endl
                 << "  <Piece NumberOfPoints=\"" << numElemGlueZoneFine*numberIntPoints
                 << "\"  NumberOfCells=\"" << numElemGlueZoneFine*numberIntPoints
                 << "\">" << std::endl;

        //WRITE NODAL COORDINATES
        output_vf << "    <Points>" << std::endl
                 << "      <DataArray type=\"Float64\" "
                 << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;


		for (int i = 0; i< numElemGlueZoneFine; i++){

        	SpecialQuadrature squad;
            double wpc[9],x1[9], x2[9];

			int *connec = elementsFine_[elementsGlueZoneFine_[i]] -> getConnectivity();
            int patch = elementsFine_[elementsGlueZoneFine_[i]] -> getPatch();
            int *inc = nodesFine_[connec[8]] -> getINC();
                	
            for (int j = 0; j < 9; j++){
                double *x = nodesFine_[connec[j]] -> getCoordinates();
                wpc[j] = nodesFine_[connec[j]] -> getWeightPC();
                x1[j] = x[0];
                x2[j] = x[1];
            }; 

            int numberIntPoints = elementsFine_[elementsGlueZoneFine_[i]] -> 
                                  getNumberOfIntegrationPointsSpecial();

            for (int ip = 0; ip < numberIntPoints; ip++){

            	double x_[2];
                x_[0] = squad.interpolateQuadraticVariable(x1,ip,wpc,inc,IsoParFine,patch);
                x_[1] = squad.interpolateQuadraticVariable(x2,ip,wpc,inc,IsoParFine,patch);
            
				output_vf << x_[0] << " " << x_[1] << " " << 0.0 << std::endl;
            }

		}

        output_vf << "      </DataArray>" << std::endl
                 << "    </Points>" << std::endl;
        
      
        //WRITE ELEMENT CONNECTIVITY
        output_vf << "    <Cells>" << std::endl
                << "      <DataArray type=\"Int32\" "
                << "Name=\"connectivity\" format=\"ascii\">" << std::endl;
        
        for (int numN = 0; numN < numElemGlueZoneFine*numberIntPoints; ++numN){
                output_vf << numN << std::endl;
        }

        output_vf << "      </DataArray>" << std::endl;
  
        //WRITE OFFSETS IN DATA ARRAY
        output_vf << "      <DataArray type=\"Int32\""
                << " Name=\"offsets\" format=\"ascii\">" << std::endl;
    
        aux = 0;
        for (int i=0; i<numElemGlueZoneFine*numberIntPoints; i++){
            output_vf << aux + 1 << std::endl;
            aux += 1;
        };

        output_vf << "      </DataArray>" << std::endl;
  
        //WRITE ELEMENT TYPES
        output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                 << "format=\"ascii\">" << std::endl;
    
        for (int i=0; i<numElemGlueZoneFine*numberIntPoints; i++){
            output_vf << 1 << std::endl;
        };

        output_vf << "      </DataArray>" << std::endl
                 << "    </Cells>" << std::endl;

        //FINALIZE OUTPUT FILE
            output_vf << "  </Piece>" << std::endl;
            output_vf << "  </UnstructuredGrid>" << std::endl
                     << "</VTKFile>" << std::endl;
    };
};


//------------------------------------------------------------------------------
//---------------------------SETS FLUID MODELS----------------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::setFluidModels(FluidMesh& coarse, FluidMesh& fine){

    coarseModel = coarse;
    fineModel = fine;

    //Gets Fine and Coarse models basic information from fluid data
    numElemCoarse = coarseModel.elements_.size();
    numElemFine   = fineModel.elements_.size();
    numNodesCoarse = coarseModel.nodes_.size();
    numNodesFine   = fineModel.nodes_.size();

    nodesCoarse_  = coarseModel.nodes_;
    nodesFine_    = fineModel.nodes_;

    elementsCoarse_ = coarseModel.elements_;
    elementsFine_   = fineModel.elements_;
 
    boundaryCoarse_ = coarseModel.boundary_;
    boundaryFine_   = fineModel.boundary_;

    numBoundElemFine = boundaryFine_.size();
    numBoundElemCoarse = boundaryCoarse_.size();

    domDecompCoarse = coarseModel.getDomainDecomposition();
    domDecompFine = fineModel.getDomainDecomposition();

    parametersFine = &fineModel.fluidParameters;
    parametersCoarse = &coarseModel.fluidParameters;

    IsoParFine = fineModel.IsoPar_;
   	IsoParCoarse = coarseModel.IsoPar_;
    
    numTimeSteps = fineModel.numTimeSteps;
    dTime = fineModel.dTime;
    integScheme = fineModel.integScheme;
    glueZoneThickness = fineModel.glueZoneThickness;
    arlequinEpsilon = fineModel.arlequinEpsilon;


    // //Define coarse elements type
    // elemType = elementsCoarse_[0] -> getElemType();
   
    //Define fine model as true
    for (int i=0; i < numElemFine; i++){
        elementsFine_[i] -> setModel(true); 
    };
    //Define coarse model elements as false
    for (int i=0; i < numElemCoarse; i++){
        elementsCoarse_[i] -> setModel(false);
    }; 

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);     

    //Boxes for the integration points search 
    setElementBoxes();
    
   	//Signaled distance from the a defined fine boundary mesh
    setSignaledDistance();

    //Sets the fine and coarse nodes/elements in the gluing zone
    //Defines the lagrange multipliers mesh
    setGluingZone();

    //Computes the Weight function for all nodes and integration points
    setWeightFunction(); 

     printResults(-1);
   
};


template<>
void Arlequin<2>::setDirichletConstrain(std::vector<int> &dofTemp) {
   
	int dim = 2;

    //Non coincidente number nodes or control points in fine and coarse mesh
    int NCNumberNodesF = fineModel.NCNumberNodes;
    int NCNumberNodesC;
    // if (elemType == 0){  //fem mesh
    //     NCNumberNodesC = numNodesCoarse;
    // } else { //IGA mesh
        NCNumberNodesC = coarseModel.NCNumberNodes;
    // };

    //Coarse mesh
    for (int ielem = 0; ielem < numElemCoarse; ielem++){
    	
        int *connec = elementsCoarse_[ielem] -> getConnectivity();
    	
        // if (elemType == 0){
        //     for (int i = 0; i < 6 ; i++){
        //     	//velocity constrain
        //     	for (int j = 0; j < dim; j++){
        //     		int constrain = nodesCoarse_[connec[i]] -> getConstrains(j);
        //     		if ((constrain == 1) || (constrain == 3)){
        //             	dofTemp.push_back(connec[i]*dim + j);
        //         	};
        //     	};
        //         // //cavity pressure constrain
        //         // if (connec[i] == 2){
        //         //  dofTemp.push_back(2*NCNumberNodesC + connec[i]);
        //         // }    
        //     };
        // } else { // IGA mesh
            for (int i = 0; i < 9 ; i++){
                int newconi = nodesCoarse_[connec[i]] -> getnewcon();
                // velocity constrain
                for (int j = 0; j < dim; j++){
                	int constrain = nodesCoarse_[connec[i]] -> getConstrains(j);
                	if ((constrain == 1) || (constrain == 3)){
                    	dofTemp.push_back(newconi*dim + j);
                	};
                };                
            };
        // }; end else   
    };

    //Fine mesh (IGA elements)
    for (int ielem = 0; ielem < numElemFine; ielem++){
        int *connec = elementsFine_[ielem] -> getConnectivity();
        for (int i = 0; i < 9 ; i++){
            int newconi = nodesFine_[connec[i]] -> getnewcon();
            // velocity constrain
            for (int j = 0; j < dim; j++){
            	int constrain = nodesFine_[connec[i]] -> getConstrains(j);
            	if ((constrain == 1) || (constrain == 3)){
                	dofTemp.push_back(newconi*dim + (dim+1)*NCNumberNodesC + j);
            	};
            };
            //cavity pressure constrain
            // if (connec[i] == 0){
            //  dofTemp.push_back(3*NCNumberNodesC + 2*NCNumberNodesF + newconi);
            // }  
        };
    };
};


template<>
int Arlequin<2>::solveArlequinProblem(int iterNumber, double tolerance) {

    Mat               A,F,A_tSUPG,A_tPSPG,A_tLSIC;
    Vec               b, u, All,b_tSUPG,b_tPSPG,b_tLSIC;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ij, Ione, iterations, *dof;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val,val1,val2,val3;
    PetscViewer       viewer;
   
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank); 

    int dim = 2;

    //Computes the Nodal correspondence between fine nodes/integration points and coarse elements
    setCorrespondenceFine();  

    printResultsIP(-1);


    //Set the degre of freedom index with boundary condition
    std::vector<int> dofTemp;
    setDirichletConstrain(dofTemp);
    PetscMalloc1(dofTemp.size(),&dof);
	for (size_t i = 0; i < dofTemp.size(); i++){
		dof[i] = dofTemp[i];
	};

    //Non coincidente number nodes or control points in fine and coarse mesh
    int NCNumberNodesF = fineModel.NCNumberNodes; //IGA mesh
    int NCNumberNodesC; //Coarse mesh (FEM mesh or IGA mesh)
    // if (elemType == 0){ //FEM coarse mesh
    // 	NCNumberNodesC = numNodesCoarse;	
    // } else { // IGA coarse mesh
    	NCNumberNodesC = coarseModel.NCNumberNodes;	
    // }

    double &alpha_f = parametersFine -> getAlphaF();
    double &alpha_m = parametersFine -> getAlphaM();
    double &gamma = parametersFine -> getGamma();

    int sysSize = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*numNodesGlueZoneFine;

    double sumtime = 0;
 
    for (iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){

        
        if (rank == 0) {std::cout << "------------------------- TIME STEP = "
                                  << iTimeStep << " -------------------------"
                                  << std::endl;}
        
        for (int i = 0; i < numNodesCoarse; i++){
            
            double accel[2], u[2];

            //Compute acceleration
            u[0] = nodesCoarse_[i] -> getVelocity(0);
            u[1] = nodesCoarse_[i] -> getVelocity(1);
            
            nodesCoarse_[i] -> setPreviousVelocity(u);

            accel[0] = nodesCoarse_[i] -> getAcceleration(0);
            accel[1] = nodesCoarse_[i] -> getAcceleration(1);

            nodesCoarse_[i] -> setPreviousAcceleration(accel);

            accel[0] *= (gamma - 1.) / gamma;
            accel[1] *= (gamma - 1.) / gamma;

            nodesCoarse_[i] -> setAcceleration(accel);

        };

        for (int i = 0; i < numNodesFine; i++){
            
            double accel[2], u[2];
           
            //Compute acceleration
            u[0] = nodesFine_[i] -> getVelocity(0);
            u[1] = nodesFine_[i] -> getVelocity(1);
            
            nodesFine_[i] -> setPreviousVelocity(u);

            accel[0] = nodesFine_[i] -> getAcceleration(0);
            accel[1] = nodesFine_[i] -> getAcceleration(1);

            nodesFine_[i] -> setPreviousAcceleration(accel);

            accel[0] *= (gamma - 1.) / gamma;
            accel[1] *= (gamma - 1.) / gamma;

            nodesFine_[i] -> setAcceleration(accel);

        };

        double duNorm=100.;
        
        for (int inewton = 0; inewton < iterNumber; inewton++){
            
            
            boost::posix_time::ptime t1 = boost::posix_time::microsec_clock::local_time();
            
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                sysSize,sysSize,1000000,NULL,1000000,NULL,&A);CHKERRQ(ierr);

            for (int i=0; i<sysSize; i++){
                double valu = 1.e-20;
                ierr = MatSetValues(A,1,&i,1,&i,&valu,ADD_VALUES);
            }
            
            // Divides the matrix between the processes
            ierr = MatGetOwnershipRange(A, &Istart, &Iend);CHKERRQ(ierr);      	

            //Create PETSc vectors
            ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
            ierr = VecSetSizes(b,PETSC_DECIDE,sysSize);CHKERRQ(ierr);
            ierr = VecSetFromOptions(b);CHKERRQ(ierr);
            ierr = VecDuplicate(b,&u);CHKERRQ(ierr);
            ierr = VecDuplicate(b,&All);CHKERRQ(ierr);

            
            //Coarse Mesh - IGA mesh
           	for (int jel = 0; jel < numElemCoarse; jel++){  

            	if (domDecompCoarse.first[jel] == rank) {  

            		int *connec = elementsCoarse_[jel] -> getConnectivity();

            		double **elemMatrix;
	                elemMatrix = new double*[27]();
	                for (int i = 0; i < 27; ++i)  elemMatrix[i] = new double[27]();
	                double elemVector[27] = {};

	            	elementsCoarse_[jel] -> getTransientNavierStokes_ISO(elemMatrix,elemVector);

	            	for (int i=0; i<9; i++){

                            int newconi = nodesCoarse_[connec[i]] -> getnewcon();

	                        for (int j=0; j<9; j++){
	                             
                                int newconj = nodesCoarse_[connec[j]] -> getnewcon();
                                int dof_i = 2 * newconi;
                                int dof_j = 2 * newconj;
                                ierr = MatSetValues(A, 1, &dof_i,1, &dof_j,
                                                    &elemMatrix[2*i  ][2*j  ],
                                                    ADD_VALUES);

                                dof_i = 2 * newconi + 1;
                                dof_j = 2 * newconj;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &elemMatrix[2*i+1][2*j  ],
                                                    ADD_VALUES);

                                dof_i = 2 * newconi;
                                dof_j = 2 * newconj + 1;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &elemMatrix[2*i  ][2*j+1],
                                                    ADD_VALUES);

                                dof_i = 2 * newconi + 1;
                                dof_j = 2 * newconj + 1;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &elemMatrix[2*i+1][2*j+1],
                                                    ADD_VALUES);

                                dof_i = 2 * newconi;
                                dof_j = 2 * NCNumberNodesC + newconj;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &elemMatrix[2*i  ][18+j],
                                                    ADD_VALUES);

                                dof_i = 2 * NCNumberNodesC + newconi;
                                dof_j = 2 * newconj;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &elemMatrix[18+i][2*j  ],
                                                    ADD_VALUES);

                                dof_i = 2 * newconi + 1;
                                dof_j = 2 * NCNumberNodesC + newconj;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &elemMatrix[2*i+1][18+j],
                                                    ADD_VALUES);

                                dof_i = 2 * NCNumberNodesC + newconi;
                                dof_j = 2 * newconj + 1;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &elemMatrix[18+i][2*j+1],
                                                    ADD_VALUES);

                                dof_i = 2 * NCNumberNodesC + newconi;
                                dof_j = 2 * NCNumberNodesC + newconj;
                                ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                    &elemMatrix[18+i][18+j],
                                                    ADD_VALUES);
	                        }; //loop j
	                        
                        	//Rhs vector
	                        int dof_i = 2 * newconi;
	                        ierr = VecSetValues(b, 1, &dof_i, &elemVector[2*i  ],
	                                            ADD_VALUES);
	                    
	                        dof_i = 2 * newconi + 1;
	                        ierr = VecSetValues(b, 1, &dof_i, &elemVector[2*i+1],
	                                            ADD_VALUES);

	                        dof_i = 2 * NCNumberNodesC + newconi;
	                        ierr = VecSetValues(b, 1, &dof_i, &elemVector[18+i],
	                                            ADD_VALUES);
	                     };// loop i
	            	
                     for (int i = 0; i < 27; ++i) delete [] elemMatrix[i];
        			 delete [] elemMatrix;

	            };
	        };


	        //Fine mesh 
	        for (int jel = 0; jel < numElemFine; jel++){   
                
                if (domDecompFine.first[jel] == rank) { 

                    int *connec = elementsFine_[jel] -> getConnectivity();
                   
                	double **elemMatrix;
            		elemMatrix = new double*[27]();
            		for (int i = 0; i < 27; ++i)  elemMatrix[i] = new double[27]();
            		double elemVector[27] = {};

                    elementsFine_[jel] -> getTransientNavierStokes_ISO(elemMatrix,elemVector);

                    //Disperse local contributions into the global matrix
                    //Matrix K and C
                    for (int i=0; i<9; i++){

                        int newconi = nodesFine_[connec[i]] -> getnewcon();
                        
                        for (int j=0; j<9; j++){
                                
                            int newconj = nodesFine_[connec[j]] -> getnewcon();
                            
                            int dof_i = 2 * newconi + 3 * NCNumberNodesC;
                            int dof_j = 2 * newconj + 3 * NCNumberNodesC;
                            ierr = MatSetValues(A, 1, &dof_i,1, &dof_j,
                                                &elemMatrix[2*i  ][2*j  ],
                                                ADD_VALUES);

                            dof_i = 2 * newconi + 1 + 3 * NCNumberNodesC;
                            dof_j = 2 * newconj + 3 * NCNumberNodesC;
                            ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                &elemMatrix[2*i+1][2*j  ],
                                                ADD_VALUES);

                            dof_i = 2 * newconi + 3 * NCNumberNodesC;
                            dof_j = 2 * newconj + 1 + 3 * NCNumberNodesC;
                            ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                &elemMatrix[2*i  ][2*j+1],
                                                ADD_VALUES);

                            dof_i = 2 * newconi + 1 + 3 * NCNumberNodesC;
                            dof_j = 2 * newconj + 1 + 3 * NCNumberNodesC;
                            ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                &elemMatrix[2*i+1][2*j+1],
                                                ADD_VALUES);
                    
                         //Matrix Q and Qt
                            dof_i = 2 * newconi + 3 * NCNumberNodesC;
                            dof_j = 2 * NCNumberNodesF + newconj + 3 * NCNumberNodesC;
                            ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                &elemMatrix[2*i  ][18+j],
                                                ADD_VALUES);

                            dof_i = 2 * NCNumberNodesF + newconi + 3 * NCNumberNodesC;
                            dof_j = 2 * newconj + 3 * NCNumberNodesC;
                            ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                &elemMatrix[18+i][2*j  ],
                                                ADD_VALUES);

                            dof_i = 2 * newconi + 1 + 3 * NCNumberNodesC;
                            dof_j = 2 * NCNumberNodesF + newconj + 3 * NCNumberNodesC;
                            ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                &elemMatrix[2*i+1][18+j],
                                                ADD_VALUES);

                            dof_i = 2 * NCNumberNodesF + newconi + 3 * NCNumberNodesC;
                            dof_j = 2 * newconj + 1 + 3 * NCNumberNodesC;
                            ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                &elemMatrix[18+i][2*j+1],
                                                ADD_VALUES);

                            dof_i = 2 * NCNumberNodesF + newconi + 3 * NCNumberNodesC;
                            dof_j = 2 * NCNumberNodesF + newconj + 3 * NCNumberNodesC;
                            ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                                &elemMatrix[18+i][18+j],
                                                ADD_VALUES);
                        }; //loop j
                        
                        //Rhs vector
                        int dof_i = 2 * newconi + 3 * NCNumberNodesC;
                        ierr = VecSetValues(b, 1, &dof_i, &elemVector[2*i  ],
                                            ADD_VALUES);

                        dof_i = 2 * newconi + 1 + 3 * NCNumberNodesC;
                        ierr = VecSetValues(b, 1, &dof_i, &elemVector[2*i+1],
                                            ADD_VALUES);

                        dof_i = 2 * NCNumberNodesF + newconi + 3 * NCNumberNodesC;
                        ierr = VecSetValues(b, 1, &dof_i, &elemVector[18+i],
                                            ADD_VALUES);

                     };// loop i
                     for (int i = 0; i < 27; ++i) delete [] elemMatrix[i];
        			 delete [] elemMatrix;  

	            }; // domain decomposition
            }; //Elements Fine

            //numElemGlueZoneFine
            for (int l = 0; l < numElemGlueZoneFine; l++){

            	int jel =  elementsGlueZoneFine_[l];

            	if (domDecompFine.first[jel] == rank) {

	            	//FINE MESH
	            	int *connec = elementsFine_[jel] -> getConnectivity();
	            	int *connecL = glueZoneFine_[l] -> getConnectivity();	

	            	//LAGRANGE MULTIPLIERS MATRIXES AND VECTORS
	        		double **elemMatrixLag1;
	        		elemMatrixLag1 = new double*[18]();
	        		for (int i = 0; i < 18; ++i)  elemMatrixLag1[i] = new double[18]();
	        		double elemVectorLag1[18] = {};
	        		double elemVectorLag2[18] = {};
	        		
	        		elementsFine_[jel] -> getLagrangeMultipliersSameMesh(elemMatrixLag1,elemVectorLag1,elemVectorLag2);

	        		//ARLEQUIN STABILIZATION MATRIX
	        		double **elemStabMatrix;
	        		elemStabMatrix = new double*[18]();
	        		for (int i = 0; i < 18; ++i)  elemStabMatrix[i] = new double[18]();
	        		double elemStabVector[18] = {};

	        		elementsFine_[jel] -> getLagrangeMultipliersArlequinSameMesh(elemStabMatrix, elemStabVector);


           //          if (inewton == 1) {

           //              for (int i = 0; i < 18; i++) std::cout << elemStabVector [i] << " " << elemVectorLag2[i] << std::endl;
           //          }
        		
        			double integ = alpha_f * gamma * dTime;
	        		for (int i = 0; i < 9; i++){
	        			for (int j = 0; j < 9; j++){

	        				int newconj = nodesFine_[connec[j]] -> getnewcon();

	        				int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	        				int dof_j = 3*NCNumberNodesC + 2*newconj;
	        				double value = integ * elemMatrixLag1[2*i][2*j];
	        				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                            &value,ADD_VALUES);
	        				ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
	                                            &elemMatrixLag1[2*i][2*j],
	                                            ADD_VALUES);

	        				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	        				dof_j = 3*NCNumberNodesC + 2*newconj;
	        				value = integ * elemMatrixLag1[2*i+1][2*j];
	        				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                            &value,ADD_VALUES);
	        				ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
	                                            &elemMatrixLag1[2*i+1][2*j],
	                                            ADD_VALUES);

	        				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	        				dof_j = 3*NCNumberNodesC + 2*newconj + 1;
	        				value = integ * elemMatrixLag1[2*i][2*j+1];
	        				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                            &value,ADD_VALUES);
	        				ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
	                                            &elemMatrixLag1[2*i][2*j+1],
	                                            ADD_VALUES);

	        				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	        				dof_j = 3*NCNumberNodesC + 2*newconj + 1;
	        				value = integ * elemMatrixLag1[2*i+1][2*j+1];
	        				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                            &value,ADD_VALUES);
	        				ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
	                                            &elemMatrixLag1[2*i+1][2*j+1],
	                                            ADD_VALUES);

	        				//Stabilization Arlequin Terms diagonal
	        				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	        				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	        				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                            &elemStabMatrix[2*i][2*j],
	                                            ADD_VALUES);

	        				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	        				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	        				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                            &elemStabMatrix[2*i+1][2*j],
	                                            ADD_VALUES);

	        				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	        				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	        				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                            &elemStabMatrix[2*i][2*j+1],
	                                            ADD_VALUES);

	        				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	        				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	        				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                            &elemStabMatrix[2*i+1][2*j+1],
	                                            ADD_VALUES);

	        			};//j

	        			int newconi = nodesFine_[connec[i]] -> getnewcon();

	        			int dof_i = 3*NCNumberNodesC + 2*newconi;
	        			ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1[2*i],ADD_VALUES);
	        			dof_i = 3*NCNumberNodesC + 2*newconi + 1;
	        			ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1[2*i+1],ADD_VALUES);


	        			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	        			ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag2[2*i  ],ADD_VALUES);
	        			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	        			ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag2[2*i+1 ],ADD_VALUES);

	        			//Stabilization Arlequin Term
	        			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	        			ierr = VecSetValues(b, 1, &dof_i, &elemStabVector[2*i  ],ADD_VALUES);
	        			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	        			ierr = VecSetValues(b, 1, &dof_i, &elemStabVector[2*i+1 ],ADD_VALUES);

	        		};//i


	        		int numberIntPoints = elementsFine_[jel] -> getNumberOfIntegrationPointsSpecial();

            	        		
	        		std::vector<int> ele, diffElem;
	            	ele.clear();
	            	diffElem.clear();

	            	//Finding coarse elements 
	            	for (int i=0; i<numberIntPoints; i++){
	                	int aux = elementsFine_[jel] -> getIntegPointCorrespondenceElement(i);
	                	ele.push_back(aux);
	            	};

	            	int numElemIntersect = 1;
	            	int flag = 0;
	            	diffElem.push_back(ele[0]);
	            
		            for (int i = 1; i<numberIntPoints; i++){
		                flag = 0;
		                for (int j = 0; j<numElemIntersect; j++){
		                    if (ele[i] == diffElem[j]) {
		                        break;
		                    }else{
		                        flag++;
		                    };
		                    if(flag == numElemIntersect){
		                        numElemIntersect++;
		                        diffElem.push_back(ele[i]);
		                    };
		                };
		            };


		            for (int ielem = 0; ielem < numElemIntersect; ielem++){

		            	for (int i = 0; i < 18; i++){
	        				elemVectorLag1[i] = 0.0;
	        				elemVectorLag2[i] = 0.0;
	        				elemStabVector[i] = 0.0;
		        			for (int j = 0; j < 18; j++){
		        				elemMatrixLag1[i][j] = 0.0;
		        				elemStabMatrix[i][j] = 0.0;
		        			}
	        			}

	                
	                	int iElemCoarse = diffElem[ielem];

	                	int *connecC = elementsCoarse_[iElemCoarse] -> getConnectivity();
	                	int patch = elementsCoarse_[iElemCoarse] -> getPatch();
	                	
	                	elementsFine_[jel] -> getLagrangeMultipliersDifferentMesh_ISO(patch,nodesCoarse_,connecC,IsoParCoarse,iElemCoarse, 
													 								  elemMatrixLag1,elemVectorLag1,elemVectorLag2);  

	                	//elementsFine_[jel] -> getLagrangeMultipliersArlequinDifferentMesh_ISO(patch,nodesCoarse_,connecC,IsoParCoarse,iElemCoarse,
                  	                                										  // elemStabMatrix,elemStabVector);



	                	for (int i = 0; i < 9; i++){
	        				
	        				for (int j = 0; j < 9; j++){

	        				int newconj = nodesCoarse_[connecC[j]] -> getnewcon();	
	        				
	        				int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	        				int dof_j = 2*newconj;
	        				double value = integ * elemMatrixLag1[2*i][2*j];
	        				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                            &value,ADD_VALUES);
	        				ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
	                                            &elemMatrixLag1[2*i][2*j],
	                                            ADD_VALUES);


	        				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i]+1;
	        				dof_j = 2*newconj;
	        				value = integ * elemMatrixLag1[2*i+1][2*j];
	        				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                            &value,ADD_VALUES);
	        				ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
	                                            &elemMatrixLag1[2*i+1][2*j],
	                                            ADD_VALUES);

	        				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	        				dof_j = 2*newconj +1 ;
	        				value = integ * elemMatrixLag1[2*i][2*j+1];
	        				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                            &value,ADD_VALUES);
	        				ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
	                                            &elemMatrixLag1[2*i][2*j+1],
	                                            ADD_VALUES);

	        				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1 ;
	        				dof_j = 2*newconj +1 ;
	        				value = integ * elemMatrixLag1[2*i+1][2*j+1];
	        				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                            &value,ADD_VALUES);
	        				ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
	                                            &elemMatrixLag1[2*i+1][2*j+1],
	                                            ADD_VALUES);

	        				//Arlequin Stabilization Terms
	        				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	        				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	        				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                            &elemStabMatrix[2*i][2*j],
	                                            ADD_VALUES);

	        				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	        				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	        				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                            &elemStabMatrix[2*i+1][2*j],
	                                            ADD_VALUES);

	        				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	        				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	        				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                            &elemStabMatrix[2*i][2*j+1],
	                                            ADD_VALUES);

	        				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	        				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	        				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                            &elemStabMatrix[2*i+1][2*j+1],
	                                            ADD_VALUES);

	        				};//i

	        				int newconi = nodesCoarse_[connecC[i]] -> getnewcon();

	        				int dof_i = 2*newconi;
	        				ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1[2*i  ],ADD_VALUES);
	        				dof_i = 2*newconi + 1 ;
	        				ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1[2*i+1],ADD_VALUES);

	        				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	        				ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag2[2*i  ],ADD_VALUES);

	        				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	        				ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag2[2*i+1],ADD_VALUES);

	        				//Stabilization Arlequin Term
		        			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
		        			ierr = VecSetValues(b, 1, &dof_i, &elemStabVector[2*i  ],ADD_VALUES);
		        			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
		        			ierr = VecSetValues(b, 1, &dof_i, &elemStabVector[2*i+1 ],ADD_VALUES);
	        				
	        			};//j

	                };//intersect
               

	                for (int i = 0; i < 18; ++i) {
	                	delete [] elemMatrixLag1[i];
	                	delete [] elemStabMatrix[i];
	                }
	        		delete [] elemMatrixLag1; 
	        		delete [] elemStabMatrix;
	        	};//decomposition

            };//gluezonefine


            //Assemble matrices and vectors
            ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
            
            ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
            ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

            
            //MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            //VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

            MatZeroRowsColumns(A,dofTemp.size(),dof,1.0,u,b);
            
            //Create KSP context to solve the linear system
            ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
            ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
            
            //Solve using GMRES
            // ierr = KSPSetTolerances(ksp,1.e-10,PETSC_DEFAULT,PETSC_DEFAULT,500); CHKERRQ(ierr);
            // ierr = KSPGetPC(ksp,&pc);
            // ierr = PCSetType(pc,PCNONE);
            // ierr = KSPSetType(ksp,KSPDGMRES); CHKERRQ(ierr);
            // ierr = KSPGMRESSetRestart(ksp, 500); CHKERRQ(ierr);

            // //Solve using Mumps
            #if defined(PETSC_HAVE_MUMPS)
            ierr = KSPSetType(ksp,KSPPREONLY);
            ierr = KSPGetPC(ksp,&pc);
            ierr = PCSetType(pc, PCLU);
            #endif          
            ierr = KSPSetFromOptions(ksp);CHKERRQ(ierr);
            ierr = KSPSetUp(ksp);


            //Solve Linear System
            ierr = KSPSolve(ksp,b,u);CHKERRQ(ierr);
            ierr = KSPGetTotalIterations(ksp, &iterations);            

             //MatView(A,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            //VecView(b,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
            //VecView(u,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

            //Gathers the solution vector to the master process
            ierr = VecScatterCreateToAll(u, &ctx, &All);CHKERRQ(ierr);
            ierr = VecScatterBegin(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
            ierr = VecScatterEnd(ctx, u, All, INSERT_VALUES, SCATTER_FORWARD);CHKERRQ(ierr);
            ierr = VecScatterDestroy(&ctx);CHKERRQ(ierr);
            
            //Updates nodal values
            
            double u_;
            double p_;
            double lag_;
            double normU = 0.;
            double normP = 0.;
            double normL = 0.;
            Ione = 1;
        
            for (int i = 0; i < numNodesCoarse; ++i){
        
             //   	if (elemType == 0){
	               
	            //     Ii = 2*i;
	            //     ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
	            //     nodesCoarse_[i] -> incrementAcceleration(0,val);
	            //     nodesCoarse_[i] -> incrementVelocity(0,val*gamma*dTime);
	            //     duNorm += val*val;
	            
	            //     Ii = 2*i+1;
	            //     ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
	            //     nodesCoarse_[i] -> incrementAcceleration(1,val);
	            //     nodesCoarse_[i] -> incrementVelocity(1,val*gamma*dTime);
	            //     duNorm += val*val;

	            //     pressure
	            
	            // } else {

	            	int newconi = nodesCoarse_[i] -> getnewcon();

	            	Ii = 2*newconi;
	                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
	                u_ = val;
	                nodesCoarse_[i] -> incrementAcceleration(0,u_);
	                nodesCoarse_[i] -> incrementVelocity(0,u_*gamma*dTime);
	                normU += val*val;
	            
	                Ii = 2*newconi+1;
	                ierr = VecGetValues(All, Ione, &Ii, &val);CHKERRQ(ierr);
	                u_ = val;
	                nodesCoarse_[i] -> incrementAcceleration(1,u_);
	                nodesCoarse_[i] -> incrementVelocity(1,u_*gamma*dTime);
	                normU += val*val;

	                Ii = 2*NCNumberNodesC + newconi;
	                ierr = VecGetValues(All,Ione,&Ii,&val);CHKERRQ(ierr);
	                p_ = val;
	                normP += val*val;
	                nodesCoarse_[i] -> incrementPressure(p_);
	                
            };	
             


            for (int i = 0; i < numNodesFine; ++i){

                int newconi = nodesFine_[i] -> getnewcon();
        
                Ii = 2*newconi + 3*NCNumberNodesC;
                ierr = VecGetValues(All, Ione, &Ii, &val);
                u_ = val;
                nodesFine_[i] -> incrementAcceleration(0,u_);
                nodesFine_[i] -> incrementVelocity(0,u_*gamma*dTime);
                normU += val*val;
            
                Ii = 2*newconi+1 + 3*NCNumberNodesC;
                ierr = VecGetValues(All, Ione, &Ii, &val);
                u_ = val;
                nodesFine_[i] -> incrementAcceleration(1,u_);
                nodesFine_[i] -> incrementVelocity(1,u_*gamma*dTime);
                normU += val*val;

                Ii = 2*NCNumberNodesF + newconi+ 3*NCNumberNodesC;
                ierr = VecGetValues(All,Ione,&Ii,&val);
                p_ = val;
                nodesFine_[i] -> incrementPressure(p_);
                normP += val*val;
            };
             
            for (int i = 0; i < numNodesGlueZoneFine; ++i){

            	Ii = 3*NCNumberNodesF + 3*NCNumberNodesC + 2*i;
            	ierr = VecGetValues(All, Ione, &Ii, &val);
            	lag_ = val;
            	normL += val * val;
            	nodesFine_[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(0,lag_);

            	Ii = 3*NCNumberNodesF + 3*NCNumberNodesC + 2*i + 1;
            	ierr = VecGetValues(All, Ione, &Ii, &val);
            	lag_ = val;
            	normL += val * val;
            	nodesFine_[nodesGlueZoneFine_[i]] -> incrementLagrangeMultiplier(1,lag_);


            };
			           
   
            boost::posix_time::ptime t2 = boost::posix_time::microsec_clock::local_time();

            boost::posix_time::time_duration diff = t2 - t1;

            sumtime += diff.total_milliseconds()/1000.;

         
            if(rank == 0){
                                
                std::cout << "Iteration = " << inewton 
                          << " (" << iterations << ")"  
                          << "   Du Norm = " << std::scientific << sqrt(normU) 
                          << " " << sqrt(normP)  << " " << sqrt(normL)
                          << "  Time (s) = " << std::fixed
                          << diff.total_milliseconds()/1000. << std::endl;
            };
                      
            ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
            ierr = VecDestroy(&b); CHKERRQ(ierr);
            ierr = VecDestroy(&u); CHKERRQ(ierr);
            ierr = VecDestroy(&All); CHKERRQ(ierr);
            ierr = MatDestroy(&A); CHKERRQ(ierr);

            if (sqrt(normU) <= tolerance) {
                break;
            };

            
        };//Newton-Raphson


        if(rank == 0){           
            std::cout << "Accumulated time = " << sumtime << std::endl;
        };


        //Computes the real velocity
        QuadShapeFunction<2> shapeQuad;
                
        for (int i = 0; i<numNodesFine; i++){
            for (int k = 0; k < dim; k++) nodesFine_[i] -> setVelocityArlequin(k,nodesFine_[i] -> getVelocity(k));
            nodesFine_[i] -> setPressureArlequin(nodesFine_[i] ->getPressure());
        };

        for (int i=0; i<numNodesCoarse; i++){
        	for (int k = 0; k < dim; k++) nodesCoarse_[i] -> setVelocityArlequin(k,nodesCoarse_[i] -> getVelocity(k));
            nodesCoarse_[i] -> setPressureArlequin(nodesCoarse_[i] -> getPressure());
        };
        
        for (int i = 0; i<numNodesGlueZoneFine; i++){
            
            int elCoarse = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalElemCorrespondence();
            double* xsiC = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalXsiCorrespondence();


            int *connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();
           
            double u_coarse[9], v_coarse[9], p_coarse[9];
            for (int j=0; j<9; j++){
                u_coarse[j] = nodesCoarse_[connecCoarse[j]] -> getVelocity(0);
                v_coarse[j] = nodesCoarse_[connecCoarse[j]] -> getVelocity(1);
                p_coarse[j] = nodesCoarse_[connecCoarse[j]] -> getPressure();
            };

            int patchC = elementsCoarse_[elCoarse] -> getPatch();
    		int *incC = nodesCoarse_[connecCoarse[8]] -> getINC();
    		double wpcC[9],phiC_[9];
    		for (int k = 0; k<9; k++) wpcC[k] = nodesCoarse_[connecCoarse[k]] -> getWeightPC();
            shapeQuad.evaluateIso(xsiC,phiC_,wpcC,incC,IsoParCoarse,patchC);
            
            double u = 0.;
            double v = 0.;
            double p = 0.;
            for (int j=0; j<9; j++){
                u += u_coarse[j] * phiC_[j];
                v += v_coarse[j] * phiC_[j];
                p += p_coarse[j] * phiC_[j];
            };
            
            double wFunc = nodesFine_[nodesGlueZoneFine_[i]] -> getWeightFunction();
            
            double u_int = nodesFine_[nodesGlueZoneFine_[i]] -> getVelocity(0) * wFunc + u * (1. - wFunc);
            double v_int = nodesFine_[nodesGlueZoneFine_[i]] -> getVelocity(1) * wFunc + v * (1. - wFunc);
            double p_int = nodesFine_[nodesGlueZoneFine_[i]] -> getPressure() * wFunc + p * (1. - wFunc);
            
            nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(0,u_int);
            nodesFine_[nodesGlueZoneFine_[i]] -> setVelocityArlequin(1,v_int);
            nodesFine_[nodesGlueZoneFine_[i]] -> setPressureArlequin(p_int);
            
        };
       
        
        
        // Compute and print drag and lift coefficients
        // if (computeDragAndLift){
        //     dragAndLiftCoefficients(dragLift);
        // };

        // Printing results
        printResults(iTimeStep);
        
    };

    PetscFree(dof);
    return 0;
    
};





#endif
