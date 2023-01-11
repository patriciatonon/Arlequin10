
#include "Arlequin.h"

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//-----------------------------SET ELEMENT BOXES--------------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::setElementBoxes() {
   
	int dim = 2;

    //Compute element boxes for coarse model (FEM/IGA coarse mesh)
    
    double &alpha_f = parametersCoarse -> getAlphaF();

    if (elemTypeCoarse == 0) { //FEM mesh

        for (int jel = 0; jel < numElemCoarse; jel++){
            
            int *connec = elementsCoarse_[jel] -> getConnectivity();
            double xk[dim], Xk[dim];

            double *xx1,*xxp1,*xx2,*xxp2,*xx3,*xxp3;
            double x1[2],x2[2],x3[2];

            xx1 = nodesCoarse_[connec[0]] -> getCoordinates();
            xx2 = nodesCoarse_[connec[1]] -> getCoordinates();
            xx3 = nodesCoarse_[connec[2]] -> getCoordinates();
            // xxp1 = nodesCoarse_[connec[0]] -> getPreviousCoordinates();
            // xxp2 = nodesCoarse_[connec[1]] -> getPreviousCoordinates();
            // xxp3 = nodesCoarse_[connec[2]] -> getPreviousCoordinates(); 

            // x1[0] = alpha_f * xx1[0] + (1. - alpha_f) * xxp1[0];
            // x1[1] = alpha_f * xx1[1] + (1. - alpha_f) * xxp1[1];
            // x2[0] = alpha_f * xx2[0] + (1. - alpha_f) * xxp2[0];
            // x2[1] = alpha_f * xx2[1] + (1. - alpha_f) * xxp2[1];
            // x3[0] = alpha_f * xx3[0] + (1. - alpha_f) * xxp3[0];
            // x3[1] = alpha_f * xx3[1] + (1. - alpha_f) * xxp3[1];

            x1[0] = xx1[0];
            x1[1] = xx1[1];
            x2[0] = xx2[0];
            x2[1] = xx2[1];
            x3[0] = xx3[0];
            x3[1] = xx3[1];
     
            xk[0] = std::min(x1[0],std::min(x2[0], x3[0]));
            xk[1] = std::min(x1[1],std::min(x2[1], x3[1]));

            Xk[0] = std::max(x1[0],std::max(x2[0], x3[0]));
            Xk[1] = std::max(x1[1],std::max(x2[1], x3[1]));        
        
            elementsCoarse_[jel] -> setIntersectionParameters(xk, Xk);

        };

    } else { //IGA mesh

        for (int jel = 0; jel < numElemCoarse; jel++){
            
            int *connec = elementsCoarse_[jel] -> getConnectivity();
            double xk[dim], Xk[dim];

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

                // double *xx,*xxp;
                // double *xx = nodesCoarse_[connec[i]]->getCoordinates();
                // double *xxp = nodesCoarse_[connec[i]]->getPreviousCoordinates();
                //for (int j = 0; j < dim; j++) coord_[i][j] = alpha_f * xx[j] + (1. - alpha_f) * xxp[j];   
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

        };

    };

        

    return;
};

//------------------------------------------------------------------------------
//--------------------------SEARCH CORRESPONDENCE-------------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::searchPointCorrespondence(double *x,std::vector<Nodes *> nodes, 
                                            std::vector<Element *> elements, 
                                            std::vector<IsoParameters* > isopar, int &elemType,
                                            int numElem, double *xsiC, int &elemC, int elSearch){

    int dim = 2;
    QuadShapeFunction<2> shapeQuad;
    double x_[dim],deltaX[dim],deltaXsi[dim],xsi[dim],xsiCC[dim+1];

    double **ainv;
    ainv = new double*[dim];
    for (int i = 0; i < dim; ++i) ainv[i] = new double[dim];

    // double &alpha_f = isopar -> getAlphaF();

    if (elemType == 0) { //FEM mesh

        for (int i = 0; i<dim; i++){
           xsi[i] = 1./3.; //central element cooordinates
           x_[i] = 0.0;
           xsiC[i] = 1.e50;
        };
        for (int i = 0; i<dim+1; i++) xsiCC[i] = 1.e10;
      
        int *connec = elements[elSearch] -> getConnectivity();
                
        //Computing basis functions 
        double   phi_[6];
        shapeQuad.evaluateFem(xsi,phi_);

        for (int i = 0; i < 6; i++){
            double *xint = nodes[connec[i]] -> getCoordinates();
            for (int j = 0; j < dim; j++) x_[j] += xint[j] * phi_[i];  

            // double *xint = nodes[connec[i]] -> getCoordinates();
            // double *xintp = nodes[connec[i]] -> getPreviousCoordinates();
            // for (int j = 0; j < dim; j++) x_[j] += (alpha_f * xint[j] + (1. - alpha_f) * xintp[j]) * phi_[i];

        };

        double error = 1.e6;
        int iterations = 0;

        while ((error > 1.e-8) && (iterations < 4)) {
            
            iterations++;

            for (int i = 0; i < dim; i++){
                deltaX[i] = x[i] - x_[i]; 
                deltaXsi[i] = 0.0;
            };

            elements[elSearch] -> getJacobianMatrixValues_FEM(xsi,ainv);

            double tempAinv[dim][dim];
            for (int i = 0; i <dim; i++){
                for (int j =0; j<dim; j++){
                    tempAinv[i][j] = ainv[j][i];
                }
            }

            for (int i = 0; i <dim; i++){
                for (int j =0; j<dim; j++){
                    ainv[i][j] = tempAinv[i][j];
                }
            }

            for (int i = 0; i <dim; i++){
                for (int j = 0; j<dim; j++){
                    deltaXsi[i] += ainv[i][j] * deltaX[j];
                };
            };
            
            for (int i = 0; i < dim; i++){
               xsi[i] += deltaXsi[i];
               x_[i] = 0.0;
            };
                               
            shapeQuad.evaluateFem(xsi,phi_);
            
            for (int i = 0; i < 6; i++){
                double *xint = nodes[connec[i]] -> getCoordinates();
                for (int j = 0; j < dim; j++) x_[j] += xint[j]* phi_[i];   

                // double *xint = nodes[connec[i]] -> getCoordinates();
                // double *xintp = nodes[connec[i]] -> getPreviousCoordinates();
                // for (int j = 0; j < dim; j++) x_[j] += (alpha_f * xint[j] + (1. - alpha_f) * xintp[j]) * phi_[i];
               
            };
                              
            error = sqrt(deltaXsi[0]*deltaXsi[0] + deltaXsi[1]*deltaXsi[1]);
        };
        
        double t1 = -1.e-2;
        double t2 =  1. - t1;
        
        xsiCC[0] = xsi[0];
        xsiCC[1] = xsi[1];       
        xsiCC[2] = 1. - xsiCC[0] - xsiCC[1];

        if ((xsiCC[0] >= t1) && (xsiCC[1] >= t1) && (xsiCC[2] >= t1) &&
            (xsiCC[0] <= t2) && (xsiCC[1] <= t2) && (xsiCC[2] <= t2)){

            xsiC[0] = xsi[0];
            xsiC[1] = xsi[1];
            elemC = elSearch;
            

        } else {

            for (int jel = 0; jel < numElem; jel++){

                int *connec = elements[jel] -> getConnectivity();
              
                //get boxes information  
                std::pair<double*,double*> XK;      
                XK = elements[jel] -> getXIntersectionParameter();

                //Chech if the node is inside the element box
                if ((x[0] < XK.first[0]) || (x[0] > XK.second[0]) ||
                    (x[1] < XK.first[1]) || (x[1] > XK.second[1])) continue;
                                            
               
                for (int i = 0; i < dim; i++) {
                    xsi[i] = 1./3.; //central element cooordinates
                    x_[i] = 0.0;
                };

                double phi_[6];
                shapeQuad.evaluateFem(xsi,phi_);
               
                for (int i = 0; i < 6; i++){
                    double *xint = nodes[connec[i]] -> getCoordinates();
                    for (int j = 0; j < dim; j++) x_[j] += xint[j] * phi_[i]; 

                    // double *xint = nodes[connec[i]] -> getCoordinates();
                    // double *xintp = nodes[connec[i]] -> getPreviousCoordinates();
                    // for (int j = 0; j < dim; j++) x_[j] += (alpha_f * xint[j] + (1. - alpha_f) * xintp[j]) * phi_[i];                 
                };

                double error = 1.e6;
                
                int iterations = 0;

                while ((error > 1.e-8) && (iterations < 4)) {
                    
                    iterations++;
                    
                    for (int i = 0; i < dim; i++){
                        deltaX[i] = x[i] - x_[i]; 
                        deltaXsi[i] = 0.0;
                    };

                    elements[jel] -> getJacobianMatrixValues_FEM(xsi,ainv);

                    double tempAinv[dim][dim];
                    for (int i = 0; i <dim; i++){
                        for (int j =0; j<dim; j++){
                            tempAinv[i][j] = ainv[j][i];
                        }
                    }

                    for (int i = 0; i <dim; i++){
                        for (int j =0; j<dim; j++){
                            ainv[i][j] = tempAinv[i][j];
                        }
                    }
                    
                    for (int i = 0; i <dim; i++){
                        for (int j = 0; j<dim; j++){
                            deltaXsi[i] += ainv[i][j] * deltaX[j];
                        };
                    };
                    
                    for (int i = 0; i < dim; i++){
                        xsi[i] += deltaXsi[i];
                        x_[i] = 0.0;
                    };
                               
                    shapeQuad.evaluateFem(xsi,phi_);
                    
                    for (int i = 0; i < 6; i++){
                        double *xint = nodes[connec[i]] -> getCoordinates();
                        for (int j = 0; j < dim; j++) x_[j] += xint[j] * phi_[i];    

                        // double *xint = nodes[connec[i]] -> getCoordinates();
                        // double *xintp = nodes[connec[i]] -> getPreviousCoordinates();
                        // for (int j = 0; j < dim; j++) x_[j] += (alpha_f * xint[j] + (1. - alpha_f) * xintp[j]) * phi_[i];                
                    };
                                      
                    error = sqrt(deltaXsi[0]*deltaXsi[0] + deltaXsi[1]*deltaXsi[1]);
                };
                
                double t1 = -1.e-2;
                double t2 =  1. - t1;
                
                xsiCC[0] = xsi[0];
                xsiCC[1] = xsi[1];       
                xsiCC[2] = 1. - xsiCC[0] - xsiCC[1];

                if ((xsiCC[0] >= t1) && (xsiCC[1] >= t1) && (xsiCC[2] >= t1) &&
                    (xsiCC[0] <= t2) && (xsiCC[1] <= t2) && (xsiCC[2] <= t2)){

                    xsiC[0] = xsi[0];
                    xsiC[1] = xsi[1];
                    elemC = jel;
                };            
            }; //loop elements
        }; //else loop elements

    } else { //IGA mesh

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
        shapeQuad.evaluateIso(xsi,phi_,wpc,INC_,isopar,patch);

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
                               
            shapeQuad.evaluateIso(xsi,phi_,wpc,INC_,isopar,patch);
            
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

                shapeQuad.evaluateIso(xsi,phi_,wpc,INC_,isopar,patch);
               
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
                               
                    shapeQuad.evaluateIso(xsi,phi_,wpc,INC_,isopar,patch);
                    
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
            }; 

        }; 
    }; //elemType
    
    
    if (fabs(xsi[0]) > 2.) std::cout << "PROBEM SEARCHING NODE CORRESPONDENCE " 
                                         << std::endl;  


    for (int i = 0; i < 2; ++i) delete [] ainv[i];
    delete [] ainv;

};

template<>
void Arlequin<2>::setCorrespondenceFine() {

    int dim = 2;

    //Node correspondence
    for (int inode = 0; inode < numNodesGlueZoneFine; inode++) {

        double* x = nodesFine_[nodesGlueZoneFine_[inode]] -> getCoordinates();

        int elemC = 0;
        double xsiC[dim] = {};

        searchPointCorrespondence(x, nodesCoarse_, elementsCoarse_, IsoParCoarse, elemTypeCoarse,
        						  elementsCoarse_.size(),xsiC,elemC,
                                  nodesFine_[nodesGlueZoneFine_[inode]] -> getNodalElemCorrespondence());
      
        nodesFine_[nodesGlueZoneFine_[inode]] -> setNodalCorrespondence(elemC,xsiC);             

    };
    
    //integration points correspondence
    for (int i = 0; i< numElemGlueZoneFine; i++){

        SpecialQuadrature squad;
        
        int *connec = elementsFine_[elementsGlueZoneFine_[i]] -> getConnectivity();

        if (elemTypeFine == 0) { //FEM mesh

	        double x1[6], x2[6];
	        for (int j = 0; j < 6; j++){
	            double *x = nodesFine_[connec[j]] -> getCoordinates();
	            x1[j] = x[0];
	            x2[j] = x[1];
	        };   

	            
	        int numberIntPoints = elementsFine_[elementsGlueZoneFine_[i]] -> 
	                              getNumberOfIntegrationPointsSpecial_FEM();

	        for (int ip = 0; ip < numberIntPoints; ip++){

	            int elemC = 0;
	            double xsiC[dim] = {};
	            double x_[dim] = {};

	            x_[0] = squad.interpolateQuadraticVariableFem(x1,ip);
	            x_[1] = squad.interpolateQuadraticVariableFem(x2,ip);

	            searchPointCorrespondence(x_,nodesCoarse_,elementsCoarse_,IsoParCoarse, elemTypeCoarse,
	                                      elementsCoarse_.size(),xsiC,elemC, 
	                                      elementsFine_[elementsGlueZoneFine_[i]] -> getIntegPointCorrespondenceElement_FEM(ip));

            
	            elementsFine_[elementsGlueZoneFine_[i]] -> setIntegrationPointCorrespondence_FEM(ip,xsiC,elemC);

	        };  

        } else { //IGA mesh

	        int patch = elementsFine_[elementsGlueZoneFine_[i]] -> getPatch();
	        int *inc = nodesFine_[connec[8]] -> getINC();

	        double x1[9], x2[9], wpc[9];
	        for (int j = 0; j < 9; j++){
	            double *x = nodesFine_[connec[j]] -> getCoordinates();
	            wpc[j] = nodesFine_[connec[j]] -> getWeightPC();
	            x1[j] = x[0];
	            x2[j] = x[1];
	        };   

	            
	        int numberIntPoints = elementsFine_[elementsGlueZoneFine_[i]] -> 
	                              getNumberOfIntegrationPointsSpecial_ISO();

	        for (int ip = 0; ip < numberIntPoints; ip++){

	            int elemC = 0;
	            double xsiC[dim] = {};
	            double x_[dim] = {};

	            x_[0] = squad.interpolateQuadraticVariableIso(x1,ip,wpc,inc,IsoParFine,patch);
	            x_[1] = squad.interpolateQuadraticVariableIso(x2,ip,wpc,inc,IsoParFine,patch);

	            searchPointCorrespondence(x_,nodesCoarse_,elementsCoarse_,IsoParCoarse, elemTypeCoarse,
	                                      elementsCoarse_.size(),xsiC,elemC, 
	                                      elementsFine_[elementsGlueZoneFine_[i]] -> getIntegPointCorrespondenceElement_ISO(ip));

	            
	            elementsFine_[elementsGlueZoneFine_[i]] -> setIntegrationPointCorrespondence_ISO(ip,xsiC,elemC);
	        };      

        };
             
    };
};

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

            // double &alpha_f = parametersFine -> getAlphaF();
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

            if (elemTypeFine == 0){ //FEM element

	            if (side == 0){                 
	                bconnec[0] = connec[1];
	                bconnec[1] = connec[2];
	                bconnec[2] = connec[4];
	            };
	            if (side == 1){
	                bconnec[0] = connec[2];
	                bconnec[1] = connec[0];
	                bconnec[2] = connec[5];
				};
	            if (side == 2){
	                bconnec[0] = connec[0];
	                bconnec[1] = connec[1];
	                bconnec[2] = connec[3];              
	            };
	 
	            //first segment
	            int no1 = bconnec[0];
	            int no2 = bconnec[2];
	            	
	            double *xx1 = nodesFine_[no1] -> getCoordinates();
	            double *xx2 = nodesFine_[no2] -> getCoordinates();
                // double *xxp1 = nodesCoarse_[no1] -> getPreviousCoordinates();
                // double *xxp2 = nodesCoarse_[no2] -> getPreviousCoordinates();

	            for (int k = 0; k < dim; k++) {
	            	x1[k] = xx1[k];
	            	x2[k] = xx2[k];
                    // x1[k] = alpha_f * xx1[k] + (1. - alpha_f) * xxp1[k];
                    // x2[k] = alpha_f * xx2[k] + (1. - alpha_f) * xxp2[k];
	            }
	            
	            double sLength = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
	                                  (x1[0] - x2[0]) * (x1[0] - x2[0]));
	            
	            n[0] = ((x2[1] - x1[1]) / sLength);
	            n[1] = ((x1[0] - x2[0]) / sLength);

	            nodesFine_[no1] -> setInnerNormal(n);
	            nodesFine_[no2] -> setInnerNormal(n);

	            //second segment
	            no1 = bconnec[2];
	            no2 = bconnec[1];
	            
	            xx1 = nodesFine_[no1] -> getCoordinates();
	            xx2 = nodesFine_[no2] -> getCoordinates();
                //xxp1 = nodesCoarse_[no1] -> getPreviousCoordinates();
                //xxp2 = nodesCoarse_[no2] -> getPreviousCoordinates();

	            for (int k = 0; k < dim; k++) {
	            	x1[k] = xx1[k];
	            	x2[k] = xx2[k];
                    // x1[k] = alpha_f * xx1[k] + (1. - alpha_f) * xxp1[k];
                    // x2[k] = alpha_f * xx2[k] + (1. - alpha_f) * xxp2[k];
	            }
	            
	            sLength = sqrt((x2[1] - x1[1]) * (x2[1] - x1[1]) +
	                           (x1[0] - x2[0]) * (x1[0] - x2[0]));

	            n[0] = ((x2[1] - x1[1]) / sLength);
	            n[1] = ((x1[0] - x2[0]) / sLength);

	            nodesFine_[no1] -> setInnerNormal(n);
	            nodesFine_[no2] -> setInnerNormal(n); 

            } else {//IGA element
            	
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

                    // double *xx = nodesFine_[connec[i]]->getCoordinates();
                    // double *xxp = nodesFine_[connec[i]]->getPreviousCoordinates();
                    // for (int j = 0; j < dim; j++) coord_[i][j] = alpha_f * xx[j] + (1. - alpha_f) * xxp[j];
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
            };
        };//constrain == 2
    }; //numBoundFine






    // lembrar que a normal do ultimo elemento com aquele ponto de controle (bezier control point) 
    // vai me dar a normal final do ponto de controle

    //Gambiarra para os cantos para o problema da cavidade com todas as bordas com malhas refinidas e malha isogeométrica fine
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


    

    //Gambiarra para malha FEM fine
    n[0] = 1.;
    n[1] = 1.;

    nodesFine_[6] -> setInnerNormal(n);

    n[0] = -1.;
    n[1] = 1.;

    nodesFine_[8] -> setInnerNormal(n);

    n[0] = -1.;
    n[1] = -1.;

    nodesFine_[17] -> setInnerNormal(n);
    
    n[0] = 1.;
    n[1] = -1.;

    nodesFine_[15] -> setInnerNormal(n);



    //Coarse mesh - closer distance for nodes or control points from defined fine boundary 
    for (int ino = 0; ino < numNodesCoarse; ino++){
        
        double *xx = nodesCoarse_[ino]->getCoordinates();
        for (int i = 0; i < dim; i++) x[i] = xx[i];

        // double &alphaC_f = parametersCoarse -> getAlphaF();
        // double *xx = nodesCoarse_[ino]->getCoordinates();
        // double *xxp = nodesCoarse_[ino]->getPreviousCoordinates();
        // for (int i = 0; i < dim; i++) x[i] = alphaC_f * xx[i] + (1. - alphaC_f) * xxp[i]; 
        
        dist=10000000000000000000000000000.;
        
        for (int ibound = 0; ibound < numBoundElemFine; ibound++){
            
            if (boundaryFine_[ibound] -> getConstrain(0) == 2){

                // double &alpha_f = parametersFine -> getAlphaF();
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

                if (elemTypeFine == 0) { //FEM MESH

	                if (side == 0){                 
		                bconnec[0] = connec[1];
		                bconnec[1] = connec[2];
		                bconnec[2] = connec[4];
		            };
		            if (side == 1){
		                bconnec[0] = connec[2];
		                bconnec[1] = connec[0];
		                bconnec[2] = connec[5];
					};
		            if (side == 2){
		                bconnec[0] = connec[0];
		                bconnec[1] = connec[1];
		                bconnec[2] = connec[3];              
		            };

		            //first segment
	                int no1,no2;
	                no1 = bconnec[0];
	                no2 = bconnec[2];
	                
	                double *xx1 = nodesFine_[no1] -> getCoordinates();
                    double *xx2 = nodesFine_[no2] -> getCoordinates();
                    // double *xxp1 = nodesCoarse_[no1] -> getPreviousCoordinates();
                    // double *xxp2 = nodesCoarse_[no2] -> getPreviousCoordinates();

                    for (int k = 0; k < dim; k++) {
                        x1[k] = xx1[k];
                        x2[k] = xx2[k];
                        // x1[k] = alpha_f * xx1[k] + (1. - alpha_f) * xxp1[k];
                        // x2[k] = alpha_f * xx2[k] + (1. - alpha_f) * xxp2[k];
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
	                
	                xx1 = nodesFine_[no1] -> getCoordinates();
                    xx2 = nodesFine_[no2] -> getCoordinates();
                    // double *xxp1 = nodesCoarse_[no1] -> getPreviousCoordinates();
                    // double *xxp2 = nodesCoarse_[no2] -> getPreviousCoordinates();

                    for (int k = 0; k < dim; k++) {
                        x1[k] = xx1[k];
                        x2[k] = xx2[k];
                        // x1[k] = alpha_f * xx1[k] + (1. - alpha_f) * xxp1[k];
                        // x2[k] = alpha_f * xx2[k] + (1. - alpha_f) * xxp2[k];
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

                } else { //IGA MESH

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
                        // double *xx = nodesFine_[connec[i]]->getCoordinates();
                        // double *xxp = nodesFine_[connec[i]]->getPreviousCoordinates();
                        // for (int j = 0; j < dim; j++) coord_[i][j] = alpha_f * xx[j] + (1. - alpha_f) * xxp[j];
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

                }; //If IGA mesh
            }; //if bf is the blend boundary
        }; // numboundaryfine

   
        if (fabs(nodesCoarse_[ino] -> getDistFunction()) < 1.e-2){

        	// //PARA O PROBLEMA DO CILINDRO FIZ ISSO PARA CORRIGIR KKKK
        	// dist *= -1.;
            nodesCoarse_[ino] -> setDistFunction(dist); 
       
        };
    
    }; //numnodescoarse


    //Fine mesh - closer distance for nodes or control points from defined fine boundary 
    for (int ino = 0; ino < numNodesFine; ino++){
        
        double &alpha_f = parametersFine -> getAlphaF();  

        double *xx=nodesFine_[ino]->getCoordinates();

        for (int i = 0; i < dim; i++) x[i] = xx[i];

        // double *xx = nodesFine_[ino]->getCoordinates();
        // double *xxp = nodesFine_[ino]->getPreviousCoordinates();
        // for (int i = 0; i < dim; i++) x[i] = alpha_f * xx[i] + (1. - alpha_f) * xxp[i]; 

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
                };

                if (elemTypeFine == 0){ //FEM mesh

                	if (side == 0){                 
		                bconnec[0] = connec[1];
		                bconnec[1] = connec[2];
		                bconnec[2] = connec[4];
		            };
		            if (side == 1){
		                bconnec[0] = connec[2];
		                bconnec[1] = connec[0];
		                bconnec[2] = connec[5];
					};
		            if (side == 2){
		                bconnec[0] = connec[0];
		                bconnec[1] = connec[1];
		                bconnec[2] = connec[3];              
		            };

		            //first segment
	                int no1,no2;
	                no1 = bconnec[0];
	                no2 = bconnec[2];
	                
	                double *xx1 = nodesFine_[no1] -> getCoordinates();
                    double *xx2 = nodesFine_[no2] -> getCoordinates();
                    // double *xxp1 = nodesCoarse_[no1] -> getPreviousCoordinates();
                    // double *xxp2 = nodesCoarse_[no2] -> getPreviousCoordinates();

                    for (int k = 0; k < dim; k++) {
                        x1[k] = xx1[k];
                        x2[k] = xx2[k];
                        // x1[k] = alpha_f * xx1[k] + (1. - alpha_f) * xxp1[k];
                        // x2[k] = alpha_f * xx2[k] + (1. - alpha_f) * xxp2[k];
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
	                
	                xx1 = nodesFine_[no1] -> getCoordinates();
                    xx2 = nodesFine_[no2] -> getCoordinates();
                    // double *xxp1 = nodesCoarse_[no1] -> getPreviousCoordinates();
                    // double *xxp2 = nodesCoarse_[no2] -> getPreviousCoordinates();

                    for (int k = 0; k < dim; k++) {
                        x1[k] = xx1[k];
                        x2[k] = xx2[k];
                        // x1[k] = alpha_f * xx1[k] + (1. - alpha_f) * xxp1[k];
                        // x2[k] = alpha_f * xx2[k] + (1. - alpha_f) * xxp2[k];
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

                } else { //IGA mesh
                	
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

                        // double *xx = nodesFine_[connec[i]]->getCoordinates();
                        // double *xxp = nodesFine_[connec[i]]->getPreviousCoordinates();
                        // for (int j = 0; j < dim; j++) coord_[i][j] = alpha_f * xx[j] + (1. - alpha_f) * xxp[j];
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
                };
            }; //if bf is the blend boundary  
        }; // numboundfine

        // //MUDANCA DE SINAL PARA PROBLEMA DO CILINDRO
        // dist *= -1.;
        if(dist < 0) dist = 0;
        nodesFine_[ino] -> setDistFunction(dist);

    };//numNodes

}; //função


//------------------------------------------------------------------------------
//------------------------------SETS GLUING ZONE--------------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::setGluingZone(){

    double glueZoneThickness = fineModel.glueZoneThickness;

	int dim = 2;
    int flag;
    int nodesCZ[numNodesFine];
    int nodesCZNC[NCNumberNodesF];
    int nodesCZ2[numNodesCoarse];
    QuadShapeFunction<2> shapeQuad;

    for (int i = 0; i < numNodesFine; i++) nodesCZ[i] = 0;   
    for (int i = 0; i < NCNumberNodesF; i++) nodesCZNC[i] = 0; 

    elementsGlueZoneFine_.reserve(numElemFine / 3);
    glueZoneFine_.reserve(numElemFine/3);
    nodesGlueZoneFine_.reserve(numNodesFine / 3);
    NCnodesGlueZoneFine_.reserve(NCNumberNodesF / 3);
    
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
       	 
       	if (elemTypeFine == 0) { //FEM elem

	        flag = 0;
	        for (int ino = 0; ino < 6; ino++){
                double dist = nodesFine_[connec[ino]] -> getDistFunction();
	            if ( dist <= glueZoneThickness + 0.00001){
	                flag += 1;
	            };
	        };

	        if (flag == 6) {
	            elementsGlueZoneFine_.push_back(jel);
	            elementsFine_[jel] -> setGlueZone();

	            GlueZone *el = new GlueZone(index++,jel);
	            glueZoneFine_.push_back(el);

	    	};
       	
       	}else { //IGA elem

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
       	}; //else IGA mesh         
    };

    //Defines which nodes are in the gluing zone
    numElemGlueZoneFine = elementsGlueZoneFine_.size();

    for (int i = 0; i < numElemGlueZoneFine; i++){

        int *connec = elementsFine_[elementsGlueZoneFine_[i]] -> getConnectivity();
        
        if (elemTypeFine == 0){ //FEM mesh
        	for (int ino = 0; ino < 6; ino++){
	            nodesCZ[connec[ino]] += 1;
	            nodesFine_[connec[ino]] -> setGlueZone();
	        };
        } else { //IGA mesh
        	for (int ino = 0; ino < 9; ino++){
               
                nodesCZNC[nodesFine_[connec[ino]] -> getnewcon()] += 1;
	            nodesCZ[connec[ino]] += 1;

	            nodesFine_[connec[ino]] -> setGlueZone();
	        };
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



    //Save the number of non-coincidents nodes in the gluing zone (new numeration and connectivity)
    if (elemTypeFine == 0){
        NCnumNodesGlueZoneFine = numNodesGlueZoneFine;
    } else {
        NCnumNodesGlueZoneFine = 0;
        for (int i = 0; i < NCNumberNodesF; i++){
            if(nodesCZNC[i] > 0) {
                NCnumNodesGlueZoneFine += 1;
                NCnodesGlueZoneFine_.push_back(i);
            };
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
        
        if (elemTypeFine == 0) { //FEM mesh

        	int *connecAux= new int[6];
	        
	        int *connec = elementsFine_[elementsGlueZoneFine_[i]] -> getConnectivity();
	   
	        for (int ino = 0; ino < numNodesGlueZoneFine; ino++)
	            for (int k = 0; k < 6; k++)
	                if (nodesGlueZoneFine_[ino] == connec[k]) connecAux[k] = ino;
	        
	        glueZoneFine_[i] -> setConnectivity(connecAux);

        } else { //IGA mesh

        	int *connecAux = new int[9];
	        
	        int *connec = elementsFine_[elementsGlueZoneFine_[i]] -> getConnectivity();
	       
	        for (int ino = 0; ino < numNodesGlueZoneFine; ino++){
                for (int k = 0; k < 9; k++){
                    if (nodesGlueZoneFine_[ino] == connec[k]) {
                        connecAux[k] = ino;
                    };
                };
      
            };
	        
	        glueZoneFine_[i] -> setConnectivity(connecAux);          

            //New connectivity without repeated nodes
            int *connecAuxNew = new int[9];
           
            for (int ino = 0; ino < NCnumNodesGlueZoneFine; ino++){
                for (int k = 0; k < 9; k++){
                    if (NCnodesGlueZoneFine_[ino] == nodesFine_[connec[k]] -> getnewcon()) {
                        connecAuxNew[k] = ino;
                    };
                };
            };
            
            glueZoneFine_[i] -> setNewConnectivity(connecAuxNew);

        };
    };

    
    // Defines a criterion to select the coarse elements that are in the gluing zone
    for (int jel = 0; jel < numElemCoarse; jel++){
        
        int *connec = elementsCoarse_[jel] -> getConnectivity();
        
        flag = 0;
        
        if (elemTypeCoarse == 0) { //FEM mesh

            for (int ino = 0; ino < 6; ino++){
                double dist = nodesCoarse_[connec[ino]]->getDistFunction();
                if ((dist  <= glueZoneThickness + 0.00001) && (dist   >= 0.00001)){
                     flag += 1;
                };
            };
            if (flag != 0) {
                elementsGlueZoneCoarse_.push_back(jel);
                elementsCoarse_[jel] -> setGlueZone();
            };


        } else { //IGA mesh

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
        };
    };


    //Defines which coarse nodes are in the gluing zone
    numElemGlueZoneCoarse = elementsGlueZoneCoarse_.size();
    for (int i = 0; i < numElemGlueZoneCoarse; i++){
        int *connec = elementsCoarse_[elementsGlueZoneCoarse_[i]] -> getConnectivity();
        
        if (elemTypeCoarse == 0) { //FEM mesh
            for (int ino = 0; ino < 6; ino++){
                nodesCZ2[connec[ino]] += 1;
                nodesCoarse_[connec[ino]] -> setGlueZone();
            };
        } else { //IGA mesh
            for (int ino = 0; ino < 9; ino++){
                nodesCZ2[connec[ino]] += 1;
                nodesCoarse_[connec[ino]] -> setGlueZone();
            };
        };
            
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

    double glueZoneThickness = fineModel.glueZoneThickness;
    double arlequinEpsilon = fineModel.arlequinEpsilon;

	glueZoneThickness *= 1.01;	// Thickness from gluing zone
	
    //IGA or FEM coarse mesh
    for (int iNode = 0; iNode < numNodesCoarse; iNode++){

    	double r = nodesCoarse_[iNode] -> getDistFunction();
        
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

    if (elemTypeCoarse == 0) {
        for (int jel = 0; jel < numElemCoarse; jel++){
            elementsCoarse_[jel] -> setIntegPointWeightFunction_FEM();        
        };
    } else {
        for (int jel = 0; jel < numElemCoarse; jel++){
            elementsCoarse_[jel] -> setIntegPointWeightFunction_ISO();        
        };
    };
        


    //IGA or FEM FINE MESH
    for (int iNode = 0; iNode < numNodesFine; iNode++){

    	double r = nodesFine_[iNode] -> getDistFunction();
        
		if (r >= glueZoneThickness){
        	wFuncValue = 1. - arlequinEpsilon;
    	} else {
        	wFuncValue = (1. - arlequinEpsilon) / glueZoneThickness * r;
        	if (wFuncValue > (1. - arlequinEpsilon)) wFuncValue = 1. - arlequinEpsilon;

    	};       
        nodesFine_[iNode] -> setWeightFunction(wFuncValue);
    };  //ielem

    if (elemTypeFine == 0) { //FEM mesh
    	for (int jel = 0; jel < numElemFine; jel++){
        	elementsFine_[jel] -> setIntegPointWeightFunction_FEM();        
    	};

    } else { //IGA mesh
    	for (int jel = 0; jel < numElemFine; jel++){
        	elementsFine_[jel] -> setIntegPointWeightFunction_ISO();        
    	};
    };
      
     
    return;

};


//------------------------------------------------------------------------------
//------------------------PRINT RESULTS IN PARAVIEW-----------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::printResults(int step) {

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);

    if (step % fineModel.printFreq == 0){

        if (rank == 0){

            std::string result;
            std::ostringstream convert;
            convert << step+100000;
            result = convert.str();

            //COARSE MESH
            int dim = 2;
            std::string s = "COARSEoutput"+result+".vtu";
            std::fstream output_v(s.c_str(), std::ios_base::out);

            if (elemTypeCoarse == 0) { //FEM mesh

                output_v << "<?xml version=\"1.0\"?>" << std::endl
                        << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                        << "  <UnstructuredGrid>" << std::endl
                        << "  <Piece NumberOfPoints=\"" << numNodesCoarse
                        << "\"  NumberOfCells=\"" << numElemCoarse
                        << "\">" << std::endl;

                //WRITE NODAL COORDINATES
                output_v << "    <Points>" << std::endl
                        << "      <DataArray type=\"Float64\" "
                        << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

                        
                for (int i = 0; i< numNodesCoarse;i++){
                    output_v << nodesCoarse_[i] -> getCoordinateValue(0) << " " 
                            << nodesCoarse_[i] -> getCoordinateValue(1) << " " << 0. << std::endl;
                }

                output_v << "      </DataArray>" << std::endl
                        << "    </Points>" << std::endl;
                
                //WRITE ELEMENT CONNECTIVITY
                output_v << "    <Cells>" << std::endl
                        << "      <DataArray type=\"Int32\" "
                        << "Name=\"connectivity\" format=\"ascii\">" << std::endl;

                for (int i = 0; i < numElemCoarse; i++){
                    int *connec = elementsCoarse_[i] -> getConnectivity();
                    int con[6];
                    for (int i = 0; i < 6; i++) con[i] = connec[i];
                    output_v << con[0] << " " << con[1] << " " << con[2] << " " 
                            << con[3] << " " << con[4] << " " << con[5] << std::endl;
                };
                output_v << "      </DataArray>" << std::endl;
            
                //WRITE OFFSETS IN DATA ARRAY
                output_v << "      <DataArray type=\"Int32\""
                        << " Name=\"offsets\" format=\"ascii\">" << std::endl;
                int aux = 0;
                for (int i = 0; i < numElemCoarse; i++){
                    output_v << aux + 6 << std::endl;
                    aux += 6;
                };
                output_v << "      </DataArray>" << std::endl;
            
                //WRITE ELEMENT TYPES
                output_v << "      <DataArray type=\"UInt8\" Name=\"types\" "
                        << "format=\"ascii\">" << std::endl;
                
                for (int i = 0; i < numElemCoarse; i++){
                    output_v << 22 << std::endl;
                };
                output_v << "      </DataArray>" << std::endl
                        << "    </Cells>" << std::endl;

                // WRITE NODAL RESULTS
                output_v << "    <PointData>" << std::endl;

                if (coarseModel.printVelocity){
                    output_v<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                            << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                    for (int i=0; i<numNodesCoarse; i++){
                        output_v << nodesCoarse_[i] -> getVelocity(0) << " "              
                                << nodesCoarse_[i] -> getVelocity(1) << " " 
                                << 0. << std::endl;
                    };
                    output_v << "      </DataArray> " << std::endl;
                };

                if (coarseModel.printRealVelocity){
                    output_v<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                            << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                    for (int i=0; i<numNodesCoarse; i++){
                        output_v << nodesCoarse_[i] -> getVelocityArlequin(0)<< " "              
                                << nodesCoarse_[i] -> getVelocityArlequin(1) << " " 
                                << 0. << std::endl;
                    };
                    output_v << "      </DataArray> " << std::endl;
                };


                if (coarseModel.printPressure){
                    output_v<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                            << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
                    for (int i=0; i<numNodesCoarse; i++){
                        output_v << 0. << " " << 0. << " " 
                                << nodesCoarse_[i] -> getPressure() << std::endl;
                    };
                    output_v << "      </DataArray> " << std::endl;
                };

                if (coarseModel.printRealPressure){
                    output_v<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                            << "Name=\"Real Pressure\" format=\"ascii\">" << std::endl;
                    for (int i=0; i<numNodesCoarse; i++){
                        output_v << 0. << " " << 0. << " " 
                                << nodesCoarse_[i] -> getPressureArlequin() << std::endl;
                    };
                    output_v << "      </DataArray> " << std::endl;
                };


                if (coarseModel.printDistFunction){
                    output_v<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                            << "Name=\"Dist Function\" format=\"ascii\">" << std::endl;
                    for (int i=0; i<numNodesCoarse; i++){
                        output_v << nodesCoarse_[i] -> getDistFunction() << std::endl;
                    };
                    output_v << "      </DataArray> " << std::endl;
                };

                if (coarseModel.printEnergyWeightFunction){
                    output_v<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                            << "Name=\"Energy Weight Function\" format=\"ascii\">" << std::endl;
                    for (int i=0; i<numNodesCoarse; i++){
                        output_v << nodesCoarse_[i] -> getWeightFunction()<< std::endl;
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

            } else { //IGA mesh


                //Interpolated variables for IGA analysis (Bezier transformation)
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

            }; //elseIGA

                


        
            //PRINT FINE MODEL RESULTS - IGA MESH or FEM mesh
            std::string f = "FINEoutput"+result+".vtu";
            
            std::fstream output_vf(f.c_str(), std::ios_base::out);

            if (elemTypeFine == 0){ //FEM elements

                output_vf << "<?xml version=\"1.0\"?>" << std::endl
                << "<VTKFile type=\"UnstructuredGrid\">" << std::endl
                << "  <UnstructuredGrid>" << std::endl
                << "  <Piece NumberOfPoints=\"" << numNodesFine
                << "\"  NumberOfCells=\"" << numElemFine
                << "\">" << std::endl;

                //WRITE NODAL COORDINATES
                output_vf << "    <Points>" << std::endl
                        << "      <DataArray type=\"Float64\" "
                        << "NumberOfComponents=\"3\" format=\"ascii\">" << std::endl;

                for (int i = 0; i < numNodesFine; i++){
                    double *x = nodesFine_[i] -> getCoordinates();
                    double x_[dim];
                    for (int i = 0; i < dim; i++) x_[i] = x[i];
                    output_vf << x_[0] << " " << x_[1] << " " << 0.1 << std::endl;
                };

                output_vf << "      </DataArray>" << std::endl
                        << "    </Points>" << std::endl;
            
                //WRITE ELEMENT CONNECTIVITY
                output_vf << "    <Cells>" << std::endl
                        << "      <DataArray type=\"Int32\" "
                        << "Name=\"connectivity\" format=\"ascii\">" << std::endl;
                
                for (int i = 0; i < numElemFine; i++){
                    int *connec = elementsFine_[i] -> getConnectivity();
                    int con[6];
                    for (int i = 0; i < 6; i++) con[i] = connec[i];
                    output_vf << con[0] << " " << con[1] << " " << con[2] << " " 
                            << con[3] << " " << con[4] << " " << con[5] << std::endl;
                };
                output_vf << "      </DataArray>" << std::endl;
            
                //WRITE OFFSETS IN DATA ARRAY
                output_vf << "      <DataArray type=\"Int32\""
                        << " Name=\"offsets\" format=\"ascii\">" << std::endl;
                
                int aux = 0;
                for (int i = 0; i < numElemFine; i++){
                    output_vf << aux + 6 << std::endl;
                    aux += 6;
                };
                output_vf << "      </DataArray>" << std::endl;
            
                //WRITE ELEMENT TYPES
                output_vf << "      <DataArray type=\"UInt8\" Name=\"types\" "
                        << "format=\"ascii\">" << std::endl;
            
                for (int i = 0; i < numElemFine; i++){
                    output_vf << 22 << std::endl;
                };

                output_vf << "      </DataArray>" << std::endl
                        << "    </Cells>" << std::endl;

                //WRITE NODAL RESULTS
                output_vf << "    <PointData>" << std::endl;

                if (fineModel.printVelocity){
                    output_vf <<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                            << "Name=\"Velocity\" format=\"ascii\">" << std::endl;
                    for (int i=0; i<numNodesFine; i++){ 
                            output_vf << nodesFine_[i] -> getVelocity(0) << " "              
                                    << nodesFine_[i] -> getVelocity(1)  << " " 
                                    << 0. << std::endl;
                    };
                    output_vf << "      </DataArray> " << std::endl;
                };

                if (fineModel. printRealVelocity){
                    output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                            << "Name=\"Real Velocity\" format=\"ascii\">" << std::endl;
                    for (int i=0; i<numNodesFine; i++){
                            output_vf << nodesFine_[i] -> getVelocityArlequin(0) << " "              
                                    << nodesFine_[i] -> getVelocityArlequin(1)  << " " 
                                    << 0. << std::endl;
                    };
                    output_vf << "      </DataArray> " << std::endl;
                };
        
            
                if (fineModel.printPressure){
                    output_vf <<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                            << "Name=\"Pressure\" format=\"ascii\">" << std::endl;
                    for (int i=0; i<numNodesFine; i++){
                        output_vf << 0. << " " << 0. << " " 
                                << nodesFine_[i] -> getPressure()<< std::endl;
                    };
                    output_vf << "      </DataArray> " << std::endl;
                };

                if (fineModel.printRealPressure){
                    output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                            << "Name=\"Real Pressure\" format=\"ascii\">" << std::endl;
                    for (int i=0; i<numNodesFine; i++){
                        output_vf << 0. << " " << 0. << " " 
                                << nodesFine_[i] -> getPressureArlequin()<< std::endl;
                    };
                    output_vf << "      </DataArray> " << std::endl;
                };

                if (fineModel.printLagrangeMultipliers){
                    output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                            << "Name=\"Lagrange Multipliers\" format=\"ascii\">" << std::endl;
                    for (int i=0; i<numNodesFine; i++){
                        output_vf << nodesFine_[i] -> getLagrangeMultiplier(0) <<  " "
                                << nodesFine_[i] -> getLagrangeMultiplier(1)<< " " 
                                << 0.0 << " " << std::endl;
                    };
                    output_vf << "      </DataArray> " << std::endl;
                };


                if (fineModel.printDistFunction){
                    output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                            << "Name=\"printDistFunction\" format=\"ascii\">" << std::endl;
                    for (int i=0; i<numNodesFine; i++){
                        output_vf << nodesFine_[i] -> getDistFunction() << std::endl;
                    };
                    output_vf << "      </DataArray> " << std::endl;
                };


                output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"3\" "
                        << "Name=\"Normal\" format=\"ascii\">" << std::endl;
                for (int i=0; i<numNodesFine; i++){

                    double *n = nodesFine_[i] -> getInnerNormal();
                    output_vf << n[0] <<  " "
                            << n[1] << " " 
                            << 0.0 << " " << std::endl;
                };
                output_vf << "      </DataArray> " << std::endl;

                

                if (fineModel.printEnergyWeightFunction){
                    output_vf<<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                            << "Name=\"EnergyWeightFunction\" format=\"ascii\">" << std::endl;
                    for (int i=0; i<numNodesFine; i++){
                        output_vf << nodesFine_[i] -> getWeightFunction() << std::endl;
                    };
                    output_vf << "      </DataArray> " << std::endl;
                };


                output_vf << "    </PointData>" << std::endl; 

                //WRITE ELEMENT RESULTS
                output_vf << "    <CellData>" << std::endl;
                
                if (fineModel.printProcess){
                    output_vf <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                            << "Name=\"Process\" format=\"ascii\">" << std::endl;
                    for (int i=0; i<numElemFine; i++){
                        output_vf << domDecompFine.first[i] << std::endl;
                    };
                    output_vf << "      </DataArray> " << std::endl;
                };

                int cont=0;
                if (fineModel.printGlueZone){
                    output_vf <<"      <DataArray type=\"Float64\" NumberOfComponents=\"1\" "
                            << "Name=\"Glue Zone\" format=\"ascii\">" << std::endl;
                    cont = 0;
                    for (int i=0; i<numElemFine; i++){
                        if (elementsGlueZoneFine_[cont] == i){
                            output_vf << 1.0 << std::endl;
                            cont++;
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

            } else { //IGA elements

                int numBezierNodes = fineModel.NumBezierNodes;

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

                //Interpolated variables for IGA analysis (Bezier transformation)
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
                
                int aux = 0;
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
            }; //else IGA

        };//rank == 0
    }; //printfreq

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

        int numberIntPoints;
        if (elemTypeFine == 0){
        	numberIntPoints = elementsFine_[0] -> getNumberOfIntegrationPointsSpecial_FEM();
        } else {
        	numberIntPoints = elementsFine_[0] -> getNumberOfIntegrationPointsSpecial_ISO();
        }

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
            	double qxsiC[2];
            	int indCoarseElem;

            	if (elemTypeFine == 0){ //FEM fine mesh

            		//correspondent integration point in the coarse mesh
	            	qxsiC[0] = elementsFine_[elementsGlueZoneFine_[i]] -> getIntegPointCoordinatesValue_FEM(ip,0);
	            	qxsiC[1] = elementsFine_[elementsGlueZoneFine_[i]] -> getIntegPointCoordinatesValue_FEM(ip,1);
	            	//coarse element index
            		indCoarseElem = elementsFine_[elementsGlueZoneFine_[i]] -> getIntegPointCorrespondenceElement_FEM(ip);
	            	
            	
            	} else { //IGA fine mesh

            		//correspondent integration point in the coarse mesh
	            	qxsiC[0] = elementsFine_[elementsGlueZoneFine_[i]] -> getIntegPointCoordinatesValue_ISO(ip,0);
	            	qxsiC[1] = elementsFine_[elementsGlueZoneFine_[i]] -> getIntegPointCoordinatesValue_ISO(ip,1);
	            	//coarse element index
	            	indCoarseElem = elementsFine_[elementsGlueZoneFine_[i]] -> getIntegPointCorrespondenceElement_ISO(ip);

            	};


                int *connec = elementsCoarse_[indCoarseElem] -> getConnectivity();

                if (elemTypeCoarse == 0) { //FEM mesh

                    double phi_[6];
                    shapeQuad.evaluateFem(qxsiC,phi_);
                                            
                    double coord[2] = {};
                    for (int j = 0; j < 6; j++){
                        double *x = nodesCoarse_[connec[j]] -> getCoordinates();
                        coord[0] += x[0]*phi_[j];
                        coord[1] += x[1]*phi_[j];
                    }

                    output_v << coord[0] << " " << coord[1] << " " << 0.2 << std::endl;

                } else { //IGA mesh

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

                    output_v << coord[0] << " " << coord[1] << " " << 0.2 << std::endl;


                }

                
	            	
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
        


    	//PRINT FINE MODEL RESULTS (IGA or FEM mesh)
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

        	if (elemTypeFine == 0) { //FEM mesh

        		
        		double x1[6], x2[6];

				int *connec = elementsFine_[elementsGlueZoneFine_[i]] -> getConnectivity();

                	
            	for (int j = 0; j < 6; j++){
                	double *x = nodesFine_[connec[j]] -> getCoordinates();
                	x1[j] = x[0];
                	x2[j] = x[1];
            	}; 

            	for (int ip = 0; ip < numberIntPoints; ip++){

            		double x_[2];
                	x_[0] = squad.interpolateQuadraticVariableFem(x1,ip);
                	x_[1] = squad.interpolateQuadraticVariableFem(x2,ip);
            
					output_vf << x_[0] << " " << x_[1] << " " << 0.2 << std::endl;
            	}



        	} else { //IGA mesh

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

            	for (int ip = 0; ip < numberIntPoints; ip++){

            		double x_[2];
                	x_[0] = squad.interpolateQuadraticVariableIso(x1,ip,wpc,inc,IsoParFine,patch);
                	x_[1] = squad.interpolateQuadraticVariableIso(x2,ip,wpc,inc,IsoParFine,patch);
            
					output_vf << x_[0] << " " << x_[1] << " " << 0.0 << std::endl;
            	}
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

template<>
void Arlequin<2>::dragAndLiftCoefficientsFEM(std::ofstream& dragLift, int& iTimeStep){

    double dragCoefficient = 0.;
    double liftCoefficient = 0.;
    double pressureDragCoefficient = 0.;
    double pressureLiftCoefficient = 0.;
    double frictionDragCoefficient = 0.;
    double frictionLiftCoefficient = 0.;
    double pitchingMomentCoefficient = 0.;
    double pMom = 0.;
    double per = 0.;
    double& velocityInf = parametersFine -> getVelocityInf(0);
    double& rhoInf = parametersFine ->getDensity();
    double& dTime = parametersFine -> getTimeStep();

    
    for (int jel = 0; jel < numBoundElemFine; jel++){   
        
        double dForce = 0.;
        double lForce = 0.;
        double pDForce = 0.;
        double pLForce = 0.;
        double fDForce = 0.;
        double fLForce = 0.;
        double aux_Mom = 0.;
        double aux_Per = 0.;
        
        for (int i=0; i< fineModel.numberOfLines; i++){
            
            if (boundaryFine_[jel] -> getBoundaryGroup() == fineModel.dragAndLiftBoundary[i]){            
                
                int iel = boundaryFine_[jel] -> getElement();

                //Recognizing the element side
                std::pair<std::vector<int>, std::vector<int> > elemBound;
                elemBound = elementsFine_[iel] -> getElemSideInBoundary();
                int numofBoundaries = elemBound.first.size(); 
                int side;
                for (int j = 0; j < numofBoundaries; j++) {
                    if (elemBound.second[j] == fineModel.dragAndLiftBoundary[i]) {
                        side = elemBound.first[j];
                    };
                };
                            
                elementsFine_[iel] -> computeDragAndLiftForces_FEM(side, pDForce, pLForce, fDForce, fLForce, 
                                                                   dForce, lForce, aux_Mom, aux_Per);

                pMom += aux_Mom;
                per += aux_Per;

            };
        };
        
        pressureDragCoefficient += pDForce / 
            (0.5 * rhoInf * velocityInf * velocityInf * 1);
        pressureLiftCoefficient += pLForce / 
            (0.5 * rhoInf * velocityInf * velocityInf * 1);
        frictionDragCoefficient += fDForce / 
            (0.5 * rhoInf * velocityInf * velocityInf * 1);
        frictionLiftCoefficient += fLForce / 
            (0.5 * rhoInf * velocityInf * velocityInf * 1);
        dragCoefficient += dForce / 
            (0.5 * rhoInf * velocityInf * velocityInf * 1);
        liftCoefficient += lForce / 
            (0.5 * rhoInf * velocityInf * velocityInf * 1);
    };

    pitchingMomentCoefficient = pMom / (rhoInf * velocityInf * velocityInf * per);
    
    if (rank == 0) {
        const int timeWidth = 15;
        const int numWidth = 15;
        dragLift << std::setprecision(5) << std::scientific;
        dragLift << std::left << std::setw(timeWidth) << iTimeStep * dTime;
        dragLift << std::setw(numWidth) << pressureDragCoefficient;
        dragLift << std::setw(numWidth) << pressureLiftCoefficient;
        dragLift << std::setw(numWidth) << frictionDragCoefficient;
        dragLift << std::setw(numWidth) << frictionLiftCoefficient;
        dragLift << std::setw(numWidth) << dragCoefficient;
        dragLift << std::setw(numWidth) << liftCoefficient;
        dragLift << std::setw(numWidth) << pitchingMomentCoefficient;
        dragLift << std::endl;
    };


return;
    
};

template<>
void Arlequin<2>::dragAndLiftCoefficientsISO(std::ofstream& dragLift, int &iTimeStep){

    
    double dragCoefficient = 0.;
    double liftCoefficient = 0.;
    double pressureDragCoefficient = 0.;
    double pressureLiftCoefficient = 0.;
    double frictionDragCoefficient = 0.;
    double frictionLiftCoefficient = 0.;
    double pitchingMomentCoefficient = 0.;
    double pMom = 0.;
    double per = 0.;
    double& velocityInf = parametersFine -> getVelocityInf(0);
    double& dTime = parametersFine -> getTimeStep();
    double& rhoInf = parametersFine ->getDensity();

    
    for (int jel = 0; jel < numBoundElemFine; jel++){   
        
        double dForce = 0.;
        double lForce = 0.;
        double pDForce = 0.;
        double pLForce = 0.;
        double fDForce = 0.;
        double fLForce = 0.;
        double aux_Mom = 0.;
        double aux_Per = 0.;
        
        for (int i=0; i< fineModel.numberOfLines; i++){
            
            if (boundaryFine_[jel] -> getBoundaryGroup() == fineModel.dragAndLiftBoundary[i]){            
                int iel = boundaryFine_[jel] -> getElement();

                //Recognizing the element side
                std::pair<std::vector<int>, std::vector<int> > elemBound;
                elemBound = elementsFine_[iel] -> getElemSideInBoundary();
                int numofBoundaries = elemBound.first.size(); 
                int side;
                for (int j = 0; j < numofBoundaries; j++) {
                    if (elemBound.second[j] == fineModel.dragAndLiftBoundary[i]) {
                        side = elemBound.first[j];
                    };
                }
                            
                elementsFine_[iel] -> computeDragAndLiftForces_ISO(side, pDForce, pLForce, fDForce, fLForce, 
                                                             dForce, lForce, aux_Mom, aux_Per);
                
                pMom += aux_Mom;
                per += aux_Per;

            };
        };
        
        pressureDragCoefficient += pDForce / 
            (0.5 * rhoInf * velocityInf * velocityInf * 1);
        pressureLiftCoefficient += pLForce / 
            (0.5 * rhoInf * velocityInf * velocityInf * 1);
        frictionDragCoefficient += fDForce / 
            (0.5 * rhoInf * velocityInf * velocityInf * 1);
        frictionLiftCoefficient += fLForce / 
            (0.5 * rhoInf * velocityInf * velocityInf * 1);
        dragCoefficient += dForce / 
            (0.5 * rhoInf * velocityInf * velocityInf * 1);
        liftCoefficient += lForce / 
            (0.5 * rhoInf * velocityInf * velocityInf * 1);
    };

    
    pitchingMomentCoefficient = pMom / (rhoInf * velocityInf * velocityInf * per);
    
    if (rank == 0) {
        const int timeWidth = 15;
        const int numWidth = 15;
        dragLift << std::setprecision(5) << std::scientific;
        dragLift << std::left << std::setw(timeWidth) << iTimeStep * dTime;
        dragLift << std::setw(numWidth) << pressureDragCoefficient;
        dragLift << std::setw(numWidth) << pressureLiftCoefficient;
        dragLift << std::setw(numWidth) << frictionDragCoefficient;
        dragLift << std::setw(numWidth) << frictionLiftCoefficient;
        dragLift << std::setw(numWidth) << dragCoefficient;
        dragLift << std::setw(numWidth) << liftCoefficient;
        dragLift << std::setw(numWidth) << pitchingMomentCoefficient;
        dragLift << std::endl;
    };

return;
    
};


//------------------------------------------------------------------------------
//---------------------------SETS FLUID MODELS----------------------------------
//------------------------------------------------------------------------------
template<>
void Arlequin<2>::setFluidModels(FluidMesh& coarse, FluidMesh& fine){

    coarseModel = coarse;
    fineModel = fine;

    //Gets Fine and Coarse models basic information from fluid data
    nodesCoarse_  = coarseModel.nodes_;
    nodesFine_    = fineModel.nodes_;

    elementsCoarse_ = coarseModel.elements_;
    elementsFine_   = fineModel.elements_;
 
    boundaryCoarse_ = coarseModel.boundary_;
    boundaryFine_   = fineModel.boundary_;

    numElemCoarse = coarseModel.elements_.size();
    numElemFine   = fineModel.elements_.size();
    numNodesCoarse = coarseModel.nodes_.size();
    numNodesFine   = fineModel.nodes_.size();

    numBoundElemFine = boundaryFine_.size();
    numBoundElemCoarse = boundaryCoarse_.size();

    domDecompCoarse = coarseModel.getDomainDecomposition();
    domDecompFine = fineModel.getDomainDecomposition();

    parametersCoarse = elementsCoarse_[0] -> getFluidParameters();
    parametersFine = elementsFine_[0] -> getFluidParameters();

    IsoParFine = fineModel.IsoPar_;
   	IsoParCoarse = coarseModel.IsoPar_;

    //starting the program with maximum dissipation
    double &IS = parametersFine -> getSpectralRadius();
    integScheme = IS;
    double dd = 0.;
    parametersCoarse -> setSpectralRadius(dd);
    parametersFine -> setSpectralRadius(dd);

    //Defines coarse elements type
    elemTypeCoarse = elementsCoarse_[0] -> getElemType();

    //Defines fine elements type
    elemTypeFine = elementsFine_[0] -> getElemType();

    //Non coincidente number controlPoints in fine and coarse mesh
    //Fine mesh (IGA or FEM)
    if (elemTypeFine == 0){ //FEM coarse mesh
    	NCNumberNodesF = numNodesFine;
    } else { // IGA coarse mesh
    	NCNumberNodesF = fineModel.NCNumberNodes;
    };

    //Coarse mesh (IGA or FEM)
    if (elemTypeCoarse == 0){ //FEM mesh
        NCNumberNodesC= numNodesCoarse;
    } else { //IGA mesh
        NCNumberNodesC= coarseModel.NCNumberNodes;
    };
   
    //Defines fine model as true
    for (int i=0; i < numElemFine; i++){
        elementsFine_[i] -> setModel(true); 
    };
    //Defines coarse model elements as false
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

    //Coarse mesh (FEM or IGA elements)
    for (int ielem = 0; ielem < numElemCoarse; ielem++){
    	
        int *connec = elementsCoarse_[ielem] -> getConnectivity();
    	
        if (elemTypeCoarse == 0) { //FEM mesh

            for (int i = 0; i < 6 ; i++){
                // velocity constrain
                for (int j = 0; j < dim; j++){
                    int constrain = nodesCoarse_[connec[i]] -> getConstrains(j);
                    if ((constrain == 1) || (constrain == 3)){
                        dofTemp.push_back(connec[i]*dim + j);
                    };
                };  

                // if (connec[i] == 2){
                //  dofTemp.push_back(2*NCNumberNodesC + connec[i]);
                // }               
            };

        } else { //IGA mesh

            for (int i = 0; i < 9 ; i++){
                int newconi = nodesCoarse_[connec[i]] -> getnewcon();
                // velocity constrain
                for (int j = 0; j < dim; j++){
                    int constrain = nodesCoarse_[connec[i]] -> getConstrains(j);
                    if ((constrain == 1) || (constrain == 3)){
                        dofTemp.push_back(newconi*dim + j);
                    };
                };
             //    if (connec[i] == 288){
             //     dofTemp.push_back(2*NCNumberNodesC + newconi);
            	// };                
            };

            

        };

        
            

    };

    //Fine mesh (FEM or IGA elements)
    for (int ielem = 0; ielem < numElemFine; ielem++){
        
        int *connec = elementsFine_[ielem] -> getConnectivity();
        
        if (elemTypeFine == 0) { //FEM mesh

        	for (int i = 0; i < 6 ; i++){
	            // velocity constrain
	            for (int j = 0; j < dim; j++){
	            	int constrain = nodesFine_[connec[i]] -> getConstrains(j);
	            	if ((constrain == 1) || (constrain == 3)){
	                	dofTemp.push_back(connec[i]*dim + (dim+1)*NCNumberNodesC + j);
	            	};
	            };
	            //cavity pressure constrain
	            // if (connec[i] == 23){
	            //  dofTemp.push_back(3*NCNumberNodesC + 2*NCNumberNodesF + connec[i]);
	            // }  
	        };

        } else {//IGA mesh

        	for (int i = 0; i < 9 ; i++){
	            int newconi = nodesFine_[connec[i]] -> getnewcon();
	            // velocity constrain
	            for (int j = 0; j < dim; j++){
	            	int constrain = nodesFine_[connec[i]] -> getConstrains(j);
	            	if ((constrain == 1) || (constrain == 3)){
	                	dofTemp.push_back(newconi*dim + (dim+1)*NCNumberNodesC + j);
	            	};
	            };
	           // cavity pressure constrain
	            // if (connec[i] == 1439){
	            //  dofTemp.push_back(3*NCNumberNodesC + 2*NCNumberNodesF + newconi);
	            // }  
	        };
        };
	        
    };
};

template<>
void Arlequin<2>::setMatVecValuesCoarseFEM(){


    for (int jel = 0; jel < numElemCoarse; jel++){  

        if (domDecompCoarse.first[jel] == rank) {  

            int *connec = elementsCoarse_[jel] -> getConnectivity();

            double **elemMatrix;
            elemMatrix = new double*[18]();
            for (int i = 0; i < 18; ++i)  elemMatrix[i] = new double[18]();
            double elemVector[18] = {};

            elementsCoarse_[jel] -> getTransientNavierStokes_FEM(elemMatrix,elemVector);

            for (int i=0; i<6; i++){

                    for (int j=0; j<6; j++){
                         
                        int dof_i = 2 * connec[i];
                        int dof_j = 2 * connec[j];
                        ierr = MatSetValues(A, 1, &dof_i,1, &dof_j,
                                            &elemMatrix[2*i  ][2*j  ],
                                            ADD_VALUES);

                        dof_i = 2 * connec[i] + 1;
                        dof_j = 2 * connec[j];
                        ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                            &elemMatrix[2*i+1][2*j  ],
                                            ADD_VALUES);

                        dof_i = 2 * connec[i];
                        dof_j = 2 * connec[j] + 1;
                        ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                            &elemMatrix[2*i  ][2*j+1],
                                            ADD_VALUES);

                        dof_i = 2 * connec[i] + 1;
                        dof_j = 2 * connec[j] + 1;
                        ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                            &elemMatrix[2*i+1][2*j+1],
                                            ADD_VALUES);

                        dof_i = 2 * connec[i];
                        dof_j = 2 * NCNumberNodesC + connec[j];
                        ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                            &elemMatrix[2*i  ][12+j],
                                            ADD_VALUES);

                        dof_i = 2 * NCNumberNodesC + connec[i];
                        dof_j = 2 * connec[j];
                        ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                            &elemMatrix[12+i][2*j  ],
                                            ADD_VALUES);

                        dof_i = 2 * connec[i] + 1;
                        dof_j = 2 * NCNumberNodesC + connec[j];
                        ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                            &elemMatrix[2*i+1][12+j],
                                            ADD_VALUES);

                        dof_i = 2 * NCNumberNodesC + connec[i];
                        dof_j = 2 * connec[j] + 1;
                        ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                            &elemMatrix[12+i][2*j+1],
                                            ADD_VALUES);

                        dof_i = 2 * NCNumberNodesC + connec[i];
                        dof_j = 2 * NCNumberNodesC + connec[j];
                        ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                            &elemMatrix[12+i][12+j],
                                            ADD_VALUES);
                    }; //loop j
                    
                    //Rhs vector
                    int dof_i = 2 * connec[i];
                    ierr = VecSetValues(b, 1, &dof_i, &elemVector[2*i  ],
                                        ADD_VALUES);
                
                    dof_i = 2 * connec[i] + 1;
                    ierr = VecSetValues(b, 1, &dof_i, &elemVector[2*i+1],
                                        ADD_VALUES);

                    dof_i = 2 * NCNumberNodesC + connec[i];
                    ierr = VecSetValues(b, 1, &dof_i, &elemVector[12+i],
                                        ADD_VALUES);
                 };// loop i
            
             for (int i = 0; i < 18; ++i) delete [] elemMatrix[i];
             delete [] elemMatrix;

        };
    };

};

template<>
void Arlequin<2>::setMatVecValuesCoarseISO(){


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

};

template<>
void Arlequin<2>::setMatVecValuesFineFEM(){
   

    for (int jel = 0; jel < numElemFine; jel++){   
        
        if (domDecompFine.first[jel] == rank) { 

    		int *connec = elementsFine_[jel] -> getConnectivity();
           
        	double **elemMatrix;
    		elemMatrix = new double*[18]();
    		for (int i = 0; i < 18; ++i)  elemMatrix[i] = new double[18]();
    		double elemVector[18] = {};

            elementsFine_[jel] -> getTransientNavierStokes_FEM(elemMatrix,elemVector);

            //Disperse local contributions into the global matrix
            //Matrix K and C
            for (int i=0; i<6; i++){	                        
                for (int j=0; j<6; j++){
                                                    
                    int dof_i = 2 * connec[i] + 3 * NCNumberNodesC;
                    int dof_j = 2 * connec[j] + 3 * NCNumberNodesC;
                    ierr = MatSetValues(A, 1, &dof_i,1, &dof_j,
                                        &elemMatrix[2*i  ][2*j  ],
                                        ADD_VALUES);

                    dof_i = 2 * connec[i] + 1 + 3 * NCNumberNodesC;
                    dof_j = 2 * connec[j] + 3 * NCNumberNodesC;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrix[2*i+1][2*j  ],
                                        ADD_VALUES);

                    dof_i = 2 * connec[i] + 3 * NCNumberNodesC;
                    dof_j = 2 * connec[j] + 1 + 3 * NCNumberNodesC;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrix[2*i  ][2*j+1],
                                        ADD_VALUES);

                    dof_i = 2 * connec[i] + 1 + 3 * NCNumberNodesC;
                    dof_j = 2 * connec[j] + 1 + 3 * NCNumberNodesC;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrix[2*i+1][2*j+1],
                                        ADD_VALUES);
            
                 //Matrix Q and Qt
                    dof_i = 2 * connec[i] + 3 * NCNumberNodesC;
                    dof_j = 2 * NCNumberNodesF + connec[j] + 3 * NCNumberNodesC;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrix[2*i  ][12+j],
                                        ADD_VALUES);

                    dof_i = 2 * NCNumberNodesF + connec[i] + 3 * NCNumberNodesC;
                    dof_j = 2 * connec[j] + 3 * NCNumberNodesC;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrix[12+i][2*j  ],
                                        ADD_VALUES);

                    dof_i = 2 * connec[i] + 1 + 3 * NCNumberNodesC;
                    dof_j = 2 * NCNumberNodesF + connec[j] + 3 * NCNumberNodesC;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrix[2*i+1][12+j],
                                        ADD_VALUES);

                    dof_i = 2 * NCNumberNodesF + connec[i] + 3 * NCNumberNodesC;
                    dof_j = 2 * connec[j] + 1 + 3 * NCNumberNodesC;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrix[12+i][2*j+1],
                                        ADD_VALUES);

                    dof_i = 2 * NCNumberNodesF + connec[i] + 3 * NCNumberNodesC;
                    dof_j = 2 * NCNumberNodesF + connec[j] + 3 * NCNumberNodesC;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemMatrix[12+i][12+j],
                                        ADD_VALUES);
                }; //loop j
                
                //Rhs vector
                int dof_i = 2 * connec[i] + 3 * NCNumberNodesC;
                ierr = VecSetValues(b, 1, &dof_i, &elemVector[2*i  ],
                                    ADD_VALUES);

                dof_i = 2 * connec[i] + 1 + 3 * NCNumberNodesC;
                ierr = VecSetValues(b, 1, &dof_i, &elemVector[2*i+1],
                                    ADD_VALUES);

                dof_i = 2 * NCNumberNodesF + connec[i] + 3 * NCNumberNodesC;
                ierr = VecSetValues(b, 1, &dof_i, &elemVector[12+i],
                                    ADD_VALUES);

             };// loop i
             for (int i = 0; i < 18; ++i) delete [] elemMatrix[i];
			 delete [] elemMatrix;  

        }; // domain decomposition
    }; //Elements Fine
};


template<>
void Arlequin<2>::setMatVecValuesFineISO(){
   

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
};

template<>
void Arlequin<2>::setMatVecValuesLagrangeFineFEM(int &iTimeStep){

	double &alpha_f = parametersFine ->getAlphaF();
	double &gamma = parametersFine ->getGamma();
    double &dTime = parametersFine -> getTimeStep();
    double integ = alpha_f * gamma * dTime;

    for (int l = 0; l < numElemGlueZoneFine; l++){

    	int jel =  elementsGlueZoneFine_[l];

    	if (domDecompFine.first[jel] == rank) {



            // //IF TARLQ COMPUTE BY ELEM
            // int nIP = elementsFine_[jel] -> getNumberOfIntegrationPointsSpecial_FEM();
            
            // double tarlq = 0.;
            // for (int ip = 0; ip< nIP; ip++){

            // 	int iElemCoarse = elementsFine_[jel] -> getIntegPointCorrespondenceElement_FEM(ip);
            // 	int *connecC = elementsCoarse_[iElemCoarse] -> getConnectivity();
            // 	int patch = elementsCoarse_[iElemCoarse] -> getPatch();
                
            //     double tarlqip;
            //     elementsFine_[jel] -> getParameterArlequinElem(tarlqip, ip, patch, nodesCoarse_,connecC,IsoParCoarse);

            //     tarlq += tarlqip;
            
            // };

            // elementsFine_[jel] -> setTarlq(tarlq);
            // // IF TARLQ COMPUTES BY ELEM



        	int *connec = elementsFine_[jel] -> getConnectivity();
        	int *connecL = glueZoneFine_[l] -> getConnectivity();	

        	//Sets global nodal equivalent velocities
            // for (int i = 0; i < 6; i++){
                
            //     int iElemCoarse = nodesFine_[connec[i]] -> getNodalElemCorrespondence();
            //     int *connecC = elementsCoarse_[iElemCoarse] -> getConnectivity();

            //     double u1[9],u2[9];
            //     for (int j = 0; j < 9; j++){
            //         u1[j] = nodesCoarse_[connecC[j]] -> getVelocity(0);
            //         u2[j] = nodesCoarse_[connecC[j]] -> getVelocity(1);
            //     };

            //     //Computing velocity in the coarse mesh
            //     double *xsi = nodesFine_[connec[i]] -> getNodalXsiCorrespondence();
            //     QuadShapeFunction<2>  shapeQuad;
            //     double phi_[9],wpc[9];
            //     int patch = elementsCoarse_[iElemCoarse] -> getPatch();
            //     int *inc = nodesCoarse_[connecC[8]] -> getINC();
            //     for (int i = 0; i < 9; i++) wpc[i] = nodesCoarse_[connecC[i]] -> getWeightPC();
            //     shapeQuad.evaluateIso(xsi,phi_,wpc,inc,IsoParCoarse,patch);
                
            //     double u_[2] = {};
            //     for (int j = 0; j < 9; j++){
            //         u_[0] += u1[j] * phi_[j];
            //         u_[1] += u2[j] * phi_[j];
            //     };

            //     nodesFine_[connec[i]] -> setGlobalVelocity(u_);
            // };


            int numberIntPoints = elementsFine_[jel] -> getNumberOfIntegrationPointsSpecial_FEM();


            for (int ip = 0; ip< numberIntPoints; ip++){

            	int iElemCoarse = elementsFine_[jel] -> getIntegPointCorrespondenceElement_FEM(ip);
            	int *connecC = elementsCoarse_[iElemCoarse] -> getConnectivity();
            	int patch = elementsCoarse_[iElemCoarse] -> getPatch();

	            //LAGRANGE MULTIPLIERS MATRIXES AND VECTORS
                double **elemMatrixLag1;
    		    elemMatrixLag1 = new double*[12]();
    		    for (int i = 0; i < 12; ++i)  elemMatrixLag1[i] = new double[18]();
    		    
                double elemVectorLag1_1[18] = {};
                double elemVectorLag1_2[12] = {};


	            //tSUPG and tPSPG STABILIZATION
                double **jacobianNRMatrix;
    		    jacobianNRMatrix = new double*[18]();
    		    for (int i = 0; i < 18; ++i)  jacobianNRMatrix[i] = new double[12]();
                double rhsVector[18] = {};

	    		//ARLEQUIN STABILIZATION MATRIXES
                 double **elemStabMatrixD;
    		    elemStabMatrixD = new double*[12]();
    		    for (int i = 0; i < 12; ++i)  elemStabMatrixD[i] = new double[12]();

                 double **elemStabMatrix1;
    		    elemStabMatrix1 = new double*[12]();
    		    for (int i = 0; i < 12; ++i)  elemStabMatrix1[i] = new double[18]();

	            double elemStabVectorD[12] = {};
                double elemStabVector1[12] = {};
	            
	    		elementsFine_[jel] -> setTimeStep(iTimeStep);
	    		elementsFine_[jel] -> getLagrangeMultipliersSameMesh_FEM(ip,elemMatrixLag1,elemVectorLag1_1,elemVectorLag1_2);
	            // elementsFine_[jel] -> getLagrangeMultipliersSameMesh_tSUPG_tPSPG_FEM(ip,jacobianNRMatrix,rhsVector);
	            elementsFine_[jel] -> getLagrangeMultipliersSameMeshArlqStab_FEM(ip,patch,nodesCoarse_,connecC,IsoParCoarse,
                                                                                 elemStabMatrixD,elemStabVectorD,
	    																	     elemStabMatrix1,elemStabVector1);

			
	    		for (int i = 0; i < 6; i++){
	    			for (int j = 0; j < 6; j++){

	    				//Lagrange Multipliers
	    				int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	    				int dof_j = 3*NCNumberNodesC + 2*connec[j];
	    				double value = integ * elemMatrixLag1[2*i][2*j];
	    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                        &value,ADD_VALUES);
	    				ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
	                                        &elemMatrixLag1[2*i][2*j],
	                                        ADD_VALUES);

	    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	    				dof_j = 3*NCNumberNodesC + 2*connec[j];
	    				value = integ * elemMatrixLag1[2*i+1][2*j];
	    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                        &value,ADD_VALUES);
	    				ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
	                                        &elemMatrixLag1[2*i+1][2*j],
	                                        ADD_VALUES);

	    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	    				dof_j = 3*NCNumberNodesC + 2*connec[j] + 1;
	    				value = integ * elemMatrixLag1[2*i][2*j+1];
	    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                        &value,ADD_VALUES);
	    				ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
	                                        &elemMatrixLag1[2*i][2*j+1],
	                                        ADD_VALUES);

	    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	    				dof_j = 3*NCNumberNodesC + 2*connec[j] + 1;
	    				value = integ * elemMatrixLag1[2*i+1][2*j+1];
	    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                        &value,ADD_VALUES);
	    				ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
	                                        &elemMatrixLag1[2*i+1][2*j+1],
	                                        ADD_VALUES);

	                    //tSUPG 
	                    dof_i = 3*NCNumberNodesC + 2*connec[i];
	                    dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
	                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                        &jacobianNRMatrix[2*i][2*j],ADD_VALUES);
	 
	                    dof_i = 3*NCNumberNodesC + 2*connec[i] + 1;
	                    dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1;
	                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                        &jacobianNRMatrix[2*i+1][2*j+1],ADD_VALUES);

	                    //tPSPG
	                    dof_i = 3*NCNumberNodesC + 2*NCNumberNodesF + connec[i];
	                    dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
	                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                        &jacobianNRMatrix[12+i][2*j],ADD_VALUES);

	                    dof_i = 3*NCNumberNodesC + 2*NCNumberNodesF + connec[i];
	                    dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1;
	                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                        &jacobianNRMatrix[12+i][2*j+1],ADD_VALUES);

	                    //Arlequin Stabilization from diagonal matrix
	    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	    				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
	    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                        &elemStabMatrixD[2*i][2*j],
	                                        ADD_VALUES);

	    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	    				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1;
	    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                        &elemStabMatrixD[2*i+1][2*j+1],
	                                        ADD_VALUES);

	    				//Stabilization Arlequin terms from fine mesh matrix
	    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	    				dof_j = 3*NCNumberNodesC + 2*connec[j];
	    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                        &elemStabMatrix1[2*i][2*j],
	                                        ADD_VALUES);

	    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i]+1;
	    				dof_j = 3*NCNumberNodesC + 2*connec[j] + 1;
	    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                        &elemStabMatrix1[2*i+1][2*j+1],
	                                        ADD_VALUES);

	    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	    				dof_j = 3*NCNumberNodesC + 2*NCNumberNodesF + connec[j];
	    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                        &elemStabMatrix1[2*i][12+j],
	                                        ADD_VALUES);

	    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	    				dof_j = 3*NCNumberNodesC + 2*NCNumberNodesF + connec[j];
	    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
	                                        &elemStabMatrix1[2*i+1][12+j],
	                                        ADD_VALUES);

	    			};//j

	    			//Lagrange Multipliers
	    			int dof_i = 3*NCNumberNodesC + 2*connec[i];
	    			ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1_1[2*i],ADD_VALUES);
	    			dof_i = 3*NCNumberNodesC + 2*connec[i] + 1;
	    			ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1_1[2*i+1],ADD_VALUES);

	    			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	    			ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1_2[2*i  ],ADD_VALUES);
	    			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	    			ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1_2[2*i+1 ],ADD_VALUES);

	                //tSUPG 
	                dof_i = 3*NCNumberNodesC + 2*connec[i];
	                ierr = VecSetValues(b, 1, &dof_i, &rhsVector[2*i],ADD_VALUES);
	                dof_i = 3*NCNumberNodesC + 2*connec[i] + 1;
	                ierr = VecSetValues(b, 1, &dof_i, &rhsVector[2*i+1],ADD_VALUES);

	                //tPSPG
	                dof_i = 3*NCNumberNodesC + 2*NCNumberNodesF + connec[i];
	                ierr = VecSetValues(b, 1, &dof_i, &rhsVector[12+i],ADD_VALUES);

	    			//Arlequin stabilization term from diagonal
	    			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	    			ierr = VecSetValues(b, 1, &dof_i, &elemStabVectorD[2*i  ],ADD_VALUES);
	    			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	    			ierr = VecSetValues(b, 1, &dof_i, &elemStabVectorD[2*i+1 ],ADD_VALUES);

	    			//Arlequin stabilization term from fine mesh
	    			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
	    			ierr = VecSetValues(b, 1, &dof_i, &elemStabVector1[2*i  ],ADD_VALUES);
	    			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
	    			ierr = VecSetValues(b, 1, &dof_i, &elemStabVector1[2*i+1 ],ADD_VALUES);
	    		};//i

                for (int i = 0; i < 12; ++i) {
                    delete [] elemMatrixLag1[i];
                    delete [] elemStabMatrixD[i];
                    delete [] elemStabMatrix1[i];
                }
                delete [] elemMatrixLag1; 
                delete [] elemStabMatrixD;
                delete [] elemStabMatrix1;

                for (int i = 0; i < 18; ++i) {
                    delete [] jacobianNRMatrix[i];
                }
                delete [] jacobianNRMatrix; 

            }//integration points		

    	};//decomposition

    };//gluezonefine

};

template<>
void Arlequin<2>::setMatVecValuesLagrangeFineISO(){

	double &alpha_f = parametersFine -> getAlphaF();
	double &gamma = parametersFine -> getGamma();
    double& dTime = parametersFine -> getTimeStep();

    double integ = alpha_f * gamma * dTime;


    for (int l = 0; l < numElemGlueZoneFine; l++){

    	int jel =  elementsGlueZoneFine_[l];

    	if (domDecompFine.first[jel] == rank) {

        	int *connec = elementsFine_[jel] -> getConnectivity();
        	int *connecL = glueZoneFine_[l] -> getNewConnectivity();	

    		//LAGRANGE MULTIPLIERS MATRIXES AND VECTORS
    		double **elemMatrixLag1;
    		elemMatrixLag1 = new double*[18]();
    		for (int i = 0; i < 18; ++i)  elemMatrixLag1[i] = new double[27]();
    		double elemVectorLag1_1[27] = {};
    		double elemVectorLag1_2[18] = {};

            //tSUPG and tPSPG STABILIZATION
            double **jacobianNRMatrix;
            jacobianNRMatrix = new double*[27]();
            for (int i = 0; i < 27; ++i)  jacobianNRMatrix[i] = new double[18]();
            double rhsVector[27] = {};

            //ARLEQUIN STABILIZATION MATRIX
            double **elemStabMatrixD;
            elemStabMatrixD = new double*[18]();
            for (int i = 0; i < 18; ++i)  elemStabMatrixD[i] = new double[18]();
            double elemStabVectorD[18] = {};

            double **elemStabMatrix1;
            elemStabMatrix1 = new double*[18]();
            for (int i = 0; i < 18; ++i)  elemStabMatrix1[i] = new double[27]();
            double elemStabVector1[18] = {};
    		
    		
            elementsFine_[jel] -> getLagrangeMultipliersSameMesh_ISO(elemMatrixLag1,elemVectorLag1_1,elemVectorLag1_2);
            elementsFine_[jel] -> getLagrangeMultipliersSameMesh_tSUPG_tPSPG_ISO(jacobianNRMatrix,rhsVector);
            elementsFine_[jel] -> getLagrangeMultipliersSameMeshArlqStab_ISO(elemStabMatrixD,elemStabVectorD,
                                                                             elemStabMatrix1,elemStabVector1);

    		for (int i = 0; i < 9; i++){

                int newconi = nodesFine_[connec[i]] -> getnewcon();

    			for (int j = 0; j < 9; j++){

    				int newconj = nodesFine_[connec[j]] -> getnewcon();
                    
    				//Lagrange Multipliers
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

                    //tSUPG 
                    dof_i = 3*NCNumberNodesC + 2*newconi;
                    dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &jacobianNRMatrix[2*i][2*j],ADD_VALUES);
 
                    dof_i = 3*NCNumberNodesC + 2*newconi + 1;
                    dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &jacobianNRMatrix[2*i+1][2*j+1],ADD_VALUES);

                    //tPSPG
                    dof_i = 3*NCNumberNodesC + 2*NCNumberNodesF + newconi;
                    dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &jacobianNRMatrix[18+i][2*j],ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 2*NCNumberNodesF + newconi;
                    dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &jacobianNRMatrix[18+i][2*j+1],ADD_VALUES);

                    //Stabilization Arlequin terms from fine mesh matrix
                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    dof_j = 3*NCNumberNodesC + 2*newconj;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix1[2*i][2*j],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i]+1;
                    dof_j = 3*NCNumberNodesC + 2*newconj + 1;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix1[2*i+1][2*j+1],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    dof_j = 3*NCNumberNodesC + 2*NCNumberNodesF + newconj;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix1[2*i][18+j],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
                    dof_j = 3*NCNumberNodesC + 2*NCNumberNodesF + newconj;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix1[2*i+1][18+j],
                                        ADD_VALUES);

    				//Stabilization Arlequin Terms diagonal
    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
    				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrixD[2*i][2*j],
                                        ADD_VALUES);

    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
    				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1;
    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrixD[2*i+1][2*j+1],
                                        ADD_VALUES);

    			};//j

    			
                //Lagrange Multipliers
    			int dof_i = 3*NCNumberNodesC + 2*newconi;
    			ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1_1[2*i],ADD_VALUES);
    			dof_i = 3*NCNumberNodesC + 2*newconi + 1;
    			ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1_1[2*i+1],ADD_VALUES);

    			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
    			ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1_2[2*i  ],ADD_VALUES);
    			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
    			ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag1_2[2*i+1 ],ADD_VALUES);

    			//Stabilization Arlequin Term
    			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
    			ierr = VecSetValues(b, 1, &dof_i, &elemStabVectorD[2*i  ],ADD_VALUES);
    			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
    			ierr = VecSetValues(b, 1, &dof_i, &elemStabVectorD[2*i+1 ],ADD_VALUES);

                //tSUPG 
                dof_i = 3*NCNumberNodesC + 2*newconi;
                ierr = VecSetValues(b, 1, &dof_i, &rhsVector[2*i],ADD_VALUES);
                dof_i = 3*NCNumberNodesC + 2*newconi + 1;
                ierr = VecSetValues(b, 1, &dof_i, &rhsVector[2*i+1],ADD_VALUES);

                //tPSPG
                dof_i = 3*NCNumberNodesC + 2*NCNumberNodesF + newconi;
                ierr = VecSetValues(b, 1, &dof_i, &rhsVector[18+i],ADD_VALUES);

                //Arlequin stabilization term from fine mesh
                dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                ierr = VecSetValues(b, 1, &dof_i, &elemStabVector1[2*i  ],ADD_VALUES);
                dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
                ierr = VecSetValues(b, 1, &dof_i, &elemStabVector1[2*i+1 ],ADD_VALUES);


    		};//i


    		for (int i = 0; i < 18; ++i) {
            	delete [] elemMatrixLag1[i];
            	delete [] elemStabMatrixD[i];
                delete [] elemStabMatrix1[i];
            }
    		delete [] elemMatrixLag1; 
    		delete [] elemStabMatrixD;
            delete [] elemStabMatrix1;

            for (int i = 0; i < 27; ++i) {
                delete [] jacobianNRMatrix[i];
            }
            delete [] jacobianNRMatrix; 

    	};//decomposition

    };//gluezonefine

};

template<>
void Arlequin<2>::setMatVecValuesLagrangeCoarseFEM_FEM(){


	double &alpha_f = parametersFine -> getAlphaF();
	double &gamma = parametersFine -> getGamma();
    double& dTime = parametersFine -> getTimeStep();
    
    double integ = alpha_f * gamma * dTime;

     //numElemGlueZoneFine
    for (int l = 0; l < numElemGlueZoneFine; l++){

    	int jel =  elementsGlueZoneFine_[l];

    	if (domDecompFine.first[jel] == rank) {

        	int *connecL = glueZoneFine_[l] -> getConnectivity();	

            int numberIntPoints = elementsFine_[jel] -> getNumberOfIntegrationPointsSpecial_FEM();

            std::vector<int> ele, diffElem;
            ele.clear();
            diffElem.clear();

            //Finding coarse elements 
            for (int i=0; i<numberIntPoints; i++){
                int aux = elementsFine_[jel] -> getIntegPointCorrespondenceElement_FEM(i);
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


        		//Lagrange multipliers matrixes and vectors
                double **elemMatrixLag0;
                elemMatrixLag0 = new double*[12]();
                for (int i = 0; i < 12; ++i)  elemMatrixLag0[i] = new double[18]();
                double elemVectorLag0_1[18] = {};
                double elemVectorLag0_2[12] = {};

                //Stabilization matrix and vector
                double **elemStabMatrix;
            	elemStabMatrix = new double*[12]();
            	for (int i = 0; i < 12; ++i)  elemStabMatrix[i] = new double[12]();
            	double elemStabVector[12] = {};

                                     
                int iElemCoarse = diffElem[ielem];

                int *connecC = elementsCoarse_[iElemCoarse] -> getConnectivity();
                
                elementsFine_[jel] -> getLagrangeMultipliersDifferentMesh_FEM_FEM(nodesCoarse_,connecC,iElemCoarse, 
                                                                                  elemMatrixLag0,elemVectorLag0_1,elemVectorLag0_2,
                                                                                  elemStabMatrix,elemStabVector);

                for (int i = 0; i < 6; i++){
                    for (int j = 0; j < 6; j++){

                    int newconj = connecC[j];  
                    
                    int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    int dof_j = 2*connecC[j];
                    double value = integ * elemMatrixLag0[2*i][2*j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &value,ADD_VALUES);
                    ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag0[2*i][2*j],
                                        ADD_VALUES);


                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i]+1;
                    dof_j = 2*connecC[j];
                    value = integ * elemMatrixLag0[2*i+1][2*j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &value,ADD_VALUES);
                    ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag0[2*i+1][2*j],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    dof_j = 2*connecC[j] +1 ;
                    value = integ * elemMatrixLag0[2*i][2*j+1];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &value,ADD_VALUES);
                    ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag0[2*i][2*j+1],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1 ;
                    dof_j = 2*connecC[j] +1 ;
                    value = integ * elemMatrixLag0[2*i+1][2*j+1];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &value,ADD_VALUES);
                    ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag0[2*i+1][2*j+1],
                                        ADD_VALUES);

                    //Stabilization Arlequin Terms diagonal
    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
    				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix[2*i][2*j],
                                        ADD_VALUES);

    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
    				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix[2*i+1][2*j],
                                        ADD_VALUES);

    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
    				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1;
    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix[2*i][2*j+1],
                                        ADD_VALUES);


    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
    				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1;
    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix[2*i+1][2*j+1],
                                        ADD_VALUES);

                    };//j

                    int dof_i = 2*connecC[i];
                    ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_1[2*i  ],ADD_VALUES);
                    dof_i = 2*connecC[i] + 1 ;
                    ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_1[2*i+1],ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_2[2*i  ],ADD_VALUES);
                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
                    ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_2[2*i+1],ADD_VALUES);


                    //Stabilization Arlequin Term
        			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
        			ierr = VecSetValues(b, 1, &dof_i, &elemStabVector[2*i  ],ADD_VALUES);

        			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
        			ierr = VecSetValues(b, 1, &dof_i, &elemStabVector[2*i+1 ],ADD_VALUES);                    

                };//j
                  
                
                for (int i = 0; i < 12; ++i) {
                delete [] elemMatrixLag0[i];
                delete [] elemStabMatrix[i];
                 }
                delete [] elemMatrixLag0; 
                delete [] elemStabMatrix;

            };//intersect
  
    	};//decomposition

    };//gluezonefine


};

// template<>
// void Arlequin<2>::setMatVecValuesLagrangeCoarseFEM_ISO(int &iTimeStep){

// 	double &alpha_f = parametersFine -> getAlphaF();
// 	double &gamma = parametersFine -> getGamma();
//     double& dTime = parametersFine -> getTimeStep();

//     double integ = alpha_f * gamma * dTime;

//     //Coarse mesh
//     for (int l = 0; l < numElemGlueZoneFine; l++){

//     	int jel =  elementsGlueZoneFine_[l];

//     	if (domDecompFine.first[jel] == rank) {

//         	int *connecL = glueZoneFine_[l] -> getConnectivity();	

//             int numberIntPoints = elementsFine_[jel] -> getNumberOfIntegrationPointsSpecial_FEM();

//             std::vector<int> ele, diffElem;
//             ele.clear();
//             diffElem.clear();

//             //Finding coarse elements 
//             for (int i=0; i<numberIntPoints; i++){
//                 int aux = elementsFine_[jel] -> getIntegPointCorrespondenceElement_FEM(i);
//                 ele.push_back(aux);
//             };

//             int numElemIntersect = 1;
//             int flag = 0;
//             diffElem.push_back(ele[0]);
        
//             for (int i = 1; i<numberIntPoints; i++){
//                 flag = 0;
//                 for (int j = 0; j<numElemIntersect; j++){
//                     if (ele[i] == diffElem[j]) {
//                         break;
//                     }else{
//                         flag++;
//                     };
//                     if(flag == numElemIntersect){
//                         numElemIntersect++;
//                         diffElem.push_back(ele[i]);
//                     };
//                 };
//             };


//             for (int ielem = 0; ielem < numElemIntersect; ielem++){

                 
//             	//LAGRANGE MULTIPLIERS MATRIXES AND VECTORS
//                 double **elemMatrixLag0;
//                 elemMatrixLag0 = new double*[12]();
//                 for (int i = 0; i < 12; ++i)  elemMatrixLag0[i] = new double[27]();
//                 double elemVectorLag0_1[27] = {};
//                 double elemVectorLag0_2[12] = {};

//                 //tSUPG and tPSPG STABILIZATION
//             	double **jacobianNRMatrix;
//             	jacobianNRMatrix = new double*[27]();
//             	for (int i = 0; i < 27; ++i)  jacobianNRMatrix[i] = new double[12]();
//             	double rhsVector[27] = {};

//                 //ARLEQUIN STABILIZATION MATRIX AND VECTOR
//                 double **elemStabMatrixD;
//                 elemStabMatrixD = new double*[12]();
//                 for (int i = 0; i < 12; ++i)  elemStabMatrixD[i] = new double[12]();
//                 double elemStabVectorD[12] = {};

//                 double **elemStabMatrix0;
//                 elemStabMatrix0 = new double*[12]();
//                 for (int i = 0; i < 12; ++i)  elemStabMatrix0[i] = new double[27]();
//                 double elemStabVector0[12] = {};

                                     
//                 int iElemCoarse = diffElem[ielem];

//                 int *connecC = elementsCoarse_[iElemCoarse] -> getConnectivity();
//                 int patch = elementsCoarse_[iElemCoarse] -> getPatch();
                
//                 elementsFine_[jel] ->setTimeStep(iTimeStep);

//                 elementsFine_[jel] -> getLagrangeMultipliersDifferentMesh_FEM_ISO(patch,nodesCoarse_,connecC,IsoParCoarse,iElemCoarse, 
//                                                                                   elemMatrixLag0,elemVectorLag0_1,elemVectorLag0_2);
//                 // elementsFine_[jel] -> getLagrangeMultipliersDifferentMesh_tSUPG_tPSPG_FEM_ISO(patch,nodesCoarse_,connecC,IsoParCoarse,
//                 //  																			  iElemCoarse, jacobianNRMatrix,rhsVector);
//                 // elementsFine_[jel] -> getLagrangeMultipliersDifferentMeshArlqStab_FEM_ISO(patch,nodesCoarse_,connecC,IsoParCoarse,iElemCoarse, 
//                 //                                                                           elemStabMatrixD,elemStabVectorD,elemStabMatrix0,
//                 //                                                                           elemStabVector0);
      
//                 for (int i = 0; i < 6; i++){
//                     for (int j = 0; j < 9; j++){

//                     int newconj = nodesCoarse_[connecC[j]] -> getnewcon();  
                    
//                     //Lagrange multipliers matrixes
//                     int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
//                     int dof_j = 2*newconj;
//                     double value = integ * elemMatrixLag0[2*i][2*j];
//                     ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
//                                         &value,ADD_VALUES);
//                     ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
//                                         &elemMatrixLag0[2*i][2*j],
//                                         ADD_VALUES);

//                     dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i]+1;
//                     dof_j = 2*newconj;
//                     value = integ * elemMatrixLag0[2*i+1][2*j];
//                     ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
//                                         &value,ADD_VALUES);
//                     ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
//                                         &elemMatrixLag0[2*i+1][2*j],
//                                         ADD_VALUES);

//                     dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
//                     dof_j = 2*newconj +1 ;
//                     value = integ * elemMatrixLag0[2*i][2*j+1];
//                     ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
//                                         &value,ADD_VALUES);
//                     ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
//                                         &elemMatrixLag0[2*i][2*j+1],
//                                         ADD_VALUES);

//                     dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1 ;
//                     dof_j = 2*newconj +1 ;
//                     value = integ * elemMatrixLag0[2*i+1][2*j+1];
//                     ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
//                                         &value,ADD_VALUES);
//                     ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
//                                         &elemMatrixLag0[2*i+1][2*j+1],
//                                         ADD_VALUES);

//                     //Arlequin Stabilization terms coarse mesh
//                     dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
//                     dof_j = 2*newconj;
//                     ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
//                                         &elemStabMatrix0[2*i][2*j],
//                                         ADD_VALUES);

//                     dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
//                     dof_j = 2*newconj + 1;
//                     ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
//                                         &elemStabMatrix0[2*i+1][2*j+1],
//                                         ADD_VALUES);

//                     dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
//                     dof_j = 2*NCNumberNodesC + newconj;
//                     ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
//                                         &elemStabMatrix0[2*i][18+j],
//                                         ADD_VALUES);

//                     dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
//                     dof_j = 2*NCNumberNodesC + newconj;
//                     ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
//                                         &elemStabMatrix0[2*i+1][18+j],
//                                         ADD_VALUES);

//                     };//j

//                     //Lagrange multipliers vector
//                     int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
//                     ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_2[2*i  ],ADD_VALUES);

//                     dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
//                     ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_2[2*i+1],ADD_VALUES);

//                     //Arlequin Stabilization terms coarse mesh
//                     dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
//                     ierr = VecSetValues(b, 1, &dof_i, &elemStabVector0[2*i  ],ADD_VALUES);

//                     dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
//                     ierr = VecSetValues(b, 1, &dof_i, &elemStabVector0[2*i+1 ],ADD_VALUES);  


//                 };//j


//                 for (int i = 0; i < 9; i++){
      				
//                     int newconi = nodesCoarse_[connecC[i]] -> getnewcon(); 
                    
//                     for (int j = 0; j < 6; j++){
                    
//                     //tSUPG stabilization
//                     int dof_i = 2*newconi;
//                     int dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
//                     ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
//                                         &jacobianNRMatrix[2*i][2*j],
//                                         ADD_VALUES);

//                     dof_i = 2*newconi + 1 ;
//                     dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1 ;
//                     ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
//                                         &jacobianNRMatrix[2*i+1][2*j+1],
//                                         ADD_VALUES);	

//                     //tPSPG stabilization
//                     dof_i = 2*NCNumberNodesC + newconi;
//                     dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
//                     ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
//                                         &jacobianNRMatrix[18+i][2*j],
//                                         ADD_VALUES);

//                     dof_i = 2*NCNumberNodesC + newconi;
//                     dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1;
//                     ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
//                                         &jacobianNRMatrix[18+i][2*j+1],
//                                         ADD_VALUES);


//                     };//j

//                     //Lagrange multipliers vector
//                     int dof_i = 2*newconi;
//                     ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_1[2*i  ],ADD_VALUES);
//                     dof_i = 2*newconi + 1 ;
//                     ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_1[2*i+1],ADD_VALUES);

//                     //tSUPG stabilization
//                     dof_i = 2*newconi;
//                     ierr = VecSetValues(b, 1, &dof_i, &rhsVector[2*i  ],ADD_VALUES);
//                     dof_i = 2*newconi + 1 ;
//                     ierr = VecSetValues(b, 1, &dof_i, &rhsVector[2*i+1],ADD_VALUES);

//                     //tPSPG stabilization
//                     dof_i = 2*NCNumberNodesC + newconi;
//                     ierr = VecSetValues(b, 1, &dof_i, &rhsVector[18+i],ADD_VALUES);
                
//                 };//j

    
//                 //Lagrange stabilization terms (Lagrange multipliers defined in the fine mesh)
//                 for (int i = 0; i < 6; i++){
//                     for (int j = 0; j < 6; j++){

//                     int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
//                     int dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
//                     ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
//                                         &elemStabMatrixD[2*i][2*j],
//                                         ADD_VALUES);

//                     dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
//                     dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1;
//                     ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
//                                         &elemStabMatrixD[2*i+1][2*j+1],
//                                         ADD_VALUES);

//                     };

//                     int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
//                     ierr = VecSetValues(b, 1, &dof_i, &elemStabVectorD[2*i  ],ADD_VALUES);

//                     dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
//                     ierr = VecSetValues(b, 1, &dof_i, &elemStabVectorD[2*i+1 ],ADD_VALUES);  

//                 };                  
                
//                 for (int i = 0; i < 12; ++i) {
//                 delete [] elemMatrixLag0[i];
//                 delete [] elemStabMatrixD[i];
//                 delete [] elemStabMatrix0[i];
//                  }
//                 delete [] elemMatrixLag0; 
//                 delete [] elemStabMatrixD;
//                 delete [] elemStabMatrix0;

//                 for (int i = 0; i < 27; ++i) delete [] jacobianNRMatrix[i];
//                 delete [] jacobianNRMatrix;

//             };//intersect
  
//     	};//decomposition

//     };//gluezonefine
// };

template<>
void Arlequin<2>::setMatVecValuesLagrangeCoarseFEM_ISO(int &iTimeStep){

	double &alpha_f = parametersFine -> getAlphaF();
	double &gamma = parametersFine -> getGamma();
    double& dTime = parametersFine -> getTimeStep();

    double integ = alpha_f * gamma * dTime;

    //Coarse mesh
    for (int l = 0; l < numElemGlueZoneFine; l++){

    	int jel =  elementsGlueZoneFine_[l];

    	if (domDecompFine.first[jel] == rank) {

        	int *connecL = glueZoneFine_[l] -> getConnectivity();	

            int numberIntPoints = elementsFine_[jel] -> getNumberOfIntegrationPointsSpecial_FEM();

            std::vector<int> ele, diffElem;
            ele.clear();
            diffElem.clear();

            //Finding coarse elements 
            for (int i=0; i<numberIntPoints; i++){
                int aux = elementsFine_[jel] -> getIntegPointCorrespondenceElement_FEM(i);
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

                 
            	//LAGRANGE MULTIPLIERS MATRIXES AND VECTORS
                double **elemMatrixLag0;
                elemMatrixLag0 = new double*[12]();
                for (int i = 0; i < 12; ++i)  elemMatrixLag0[i] = new double[27]();
                double elemVectorLag0_1[27] = {};
                double elemVectorLag0_2[12] = {};

                //tSUPG and tPSPG STABILIZATION
            	double **jacobianNRMatrix;
            	jacobianNRMatrix = new double*[27]();
            	for (int i = 0; i < 27; ++i)  jacobianNRMatrix[i] = new double[12]();
            	double rhsVector[27] = {};

                //ARLEQUIN STABILIZATION MATRIX AND VECTOR
                double **elemStabMatrixD;
                elemStabMatrixD = new double*[12]();
                for (int i = 0; i < 12; ++i)  elemStabMatrixD[i] = new double[12]();
                double elemStabVectorD[12] = {};

                double **elemStabMatrix0;
                elemStabMatrix0 = new double*[12]();
                for (int i = 0; i < 12; ++i)  elemStabMatrix0[i] = new double[27]();
                double elemStabVector0[12] = {};

                                     
                int iElemCoarse = diffElem[ielem];

                int *connecC = elementsCoarse_[iElemCoarse] -> getConnectivity();
                int patch = elementsCoarse_[iElemCoarse] -> getPatch();
                
                elementsFine_[jel] ->setTimeStep(iTimeStep);

                elementsFine_[jel] -> getLagrangeMultipliersDifferentMesh_FEM_ISO(patch,nodesCoarse_,connecC,IsoParCoarse,iElemCoarse, 
                                                                                  elemMatrixLag0,elemVectorLag0_1,elemVectorLag0_2);
                // elementsFine_[jel] -> getLagrangeMultipliersDifferentMesh_tSUPG_tPSPG_FEM_ISO(patch,nodesCoarse_,connecC,IsoParCoarse,
                //  																			  iElemCoarse, jacobianNRMatrix,rhsVector);
                elementsFine_[jel] -> getLagrangeMultipliersDifferentMeshArlqStab_FEM_ISO(patch,nodesCoarse_,connecC,IsoParCoarse,iElemCoarse, 
                                                                                          elemStabMatrixD,elemStabVectorD,elemStabMatrix0,
                                                                                          elemStabVector0);
      
                for (int i = 0; i < 6; i++){
                    for (int j = 0; j < 9; j++){

                    int newconj = nodesCoarse_[connecC[j]] -> getnewcon();  
                    
                    //Lagrange multipliers matrixes
                    int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    int dof_j = 2*newconj;
                    double value = integ * elemMatrixLag0[2*i][2*j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &value,ADD_VALUES);
                    ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag0[2*i][2*j],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i]+1;
                    dof_j = 2*newconj;
                    value = integ * elemMatrixLag0[2*i+1][2*j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &value,ADD_VALUES);
                    ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag0[2*i+1][2*j],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    dof_j = 2*newconj +1 ;
                    value = integ * elemMatrixLag0[2*i][2*j+1];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &value,ADD_VALUES);
                    ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag0[2*i][2*j+1],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1 ;
                    dof_j = 2*newconj +1 ;
                    value = integ * elemMatrixLag0[2*i+1][2*j+1];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &value,ADD_VALUES);
                    ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag0[2*i+1][2*j+1],
                                        ADD_VALUES);

                    //Arlequin Stabilization terms coarse mesh
                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    dof_j = 2*newconj;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix0[2*i][2*j],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
                    dof_j = 2*newconj + 1;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix0[2*i+1][2*j+1],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    dof_j = 2*NCNumberNodesC + newconj;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix0[2*i][18+j],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
                    dof_j = 2*NCNumberNodesC + newconj;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix0[2*i+1][18+j],
                                        ADD_VALUES);

                    };//j

                    //Lagrange multipliers vector
                    int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_2[2*i  ],ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
                    ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_2[2*i+1],ADD_VALUES);

                    //Arlequin Stabilization terms coarse mesh
                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    ierr = VecSetValues(b, 1, &dof_i, &elemStabVector0[2*i  ],ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
                    ierr = VecSetValues(b, 1, &dof_i, &elemStabVector0[2*i+1 ],ADD_VALUES);  


                };//j


                for (int i = 0; i < 9; i++){
      				
                    int newconi = nodesCoarse_[connecC[i]] -> getnewcon(); 
                    
                    for (int j = 0; j < 6; j++){
                    
                    //tSUPG stabilization
                    int dof_i = 2*newconi;
                    int dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &jacobianNRMatrix[2*i][2*j],
                                        ADD_VALUES);

                    dof_i = 2*newconi + 1 ;
                    dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1 ;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &jacobianNRMatrix[2*i+1][2*j+1],
                                        ADD_VALUES);	

                    //tPSPG stabilization
                    dof_i = 2*NCNumberNodesC + newconi;
                    dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &jacobianNRMatrix[18+i][2*j],
                                        ADD_VALUES);

                    dof_i = 2*NCNumberNodesC + newconi;
                    dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &jacobianNRMatrix[18+i][2*j+1],
                                        ADD_VALUES);


                    };//j

                    //Lagrange multipliers vector
                    int dof_i = 2*newconi;
                    ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_1[2*i  ],ADD_VALUES);
                    dof_i = 2*newconi + 1 ;
                    ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_1[2*i+1],ADD_VALUES);

                    //tSUPG stabilization
                    dof_i = 2*newconi;
                    ierr = VecSetValues(b, 1, &dof_i, &rhsVector[2*i  ],ADD_VALUES);
                    dof_i = 2*newconi + 1 ;
                    ierr = VecSetValues(b, 1, &dof_i, &rhsVector[2*i+1],ADD_VALUES);

                    //tPSPG stabilization
                    dof_i = 2*NCNumberNodesC + newconi;
                    ierr = VecSetValues(b, 1, &dof_i, &rhsVector[18+i],ADD_VALUES);
                
                };//j

    
                //Lagrange stabilization terms (Lagrange multipliers defined in the fine mesh)
                for (int i = 0; i < 6; i++){
                    for (int j = 0; j < 6; j++){

                    int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    int dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrixD[2*i][2*j],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
                    dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrixD[2*i+1][2*j+1],
                                        ADD_VALUES);

                    };

                    int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    ierr = VecSetValues(b, 1, &dof_i, &elemStabVectorD[2*i  ],ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
                    ierr = VecSetValues(b, 1, &dof_i, &elemStabVectorD[2*i+1 ],ADD_VALUES);  

                };                  
                
                for (int i = 0; i < 12; ++i) {
                delete [] elemMatrixLag0[i];
                delete [] elemStabMatrixD[i];
                delete [] elemStabMatrix0[i];
                 }
                delete [] elemMatrixLag0; 
                delete [] elemStabMatrixD;
                delete [] elemStabMatrix0;

                for (int i = 0; i < 27; ++i) delete [] jacobianNRMatrix[i];
                delete [] jacobianNRMatrix;

            };//intersect
  
    	};//decomposition

    };//gluezonefine
};


template<>
void Arlequin<2>::setMatVecValuesLagrangeCoarseISO_FEM(){

    double &alpha_f = parametersFine -> getAlphaF();
    double &gamma = parametersFine -> getGamma();
    double& dTime = parametersFine -> getTimeStep();

    double integ = alpha_f * gamma * dTime;

    //Coarse mesh
    for (int l = 0; l < numElemGlueZoneFine; l++){

        int jel =  elementsGlueZoneFine_[l];

        if (domDecompFine.first[jel] == rank) {

            int *connecL = glueZoneFine_[l] -> getNewConnectivity();   

            int numberIntPoints = elementsFine_[jel] -> getNumberOfIntegrationPointsSpecial_ISO();

            std::vector<int> ele, diffElem;
            ele.clear();
            diffElem.clear();

            //Finding coarse elements 
            for (int i=0; i<numberIntPoints; i++){
                int aux = elementsFine_[jel] -> getIntegPointCorrespondenceElement_ISO(i);
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

                //LAGRANGE MULTIPLIERS MATRIXES AND VECTORS
                double **elemMatrixLag0;
                elemMatrixLag0 = new double*[18]();
                for (int i = 0; i < 18; ++i)  elemMatrixLag0[i] = new double[18]();
                double elemVectorLag0_1[18] = {};
                double elemVectorLag0_2[18] = {};

                //tSUPG and tPSPG STABILIZATION
                double **jacobianNRMatrix;
                jacobianNRMatrix = new double*[18]();
                for (int i = 0; i < 18; ++i)  jacobianNRMatrix[i] = new double[18]();
                double rhsVector[18] = {};

                //ARLEQUIN STABILIZATION MATRIX AND VECTOR
                double **elemStabMatrixD;
                elemStabMatrixD = new double*[18]();
                for (int i = 0; i < 18; ++i)  elemStabMatrixD[i] = new double[18]();
                double elemStabVectorD[18] = {};

                double **elemStabMatrix0;
                elemStabMatrix0 = new double*[18]();
                for (int i = 0; i < 18; ++i)  elemStabMatrix0[i] = new double[18]();
                double elemStabVector0[18] = {};

                                     
                int iElemCoarse = diffElem[ielem];

                int *connecC = elementsCoarse_[iElemCoarse] -> getConnectivity();

                
                elementsFine_[jel] -> getLagrangeMultipliersDifferentMesh_ISO_FEM(nodesCoarse_,connecC,iElemCoarse,elemMatrixLag0,
                                                                                 elemVectorLag0_1,elemVectorLag0_2);
                elementsFine_[jel] -> getLagrangeMultipliersDifferentMesh_tSUPG_tPSPG_ISO_FEM(nodesCoarse_,connecC,iElemCoarse, 
                                                                                              jacobianNRMatrix,rhsVector);
                elementsFine_[jel] -> getLagrangeMultipliersDifferentMeshArlqStab_ISO_FEM(nodesCoarse_,connecC,iElemCoarse,elemStabMatrixD, 
                                                                                          elemStabVectorD,elemStabMatrix0,elemStabVector0);
      
                for (int i = 0; i < 9; i++){
                    for (int j = 0; j < 6; j++){
                    
                    //Lagrange multipliers matrixes
                    int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    int dof_j = 2*connecC[j];
                    double value = integ * elemMatrixLag0[2*i][2*j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &value,ADD_VALUES);
                    ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag0[2*i][2*j],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i]+1;
                    dof_j = 2*connecC[j];
                    value = integ * elemMatrixLag0[2*i+1][2*j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &value,ADD_VALUES);
                    ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag0[2*i+1][2*j],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    dof_j = 2*connecC[j] +1 ;
                    value = integ * elemMatrixLag0[2*i][2*j+1];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &value,ADD_VALUES);
                    ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag0[2*i][2*j+1],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1 ;
                    dof_j = 2*connecC[j] +1 ;
                    value = integ * elemMatrixLag0[2*i+1][2*j+1];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &value,ADD_VALUES);
                    ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag0[2*i+1][2*j+1],
                                        ADD_VALUES);

                    //Arlequin Stabilization terms coarse mesh
                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    dof_j = 2*connecC[j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix0[2*i][2*j],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
                    dof_j = 2*connecC[j] + 1;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix0[2*i+1][2*j+1],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    dof_j = 2*NCNumberNodesC + connecC[j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix0[2*i][12+j],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
                    dof_j = 2*NCNumberNodesC + connecC[j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix0[2*i+1][12+j],
                                        ADD_VALUES);

                    };//j

                    //Lagrange multipliers vector
                    int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_2[2*i  ],ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
                    ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_2[2*i+1],ADD_VALUES);

                    //Arlequin Stabilization terms coarse mesh
                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    ierr = VecSetValues(b, 1, &dof_i, &elemStabVector0[2*i  ],ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
                    ierr = VecSetValues(b, 1, &dof_i, &elemStabVector0[2*i+1 ],ADD_VALUES);  



                    
                };//j


                for (int i = 0; i < 6; i++){
                    for (int j = 0; j < 9; j++){
                    
                    //tSUPG stabilization
                    int dof_i = 2*connecC[i];
                    int dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &jacobianNRMatrix[2*i][2*j],
                                        ADD_VALUES);

                    dof_i = 2*connecC[i] + 1 ;
                    dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1 ;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &jacobianNRMatrix[2*i+1][2*j+1],
                                        ADD_VALUES);    

                    //tPSPG stabilization
                    dof_i = 2*NCNumberNodesC + connecC[i];
                    dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &jacobianNRMatrix[12+i][2*j],
                                        ADD_VALUES);

                    dof_i = 2*NCNumberNodesC + connecC[i];
                    dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &jacobianNRMatrix[12+i][2*j+1],
                                        ADD_VALUES);


                    };//j

                    //Lagrange multipliers vector
                    int dof_i = 2*connecC[i];
                    ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_1[2*i  ],ADD_VALUES);
                    dof_i = 2*connecC[i] + 1 ;
                    ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_1[2*i+1],ADD_VALUES);

                    //tSUPG stabilization
                    dof_i = 2*connecC[i];
                    ierr = VecSetValues(b, 1, &dof_i, &rhsVector[2*i  ],ADD_VALUES);
                    dof_i = 2*connecC[i]+ 1 ;
                    ierr = VecSetValues(b, 1, &dof_i, &rhsVector[2*i+1],ADD_VALUES);

                    //tPSPG stabilization
                    dof_i = 2*NCNumberNodesC + connecC[i];
                    ierr = VecSetValues(b, 1, &dof_i, &rhsVector[12+i],ADD_VALUES);
                
                };//j

    
                //Lagrange stabilization terms (Lagrange multipliers defined in the fine mesh)
                for (int i = 0; i < 9; i++){
                    for (int j = 0; j < 9; j++){

                    int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    int dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrixD[2*i][2*j],
                                        ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
                    dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1;
                    ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrixD[2*i+1][2*j+1],
                                        ADD_VALUES);

                    };

                    int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
                    ierr = VecSetValues(b, 1, &dof_i, &elemStabVectorD[2*i  ],ADD_VALUES);

                    dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
                    ierr = VecSetValues(b, 1, &dof_i, &elemStabVectorD[2*i+1 ],ADD_VALUES);  

                };                  
                
                for (int i = 0; i < 18; ++i) {
                delete [] elemMatrixLag0[i];
                delete [] elemStabMatrixD[i];
                delete [] elemStabMatrix0[i];
                delete [] jacobianNRMatrix[i];
                 }
                delete [] elemMatrixLag0; 
                delete [] elemStabMatrixD;
                delete [] elemStabMatrix0;
                delete [] jacobianNRMatrix;
                

            };//intersect
  
        };//decomposition

    };//gluezonefine
};

template<>
void Arlequin<2>::setMatVecValuesLagrangeCoarseISO_ISO(){

	double &alpha_f = parametersFine -> getAlphaF();
	double &gamma = parametersFine -> getGamma();
    double& dTime = parametersFine -> getTimeStep();

    double integ = alpha_f * gamma * dTime;

     //numElemGlueZoneFine
    for (int l = 0; l < numElemGlueZoneFine; l++){

    	int jel =  elementsGlueZoneFine_[l];

    	if (domDecompFine.first[jel] == rank) {

        	int *connecL = glueZoneFine_[l] -> getConnectivity();	


    		int numberIntPoints = elementsFine_[jel] -> getNumberOfIntegrationPointsSpecial_ISO();

    		std::vector<int> ele, diffElem;
        	ele.clear();
        	diffElem.clear();

        	//Finding coarse elements 
        	for (int i=0; i<numberIntPoints; i++){
            	int aux = elementsFine_[jel] -> getIntegPointCorrespondenceElement_ISO(i);
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

 
            	//LAGRANGE MULTIPLIERS MATRIXES AND VECTORS
    			double **elemMatrixLag0;
    			elemMatrixLag0 = new double*[18]();
    			for (int i = 0; i < 18; ++i)  elemMatrixLag0[i] = new double[27]();
    			double elemVectorLag0_1[18] = {};
    			double elemVectorLag0_2[27] = {};

                double **elemStabMatrix;
            	elemStabMatrix = new double*[18]();
           		for (int i = 0; i < 18; ++i)  elemStabMatrix[i] = new double[18]();
            	double elemStabVector[18] = {};

            
            	int iElemCoarse = diffElem[ielem];

            	int *connecC = elementsCoarse_[iElemCoarse] -> getConnectivity();
            	int patch = elementsCoarse_[iElemCoarse] -> getPatch();
            	
            	elementsFine_[jel] -> getLagrangeMultipliersDifferentMesh_ISO_ISO(patch,nodesCoarse_,connecC,IsoParCoarse,iElemCoarse, 
											 								     elemMatrixLag0,elemVectorLag0_1,elemVectorLag0_2,
                                                                                 elemStabMatrix,elemStabVector);

      
            	for (int i = 0; i < 9; i++){
    				for (int j = 0; j < 9; j++){

    				int newconj = nodesCoarse_[connecC[j]] -> getnewcon();	
    				
    				int dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
    				int dof_j = 2*newconj;
    				double value = integ * elemMatrixLag0[2*i][2*j];
    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &value,ADD_VALUES);
    				ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag0[2*i][2*j],
                                        ADD_VALUES);


    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i]+1;
    				dof_j = 2*newconj;
    				value = integ * elemMatrixLag0[2*i+1][2*j];
    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &value,ADD_VALUES);
    				ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag0[2*i+1][2*j],
                                        ADD_VALUES);

    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
    				dof_j = 2*newconj +1 ;
    				value = integ * elemMatrixLag0[2*i][2*j+1];
    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &value,ADD_VALUES);
    				ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag0[2*i][2*j+1],
                                        ADD_VALUES);

    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1 ;
    				dof_j = 2*newconj +1 ;
    				value = integ * elemMatrixLag0[2*i+1][2*j+1];
    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &value,ADD_VALUES);
    				ierr = MatSetValues(A, 1, &dof_j, 1, &dof_i,
                                        &elemMatrixLag0[2*i+1][2*j+1],
                                        ADD_VALUES);

    				// //Arlequin Stabilization Terms
    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
    				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix[2*i][2*j],
                                        ADD_VALUES);

    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
    				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j];
    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix[2*i+1][2*j],
                                        ADD_VALUES);

    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
    				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1;
    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix[2*i][2*j+1],
                                        ADD_VALUES);

    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
    				dof_j = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[j] + 1;
    				ierr = MatSetValues(A, 1, &dof_i, 1, &dof_j,
                                        &elemStabMatrix[2*i+1][2*j+1],
                                        ADD_VALUES);

    				};//i

    				int newconi = nodesCoarse_[connecC[i]] -> getnewcon();

    				int dof_i = 2*newconi;
    				ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_1[2*i  ],ADD_VALUES);
    				dof_i = 2*newconi + 1 ;
    				ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_1[2*i+1],ADD_VALUES);

    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
    				ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_2[2*i  ],ADD_VALUES);

    				dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
    				ierr = VecSetValues(b, 1, &dof_i, &elemVectorLag0_2[2*i+1],ADD_VALUES);

    				//Stabilization Arlequin Term
        			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i];
        			ierr = VecSetValues(b, 1, &dof_i, &elemStabVector[2*i  ],ADD_VALUES);
        			dof_i = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*connecL[i] + 1;
        			ierr = VecSetValues(b, 1, &dof_i, &elemStabVector[2*i+1 ],ADD_VALUES);
    				
    			};//j
    			
    			for (int i = 0; i < 18; ++i) {
            	delete [] elemMatrixLag0[i];
                delete [] elemStabMatrix[i];
           		 }
    			delete [] elemMatrixLag0; 
                delete [] elemStabMatrix;
            };//intersect
  
    	};//decomposition

    };//gluezonefine
}

template<>
int Arlequin<2>::solveArlequinProblem(int iterNumber, double tolerance) {

    std::ofstream dragLift;
    dragLift.open("dragLift.dat", std::ofstream::out | std::ofstream::app);
    if (rank == 0) {
        dragLift << "Time   Pressure Drag   Pressure Lift " 
                 << "Friction Drag  Friction Lift Drag    Lift " 
                 << std::endl;
    };


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


    int sysSize = 3*NCNumberNodesC + 3*NCNumberNodesF + 2*NCnumNodesGlueZoneFine;

    double sumtime = 0;

    int numTimeSteps = fineModel.numTimeSteps;
    double& dTime = parametersFine -> getTimeStep();
    int iTimeStep;

    
 
    for (iTimeStep = 0; iTimeStep < numTimeSteps; iTimeStep++){
 
        
        if (rank == 0) {std::cout << "------------------------- TIME STEP = "
                                  << iTimeStep << " -------------------------"
                                  << std::endl;}


        if (iTimeStep == 10){
        	parametersCoarse -> setSpectralRadius(integScheme);
        	parametersFine -> setSpectralRadius(integScheme);
        }; 



        double &gamma = parametersCoarse ->getGamma();
      
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
            
            
            std::clock_t t1 = std::clock();
            
            ierr = MatCreateAIJ(PETSC_COMM_WORLD, PETSC_DECIDE, PETSC_DECIDE,
                                sysSize,sysSize,10000,NULL,10000,NULL,&A);CHKERRQ(ierr);

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

            //Matrix and vectors - COARSE MESH
            if (elemTypeCoarse == 0) { //FEM mesh
            	setMatVecValuesCoarseFEM();
            } else { //IGA mesh
            	setMatVecValuesCoarseISO();
            }
          	
          	// //Matrix and vectors -FINE MESH  
            if (elemTypeFine == 0) { //FEM mesh
            	setMatVecValuesFineFEM();
            } else { //IGA mesh
            	setMatVecValuesFineISO();
            }
            
            //Matrix and vectors - Lagrange multiplieres - FINE MESH
            if (elemTypeFine == 0){ //FEM fine mesh
            	setMatVecValuesLagrangeFineFEM(iTimeStep);
            } else { //IGA fine mesh
            	setMatVecValuesLagrangeFineISO();
            }
            
            //Matrix and vectors - Lagrange multiplieres - COARSE MESH
            if (elemTypeFine == 0){ //FEM fine mesh
            	if (elemTypeCoarse == 0){ //FEM coarse mesh
            		setMatVecValuesLagrangeCoarseFEM_FEM();
            	} else { //IGA coarse mesh
            		setMatVecValuesLagrangeCoarseFEM_ISO(iTimeStep);
            	};
            } else { //IGA fine mesh
            	if (elemTypeCoarse == 0){ //IGA coarse mesh
            		setMatVecValuesLagrangeCoarseISO_FEM();
            	} else { //IGA coarse mesh
            		setMatVecValuesLagrangeCoarseISO_ISO();
            	};
            };          
	        

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
        
            	
                if (elemTypeCoarse == 0) { //FEM mesh

                    int newconi = i;

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

                } else { //IGA mesh

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
 
            };	
             


            for (int i = 0; i < numNodesFine; ++i){


            	if (elemTypeFine == 0){ //FEM mesh

	                Ii = 2*i + 3*NCNumberNodesC;
	                ierr = VecGetValues(All, Ione, &Ii, &val);
	                u_ = val;
	                nodesFine_[i] -> incrementAcceleration(0,u_);
	                nodesFine_[i] -> incrementVelocity(0,u_*gamma*dTime);
	                normU += val*val;
	            
	                Ii = 2*i + 1 + 3*NCNumberNodesC;
	                ierr = VecGetValues(All, Ione, &Ii, &val);
	                u_ = val;
	                nodesFine_[i] -> incrementAcceleration(1,u_);
	                nodesFine_[i] -> incrementVelocity(1,u_*gamma*dTime);
	                normU += val*val;

	                Ii = 2*NCNumberNodesF + i + 3*NCNumberNodesC;
	                ierr = VecGetValues(All,Ione,&Ii,&val);
	                p_ = val;
	                nodesFine_[i] -> incrementPressure(p_);
	                normP += val*val;

            	} else { //IGA mesh

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
	                
            };

            

            for (int i = 0; i < NCnumNodesGlueZoneFine; ++i){

                if (elemTypeFine == 0){ //fem fine mesh

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

                } else { //IGA fine mesh

                    
                    for (int j = 0; j < numNodesGlueZoneFine; j++){

                        if (nodesFine_[nodesGlueZoneFine_[j]] -> getnewcon() == NCnodesGlueZoneFine_[i]){
                            Ii = 3*NCNumberNodesF + 3*NCNumberNodesC + 2*i;
                            ierr = VecGetValues(All, Ione, &Ii, &val);
                            lag_ = val;
                            normL += val * val;
                            nodesFine_[nodesGlueZoneFine_[j]] -> incrementLagrangeMultiplier(0,lag_);

                            Ii = 3*NCNumberNodesF + 3*NCNumberNodesC + 2*i + 1;
                            ierr = VecGetValues(All, Ione, &Ii, &val);
                            lag_ = val;
                            normL += val * val;
                            nodesFine_[nodesGlueZoneFine_[j]] -> incrementLagrangeMultiplier(1,lag_);
                        };
                    };

                };

            };
			           
   
            std::clock_t t2 = std::clock();


            sumtime += 1000.*(t2-t1)/CLOCKS_PER_SEC/1000.;

         
            if(rank == 0){
                                
                std::cout << "Iteration = " << inewton 
                          << " (" << iterations << ")"  
                          << "   Du Norm = " << std::scientific << sqrt(normU) 
                          << " " << sqrt(normP)  << " " << sqrt(normL)
                          << "  Time (s) = " << std::fixed
                          << 1000.*(t2-t1)/CLOCKS_PER_SEC/1000. << std::endl;
            };
                      
            ierr = KSPDestroy(&ksp); CHKERRQ(ierr);
            ierr = VecDestroy(&b); CHKERRQ(ierr);
            ierr = VecDestroy(&u); CHKERRQ(ierr);
            ierr = VecDestroy(&All); CHKERRQ(ierr);
            ierr = MatDestroy(&A); CHKERRQ(ierr);
            if(iTimeStep<5){normU=0.;}
            if(iTimeStep<100){if(inewton>1)normU=0.;}
            if (sqrt(normU) <= tolerance) {
                break;
            };

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

      
        
        for (int i = 0; i<numNodesGlueZoneFine; i++){

            int elCoarse = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalElemCorrespondence();
            double* xsiC = nodesFine_[nodesGlueZoneFine_[i]] -> getNodalXsiCorrespondence();


            int *connecCoarse = elementsCoarse_[elCoarse] -> getConnectivity();


            if (elemTypeCoarse == 0){

            	double u_coarse[6], v_coarse[6], p_coarse[6];
	            for (int j=0; j<6; j++){
	                u_coarse[j] = nodesCoarse_[connecCoarse[j]] -> getVelocity(0);
	                v_coarse[j] = nodesCoarse_[connecCoarse[j]] -> getVelocity(1);
	                p_coarse[j] = nodesCoarse_[connecCoarse[j]] -> getPressure();
	            };

	            double phiC_[6];
	            shapeQuad.evaluateFem(xsiC,phiC_);
	            
	            double u = 0.;
	            double v = 0.;
	            double p = 0.;
	            for (int j=0; j<6; j++){
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

            } else {

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

            for (int i=0; i<numNodesCoarse; i++){
        		for (int k = 0; k < dim; k++) nodesCoarse_[i] -> setVelocityArlequin(k,nodesCoarse_[i] -> getVelocity(k));
            	nodesCoarse_[i] -> setPressureArlequin(nodesCoarse_[i] -> getPressure());
        	};
        };
       
        
        //Compute and print drag and lift coefficients
        if (fineModel.computeDragAndLift){
            if (elemTypeFine == 0){
                dragAndLiftCoefficientsFEM(dragLift,iTimeStep);
            } else {
                dragAndLiftCoefficientsISO(dragLift,iTimeStep);
            };
        };

        // Printing results
        printResults(iTimeStep);
        
    };

    PetscFree(dof);
    return 0;
    
};
