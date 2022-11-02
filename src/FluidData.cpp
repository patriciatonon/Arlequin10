#include "FluidData.h"

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//---------------------SUBDIVIDES THE FINITE ELEMENT DOMAIN---------------------
//------------------------------------------------------------------------------
template<>
void FluidData<2>::domainDecompositionMETIS_FEM(std::vector<Elements *> &elem_, int &numFemElem, int &numFemNodes) {
    
    std::string mirror2;
    mirror2 = "domain_decomposition.txt";
    std::ofstream mirrorData(mirror2.c_str());
    
    int size;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    idx_t objval;
    idx_t numEl = numFemElem;
    idx_t numNd = numFemNodes;
    idx_t dd = 2;
    idx_t ssize = size;
    idx_t three = 3;
    idx_t one = 1;
    idx_t elem_start[numEl+1], elem_connec[6*numEl];
   

    MPI_Bcast(&numEl,1,MPI_INT,0,PETSC_COMM_WORLD);
    MPI_Bcast(&numNd,1,MPI_INT,0,PETSC_COMM_WORLD);

    part_elem = new idx_t[numEl];
    part_nodes = new idx_t[numNd];

    if (rank == 0){

        for (idx_t i = 0; i < numEl+1; i++){
            elem_start[i]=6*i;
        };
        for (idx_t jel = 0; jel < numEl; jel++){
            int *connec;
            connec = elem_[jel] -> getConnectivity();      
            for (idx_t i=0; i<6; i++){
                elem_connec[6*jel+i] = connec[i];
            };
        };
    

        //Performs the domain decomposition
        if (size == 1){

            for(int i = 0; i < numFemElem; i++){
                part_elem[i] = 0;
            };

            for(int i = 0; i < numFemNodes; i++){
                part_nodes[i] = 0;
            };
        
        } else {
            
            METIS_PartMeshDual(&numEl, &numNd, elem_start, elem_connec,
                               NULL, NULL, &one, &ssize, NULL, NULL,
                               &objval, part_elem, part_nodes);
        };
        

        mirrorData << std::endl 
                   << "FLUID MESH DOMAIN DECOMPOSITION - ELEMENTS" << std::endl;
        for(int i = 0; i < numFemElem; i++){
            mirrorData << "process = " << part_elem[i]
                       << ", element = " << i << std::endl;
        };

        mirrorData << std::endl 
                   << "FLUID MESH DOMAIN DECOMPOSITION - NODES" << std::endl;
        for(int i = 0; i < numFemNodes; i++){
            mirrorData << "process = " << part_nodes[i]
                       << ", node = " << i << std::endl;
        };
        
    }

    MPI_Bcast(part_elem,numEl,MPI_INT,0,PETSC_COMM_WORLD);
    MPI_Bcast(part_nodes,numNd,MPI_INT,0,PETSC_COMM_WORLD);

    return;

};

template<>
void FluidData<2>::domainDecompositionMETIS_ISO(std::vector<Elements *> &elem_, int &numIsoElem, int &numCP) {
    
    std::string mirror2;
    mirror2 = "domain_decomposition.txt";
    std::ofstream mirrorData(mirror2.c_str());
    
    int size;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    idx_t objval;
    idx_t numEl = numIsoElem;
    idx_t numNd = numCP;
    idx_t dd = 2;
    idx_t ssize = size;
    idx_t three = 3;
    idx_t one = 1;
    idx_t elem_start[numEl+1], elem_connec[9*numEl];
    

    MPI_Bcast(&numEl,1,MPI_INT,0,PETSC_COMM_WORLD);
    MPI_Bcast(&numNd,1,MPI_INT,0,PETSC_COMM_WORLD);

    part_elem = new idx_t[numEl];
    part_nodes = new idx_t[numNd];

    if (rank == 0){

        for (idx_t i = 0; i < numEl+1; i++){
            elem_start[i]=9*i;
        };
        for (idx_t jel = 0; jel < numEl; jel++){
            int *connec;
            connec = elem_[jel] -> getConnectivity();      
            for (idx_t i=0; i<9; i++){
                elem_connec[9*jel+i] = connec[i];
            };
        };
    

        //Performs the domain decomposition

        if (size == 1){

            for(int i = 0; i < numIsoElem; i++){
                 part_elem[i] = 0;
            };

            for(int i = 0; i < numCP; i++){
                part_nodes[i] = 0;
            };
            
        } else {
            
            METIS_PartMeshDual(&numEl, &numNd, elem_start, elem_connec,
                               NULL, NULL, &one, &ssize, NULL, NULL,
                               &objval, part_elem, part_nodes);
        };
        

        mirrorData << std::endl 
                   << "FLUID MESH DOMAIN DECOMPOSITION - ELEMENTS" << std::endl;
        for(int i = 0; i < numIsoElem; i++){
            mirrorData << "process = " << part_elem[i]
                       << ", element = " << i << std::endl;
        };

        mirrorData << std::endl 
                   << "FLUID MESH DOMAIN DECOMPOSITION - NODES" << std::endl;
        for(int i = 0; i < numCP; i++){
            mirrorData << "process = " << part_nodes[i]
                       << ", node = " << i << std::endl;
        };
        
    }

    MPI_Bcast(part_elem,numEl,MPI_INT,0,PETSC_COMM_WORLD);
    MPI_Bcast(part_nodes,numNd,MPI_INT,0,PETSC_COMM_WORLD);

    return;

};
//------------------------------------------------------------------------------
//----------------------------BEZIER CONNECTIVITY-------------------------------
//------------------------------------------------------------------------------
template<>
void FluidData<2>::BezierConnectivity(int &numPatches) {

    int numNodes = 0;
    int k = 0;
    for (int ipatch = 0; ipatch < numPatches; ipatch++){
        int numElemx = IsoPar_[ipatch] -> getNumElemPatchx();
        int numElemy = IsoPar_[ipatch] -> getNumElemPatchy();
        int  Beconnec[numElemx*numElemy][9];
        for (int j = 0; j< 2*numElemy; j++){    
            for (int i = 0; i< numElemx; i++){
                int rest = j%2;            
                if (j == 0) {
                    if (i == 0) Beconnec[j*numElemx+i][0] = numNodes++;  
                    Beconnec[j*numElemx+i][1] = numNodes++;
                    Beconnec[j*numElemx+i][2] = numNodes;
                    if (i != numElemx-1) Beconnec[j*numElemx+i+1][0] = numNodes;
                    numNodes ++;
                }
                if (j != 0) {
                    if (rest == 0) {
                        int kp = (j/2);
                        if (i == 0) {
                            Beconnec[kp*numElemx+i][0] = numNodes; 
                            Beconnec[(kp-1)*numElemx+i][6] = numNodes++;  
                        }
                        Beconnec[kp*numElemx+i][1] = numNodes;
                        Beconnec[(kp-1)*numElemx+i][7] = numNodes++;
                        Beconnec[kp*numElemx+i][2] = numNodes;
                        Beconnec[(kp-1)*numElemx+i][8] = numNodes;  
                        if (i != numElemx-1) {
                            Beconnec[kp*numElemx+i+1][0] = numNodes; 
                            Beconnec[(kp-1)*numElemx+i+1][6] = numNodes;       
                        }       
                        numNodes ++;  
                    } else {
                        int kp = ((j-1)/2);
                        if (i == 0) Beconnec[kp*numElemx+i][3] = numNodes++;
                        Beconnec[kp*numElemx+i][4] = numNodes++;
                        Beconnec[kp*numElemx+i][5] = numNodes;
                        if (i != numElemx -1) Beconnec[kp*numElemx+i+1][3] = numNodes;
                        numNodes ++;
                    }
                } //j!=0
            } //elemx
        }//elemy
        for (int i = 0; i < numElemx; i++) {
            if (i == 0) Beconnec[(numElemy-1)*numElemx+i][6] = numNodes++;  
            Beconnec[(numElemy-1)*numElemx+i][7] = numNodes++;
            Beconnec[(numElemy-1)*numElemx+i][8] = numNodes;
            if (i != numElemx-1) Beconnec[(numElemy-1)*numElemx+i+1][6] = numNodes;
            numNodes ++;
        }
        int contel = 0;
        for (int aux = 0; aux < ipatch; aux++){
             int numElemxx = IsoPar_[aux] -> getNumElemPatchx();
             int numElemyy = IsoPar_[aux] -> getNumElemPatchy();
             contel +=  numElemxx*numElemyy;
        }   

        for (int iElem = 0; iElem < numElemx*numElemy; iElem++){
            int *Becon_;
            Becon_ = new int[9];
            for (int j=0; j<9; j++){
                Becon_[j] = Beconnec[iElem][j];
            }
            elements_[iElem+contel] -> setBezierConnectivity(Becon_);
        }
    }//patch
    NumBezierNodes = numNodes;
};

template<>
void FluidData<3>::BezierConnectivity(int &numPatches) {

  //   int numNodes = 0;
  //   int kne = 0;
  //   for (int ipatch = 0; ipatch < numPatches; ipatch++){
  //       int numElemx = IsoPar_[ipatch] -> getNumElemPatchx();
  //       int numElemy = IsoPar_[ipatch] -> getNumElemPatchy();
  //       int numElemz = IsoPar_[ipatch] -> getNumElemPatchz();
  //       int  Beconnec[numElemx*numElemy*numElemz][27];
        
  //       for (int k = 0; k <= 2*numElemz; k++){
        	
  //       	if (k == 0){
		//         for (int j = 0; j< 2*numElemy; j++){    
		//             for (int i = 0; i< numElemx; i++){
		//                 int rest = j%2;                
		//                 if (j == 0) {
		//                     if (i == 0) Beconnec[j*numElemx+i][0] = numNodes++;  
		//                     Beconnec[j*numElemx+i][1] = numNodes++;
		//                     Beconnec[j*numElemx+i][2] = numNodes;
		//                     if (i != numElemx-1) Beconnec[j*numElemx+i+1][0] = numNodes;
		//                     numNodes ++;
		//                 }
		//                 if (j != 0) {
		//                     if (rest == 0) {
		//                         int kp = (j/2);
		//                         if (i == 0) {
		//                             Beconnec[kp*numElemx+i][0] = numNodes; 
		//                             Beconnec[(kp-1)*numElemx+i][6] = numNodes++;  
		//                         }
		//                         Beconnec[kp*numElemx+i][1] = numNodes;
		//                         Beconnec[(kp-1)*numElemx+i][7] = numNodes++;
		//                         Beconnec[kp*numElemx+i][2] = numNodes;
		//                         Beconnec[(kp-1)*numElemx+i][8] = numNodes;   
		//                         if (i != numElemx-1) {
		//                             Beconnec[kp*numElemx+i+1][0] = numNodes; 
		//                             Beconnec[(kp-1)*numElemx+i+1][6] = numNodes;       
		//                         }       
		//                         numNodes ++;
		//                     } else {
		//                         int kp = ((j-1)/2);
		//                         if (i == 0) Beconnec[kp*numElemx+i][3] = numNodes++;
		//                         Beconnec[kp*numElemx+i][4] = numNodes++;
		//                         Beconnec[kp*numElemx+i][5] = numNodes;
		//                         if (i != numElemx -1) Beconnec[kp*numElemx+i+1][3] = numNodes;
		//                         numNodes ++;
		//                     }
		//                 } //j!=0
		//             } //elemx
		//         }//elemy
		//         for (int i = 0; i < numElemx; i++) {
		//             if (i == 0) Beconnec[(numElemy-1)*numElemx+i][6] = numNodes++;  
		//             Beconnec[(numElemy-1)*numElemx+i][7] = numNodes++;
		//             Beconnec[(numElemy-1)*numElemx+i][8] = numNodes;
		//             if (i != numElemx-1) Beconnec[(numElemy-1)*numElemx+i+1][6] = numNodes;
		//             numNodes ++;
		//         }
  //       	}


  //       	if ((k != 0) && (k != 2*numElemz)){
  //       		int restk = k%2;
        		
  //       		if (restk != 0){
  //       			int ind = (k-1)/2;
		//     		for (int j = 0; j< 2*numElemy; j++){    
		// 	            for (int i = 0; i< numElemx; i++){
		// 	                int rest = j%2;               
		// 	                if (j == 0) {
		// 	                    if (i == 0) Beconnec[ind*numElemx*numElemy + j*numElemx + i][9] = numNodes++;  
		// 	                    Beconnec[ind*numElemx*numElemy + j*numElemx+i][10] = numNodes++;
		// 	                    Beconnec[ind*numElemx*numElemy + j*numElemx+i][11] = numNodes;
		// 	                    if (i != numElemx-1) Beconnec[ind*numElemx*numElemy + j*numElemx+i+1][9] = numNodes;
		// 	                    numNodes ++;
		// 	                }
		// 	                if (j != 0) {
		// 	                    if (rest == 0) {
		// 	                        int kp = (j/2);
		// 	                        if (i == 0) {
		// 	                            Beconnec[ind*numElemx*numElemy + kp*numElemx+i][9] = numNodes; 
		// 	                            Beconnec[ind*numElemx*numElemy + (kp-1)*numElemx+i][15] = numNodes++;  
		// 	                        }
		// 	                        Beconnec[ind*numElemx*numElemy + kp*numElemx+i][10] = numNodes;
		// 	                        Beconnec[ind*numElemx*numElemy + (kp-1)*numElemx+i][16] = numNodes++;
		// 	                        Beconnec[ind*numElemx*numElemy + kp*numElemx+i][11] = numNodes;
		// 	                        Beconnec[ind*numElemx*numElemy + (kp-1)*numElemx+i][17] = numNodes;			                            
		// 	                        if (i != numElemx-1) {
		// 	                            Beconnec[ind*numElemx*numElemy + kp*numElemx+i+1][9] = numNodes; 
		// 	                            Beconnec[ind*numElemx*numElemy + (kp-1)*numElemx+i+1][15] = numNodes;       
		// 	                        }       
		// 	                        numNodes ++;			                        
		// 	                    } else {
		// 	                        int kp = ((j-1)/2);
		// 	                        if (i == 0) Beconnec[ind*numElemx*numElemy + kp*numElemx+i][12] = numNodes++;
		// 	                        Beconnec[ind*numElemx*numElemy + kp*numElemx+i][13] = numNodes++;
		// 	                        Beconnec[ind*numElemx*numElemy + kp*numElemx+i][14] = numNodes;
		// 	                        if (i != numElemx -1) Beconnec[ind*numElemx*numElemy + kp*numElemx+i+1][12] = numNodes;
		// 	                        numNodes ++;
		// 	                    }			                
		// 	                } //j!=0
		// 	            } //elemx
		// 	        }//elemy			      
		// 	        for (int i = 0; i < numElemx; i++) {
		// 	            if (i == 0) Beconnec[ind*numElemx*numElemy + (numElemy-1)*numElemx+i][15] = numNodes++;  
		// 	            Beconnec[ind*numElemx*numElemy + (numElemy-1)*numElemx+i][16] = numNodes++;
		// 	            Beconnec[ind*numElemx*numElemy + (numElemy-1)*numElemx+i][17] = numNodes;
		// 	            if (i != numElemx-1) Beconnec[ind*numElemx*numElemy + (numElemy-1)*numElemx+i+1][15] = numNodes;
		// 	            numNodes ++;
		// 	        }
  //       		} // if restz !=0
        		
  //       		if (restk == 0){
  //       			int ind = (k/2);
		//     		for (int j = 0; j< 2*numElemy; j++){    
		// 	            for (int i = 0; i< numElemx; i++){
		// 	                int rest = j%2;			                                
		// 	                if (j == 0) {
		// 	                    if (i == 0) {
		// 	                    	Beconnec[ind*numElemx*numElemy + j*numElemx+i][0] = numNodes;  
		// 	                    	Beconnec[(ind-1)*numElemx*numElemy + j*numElemx+i][18] = numNodes++; 
		// 	                    }
		// 	                    Beconnec[ind*numElemx*numElemy + j*numElemx+i][1] = numNodes;
		// 	                    Beconnec[(ind-1)*numElemx*numElemy + j*numElemx+i][19] = numNodes++;
		// 	                    Beconnec[ind*numElemx*numElemy + j*numElemx+i][2] = numNodes;
		// 	                    Beconnec[(ind-1)*numElemx*numElemy + j*numElemx+i][20] = numNodes;
		// 	                    if (i != numElemx-1) {
		// 	                    	Beconnec[ind*numElemx*numElemy + j*numElemx+i+1][0] = numNodes;
		// 	                    	Beconnec[(ind-1)*numElemx*numElemy + j*numElemx+i+1][18] = numNodes;
		// 	                    }
		// 	                    numNodes ++;
		// 	                }
		// 	                if (j != 0){
		// 	                    if (rest == 0) {
		// 	                        int kp = (j/2);
		// 	                        if (i == 0) {
		// 	                            Beconnec[ind*numElemx*numElemy + kp*numElemx+i][0] = numNodes; 
		// 	                            Beconnec[(ind-1)*numElemx*numElemy + kp*numElemx+i][18] = numNodes;
		// 	                            Beconnec[ind*numElemx*numElemy +(kp-1)*numElemx+i][6] = numNodes; 
		// 	                            Beconnec[(ind-1)*numElemx*numElemy + (kp-1)*numElemx+i][24] = numNodes++;  
		// 	                        }
		// 	                        Beconnec[ind*numElemx*numElemy + kp*numElemx+i][1] = numNodes;
		// 	                        Beconnec[(ind-1)*numElemx*numElemy + kp*numElemx+i][19] = numNodes;
		// 	                        Beconnec[ind*numElemx*numElemy +(kp-1)*numElemx+i][7] = numNodes;
		// 	                        Beconnec[(ind-1)*numElemx*numElemy +(kp-1)*numElemx+i][25] = numNodes++;
		// 	                        Beconnec[ind*numElemx*numElemy +kp*numElemx+i][2] = numNodes;
		// 	                        Beconnec[(ind-1)*numElemx*numElemy + kp*numElemx+i][20] = numNodes;
		// 	                        Beconnec[ind*numElemx*numElemy +(kp-1)*numElemx+i][8] = numNodes;
		// 	                        Beconnec[(ind-1)*numElemx*numElemy +(kp-1)*numElemx+i][26] = numNodes;			                            
		// 	                        if (i != numElemx-1) {
		// 	                            Beconnec[ind*numElemx*numElemy +kp*numElemx+i+1][0] = numNodes; 
		// 	                            Beconnec[(ind-1)*numElemx*numElemy + kp*numElemx+i+1][18] = numNodes; 
		// 	                            Beconnec[ind*numElemx*numElemy +(kp-1)*numElemx+i+1][6] = numNodes; 
		// 	                            Beconnec[(ind-1)*numElemx*numElemy + (kp-1)*numElemx+i+1][24] = numNodes;       
		// 	                        }       
		// 	                        numNodes ++;			                        
		// 	                    } else {
		// 	                        int kp = ((j-1)/2);
		// 	                        if (i == 0) {
		// 								Beconnec[ind*numElemx*numElemy +kp*numElemx+i][3] = numNodes;
		// 								Beconnec[(ind-1)*numElemx*numElemy +kp*numElemx+i][21] = numNodes++;
		// 	                        } 
		// 	                        Beconnec[ind*numElemx*numElemy + kp*numElemx+i][4] = numNodes;
		// 	                        Beconnec[(ind-1)*numElemx*numElemy + kp*numElemx+i][22] = numNodes++;
		// 	                        Beconnec[ind*numElemx*numElemy + kp*numElemx+i][5] = numNodes;
		// 	                        Beconnec[(ind-1)*numElemx*numElemy + kp*numElemx+i][23]= numNodes;
		// 	                        if (i != numElemx -1) {
		// 	                        	Beconnec[ind*numElemx*numElemy + kp*numElemx+i+1][3] = numNodes;
		// 	                        	Beconnec[(ind-1)*numElemx*numElemy + kp*numElemx+i+1][21] = numNodes;
		// 	                        }
		// 	                        numNodes ++;
		// 	                    } //else
		// 	                } //j!=0
		// 	            } //elemx
		// 	        }//elemy
		// 	        for (int i = 0; i < numElemx; i++) {
		// 	            if (i == 0) {
		// 	            	Beconnec[ind*numElemx*numElemy + (numElemy-1)*numElemx+i][6] = numNodes;  
		// 	            	Beconnec[(ind-1)*numElemx*numElemy +(numElemy-1)*numElemx+i][24] = numNodes++;  
		// 	            }
		// 	            Beconnec[ind*numElemx*numElemy + (numElemy-1)*numElemx+i][7] = numNodes;
		// 	            Beconnec[(ind-1)*numElemx*numElemy + (numElemy-1)*numElemx+i][25] = numNodes++;
		// 	            Beconnec[ind*numElemx*numElemy + (numElemy-1)*numElemx+i][8] = numNodes;
		// 	            Beconnec[(ind-1)*numElemx*numElemy + (numElemy-1)*numElemx+i][26] = numNodes;
		// 	            if (i != numElemx-1) {
		// 	            	Beconnec[ind*numElemx*numElemy + (numElemy-1)*numElemx+i+1][6] = numNodes;
		// 	            	Beconnec[(ind-1)*numElemx*numElemy + (numElemy-1)*numElemx+i+1][24] = numNodes;
		// 	            }
		// 	            numNodes ++;
		// 	        }    
		// 	    } // if restz == 0	

	 //       	} // k != 0 or k != 2*numElemz
        	


  //       	if (k == 2*numElemz){
  //       		int ind = numElemz - 1;
  //       		for (int j = 0; j< 2*numElemy; j++){    
		//             for (int i = 0; i< numElemx; i++){
		//                 int rest = j%2;               
		//                 if (j == 0) {
		//                     if (i == 0) Beconnec[ind*numElemx*numElemy +j*numElemx+i][18] = numNodes++;  
		//                     Beconnec[ind*numElemx*numElemy +j*numElemx+i][19] = numNodes++;
		//                     Beconnec[ind*numElemx*numElemy +j*numElemx+i][20] = numNodes;
		//                     if (i != numElemx-1) Beconnec[ind*numElemx*numElemy +j*numElemx+i+1][18] = numNodes;
		//                     numNodes ++;
		//                 }
		//                 if (j != 0) {
		//                     if (rest == 0) {
		//                         int kp = (j/2);
		//                         if (i == 0) {
		//                             Beconnec[ind*numElemx*numElemy +kp*numElemx+i][18] = numNodes; 
		//                             Beconnec[ind*numElemx*numElemy +(kp-1)*numElemx+i][24] = numNodes++;  
		//                         }
		//                         Beconnec[ind*numElemx*numElemy +kp*numElemx+i][19] = numNodes;
		//                         Beconnec[ind*numElemx*numElemy +(kp-1)*numElemx+i][25] = numNodes++;
		//                         Beconnec[ind*numElemx*numElemy +kp*numElemx+i][20] = numNodes;
		//                         Beconnec[ind*numElemx*numElemy +(kp-1)*numElemx+i][26] = numNodes;
		//                         if (i != numElemx-1) {
		//                             Beconnec[ind*numElemx*numElemy +kp*numElemx+i+1][18] = numNodes; 
		//                             Beconnec[ind*numElemx*numElemy +(kp-1)*numElemx+i+1][24] = numNodes;       
		//                         }       
		//                         numNodes ++;    
		//                     } else {
		//                         int kp = ((j-1)/2);
		//                         if (i == 0) Beconnec[ind*numElemx*numElemy +kp*numElemx+i][21] = numNodes++;
		//                         Beconnec[ind*numElemx*numElemy +kp*numElemx+i][22] = numNodes++;
		//                         Beconnec[ind*numElemx*numElemy +kp*numElemx+i][23] = numNodes;
		//                         if (i != numElemx -1) Beconnec[ind*numElemx*numElemy +kp*numElemx+i+1][21] = numNodes;
		//                         numNodes ++;
		//                     }
		//                 } //j!=0
		//             } //elemx
		//         }//elemy
		//         for (int i = 0; i < numElemx; i++) {
		//             if (i == 0) Beconnec[ind*numElemx*numElemy +(numElemy-1)*numElemx+i][24] = numNodes++;  
		//             Beconnec[ind*numElemx*numElemy +(numElemy-1)*numElemx+i][25] = numNodes++;
		//             Beconnec[ind*numElemx*numElemy +(numElemy-1)*numElemx+i][26] = numNodes;
		//             if (i != numElemx-1) Beconnec[ind*numElemx*numElemy +(numElemy-1)*numElemx+i+1][24] = numNodes;
		//             numNodes ++;
		//         }
  //       	}//if k = 2*numElemz

  //       } // loop in z direction
		
		// int contel = 0;
  //       for (int aux = 0; aux < ipatch; aux++){
  //            int numElemxx = IsoPar_[ipatch] -> getNumElemPatchx();
  //            int numElemyy = IsoPar_[ipatch] -> getNumElemPatchy();
  //            int numElemzz = IsoPar_[ipatch] -> getNumElemPatchz();
  //            contel +=  numElemxx*numElemyy*numElemzz;
  //       } 
      

  //       for (int iElem = 0; iElem < numElemx*numElemy*numElemz; ++iElem){
  //           int *Becon_;
  //           Becon_ = new int[27];
  //           for (int i=0; i<27; i++) Becon_[i] = 0;
  //           if (part_elem[iElem+contel] == rank){
  //               for (int j=0; j<27; j++){
  //                   Becon_[j] = Beconnec[iElem][j];
  //               }
  //               elements_[kne] -> setBezierConnectivity(Becon_);
  //               kne ++;
  //           }
  //       }
    
  //   }//patch



  //   NumBezierNodes = numNodes;
};


//------------------------------------------------------------------------------
//----------------------------READS FLUID INPUT FILE----------------------------
//------------------------------------------------------------------------------
template<>
void FluidData<2>::dataReading_FEM(const std::string& inputFile,const std::string& inputMeshFem, 
                                   const std::string& mirror,const bool& deleteFiles){

    //variables
    //Fluid data
    double pressInf;       			     //Undisturbed pressure 
    double rhoInf;         			     //Density
    double viscInf;        				 //Viscosity
    double fieldForces[3]; 				 //Field forces (constant) 
    double velocityInf[3];               //Infinity velocity
    int numTimeSteps;                    //Number of Time Steps
    double dTime;                        //Time Step
    int printFreq;         				 //Printing frequence of output files
    //Arlequin Problem
    double arlequinK1;                   //L2
    double arlequinK2;                   //H1
    double glueZoneThickness;			 //Thickness from gluing zone
    double arlequinEpsilon;				 //Constant > 0
    //Time integration parameters
    double integScheme;                  //Time Integration Scheme (0 - max. dissipation; 1 - no dissipation)
    double alpha_f;						 //Alphageneralized parameter - matrix except mass
    double alpha_m;						 //Alphageneralized parameter - mass
    double gamma;						 //Alphageneralized parameter 
    //Data meshes
    int numFemNodes;                     //Number of nodes in FEM mesh
    int numFemElem;                      //Number of FEM elements 
    int numFemBoundaries;                //Number of fluid boundaries (FEM MESH)
    int numFemBoundElem;                 //Number of elements in fluid boundaries (FEM MESH)
                           
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);      
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    
    int dim = 2;
       
    if (rank == 0)std::cout << "Reading fluid data from \"" 
                            << inputFile << "\"" << std::endl;
                    

    //Defines input and output files
    std::ifstream inputData(inputFile.c_str());
    std::ifstream file(inputMeshFem);
    std::ofstream mirrorData(mirror.c_str());
    std::string line;
    std::string remove2 = "rm ";
      
    //Start to read the data problems  
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
    //Read and print out time steps and printing frequence
    inputData >> numTimeSteps >> printFreq;
    
    mirrorData << "Number of Time Steps   = " << numTimeSteps << std::endl;
    mirrorData << "Printing Frequence     = " << printFreq << std::endl;
    
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read and print out data from Arlequin problem
    inputData >> glueZoneThickness;
    mirrorData << "Thickness of glue zone   = " << glueZoneThickness << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    inputData >> arlequinEpsilon;
    mirrorData << "Epsilon from Arlequin formulation   = " << arlequinEpsilon << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    inputData >> arlequinK1 >> arlequinK2;
    mirrorData << "Arlequin k1   = " << arlequinK1 << std::endl;
    mirrorData << "Arlequin k2   = " << arlequinK2 << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
    //Read and print out undisturbed velocity and pressure components
    inputData >> velocityInf[0] >> velocityInf[1] >> velocityInf[2] >> pressInf;
    
    mirrorData << "Undisturbed Velocity x = " << velocityInf[0] << std::endl;
    mirrorData << "Undisturbed Velocity y = " << velocityInf[1] << std::endl;
    mirrorData << "Undisturbed Velocity z = " << velocityInf[2] << std::endl;
    mirrorData << "Undisturbed Pressure   = " << pressInf << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read and print out undisturbed density and viscosity
    inputData >> rhoInf >> viscInf;
    
    mirrorData << "Undisturbed Density    = " << rhoInf << std::endl;
    mirrorData << "Undisturbed Viscosity  = " << viscInf << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read and print out time step and time integration scheme
    inputData >> dTime >> integScheme;
    
    mirrorData << "Time Step              = " << dTime << std::endl;
    mirrorData << "Time Integration Scheme= " << integScheme << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);

    //Read and print out field forces
    inputData >> fieldForces[0] >> fieldForces[1] >> fieldForces[2];
    
    mirrorData << "Field Forces x         = " << fieldForces[0] << std::endl;
    mirrorData << "Field Forces y         = " << fieldForces[1] << std::endl;
    mirrorData << "Field Forces z         = " << fieldForces[2] << std::endl 
               << std::endl;
    
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    fluidParameters.setViscosity(viscInf);
    fluidParameters.setDensity(rhoInf);
    fluidParameters.setTimeStep(dTime);
    fluidParameters.setSpectralRadius(integScheme);
    fluidParameters.setFieldForce(fieldForces);
    fluidParameters.setArlequink1(arlequinK1);
    fluidParameters.setArlequink2(arlequinK2);
    fluidParameters.setVelocityInf(velocityInf);
    fluidParameters.setNumTimeSteps(numTimeSteps);
    fluidParameters.setFreqPrint(printFreq);
    fluidParameters.setArlequinEpsilon(arlequinEpsilon);
    fluidParameters.setGlueZoneThickness(glueZoneThickness);


    //Read and print out Drag and lift coeficients
    inputData >> computeDragAndLift >> numberOfLines; 

    dragAndLiftBoundary.reserve(numberOfLines);
    for (int i = 0; i < numberOfLines; ++i)
    {
        int aux;
        inputData >> aux; 
        dragAndLiftBoundary.push_back(aux);
    }


    mirrorData << "Compute Drag and Lift  = " << computeDragAndLift<< std::endl;
    mirrorData << "Number of Lines  = " << numberOfLines << std::endl;
    for (int i = 0; i < numberOfLines; ++i)
    {
        mirrorData << "Lines  = " << dragAndLiftBoundary[i] << std::endl;
    }

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read and print out Printing results
    inputData >> printVelocity;              getline(inputData,line);
    inputData >> printRealVelocity;          getline(inputData,line);
    inputData >> printPressure;              getline(inputData,line);
    inputData >> printRealPressure;          getline(inputData,line);
    inputData >> printMeshVelocity;          getline(inputData,line);
    inputData >> printMeshDisplacement;      getline(inputData,line);
    inputData >> printDistFunction;          getline(inputData,line);
    inputData >> printGlueZone;              getline(inputData,line);
    inputData >> printEnergyWeightFunction;  getline(inputData,line);
    inputData >> printLagrangeMultipliers;   getline(inputData,line);
    inputData >> printNodalCorrespondence;   getline(inputData,line);
    inputData >> printProcess;          

    mirrorData << "PrintVelocity              = " << printVelocity << std::endl;
    mirrorData << "PrintRealVelocity          = " << printRealVelocity << std::endl;
    mirrorData << "PrintPressure              = " << printPressure << std::endl;
    mirrorData << "PrintRealPressure          = " << printRealPressure << std::endl;
    mirrorData << "PrintMeshVelocity          = " << printMeshVelocity << std::endl;
    mirrorData << "PrintMeshDisplacement      = " << printMeshDisplacement << std::endl;
    mirrorData << "printDistFunction          = " << printDistFunction << std::endl;
    mirrorData << "printGlueZone              = " << printGlueZone << std::endl;
    mirrorData << "printEnergyWeightFunction  = " << printEnergyWeightFunction << std::endl;
    mirrorData << "printLagrangeMultipliers   = " << printLagrangeMultipliers << std::endl;
    mirrorData << "printNodalCorrespondence   = " << printNodalCorrespondence << std::endl;
    mirrorData << "PrintProcess               = " << printProcess << std::endl << std::endl;

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++READING FEM MESH++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    getline(file, line); getline(file, line); getline(file, line);

    // Read and print out number of FEM elements and FEM nodes
    file >> numFemElem >> numFemNodes; 

	mirrorData << "FEM nodes " << numFemNodes << std::endl;
    mirrorData << "FEM elements " << numFemElem << std::endl;

    nodes_.reserve(numFemNodes);
    elements_.reserve(numFemElem);

    getline(file, line); getline(file, line); getline(file, line);
    getline(file, line); getline(file, line);

    // Read and print out nodes from the mesh
    int index = 0;
    for (int i = 0; i < numFemNodes; i++) {        
        
        double wei, x3;
        double x[dim];
        
        file >> x[0] >> x[1] >> x3 >> wei;
        std::getline(file, line);
        
        Node *node = new Node(x, index++, wei); 
        nodes_.push_back(node);
    };

    mirrorData << "Nodal Coordinates " << std::endl;
    
    for (int i = 0 ; i<numFemNodes; i++){       
        
        double *x = nodes_[i]->getCoordinates();       
        double x_[dim];
        for (int j = 0; j < dim; j++) x_[j] = x[j];
        for (int j = 0; j < dim; j++) mirrorData << x_[j] << " ";
        mirrorData << std::endl;
        
        //Setting the initial nodal values variables
        double u[dim];
        for (int j = 0; j < dim; j++) u[j] = 0.;
        nodes_[i] -> setVelocity(velocityInf);
        nodes_[i] -> setPreviousVelocity(velocityInf);
        nodes_[i] -> setPreviousMeshVelocity(u);
        nodes_[i] -> setMeshVelocity(u);
        nodes_[i] -> setPreviousCoordinates(x_);

    };


    std::getline(file, line); std::getline(file, line); std::getline(file, line);
    std::getline(file, line);
    

    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++ELEMENTS++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //Read connectivity of the elements
    index = 0;
    for (int i = 0; i < numFemElem; i++){
        
        int *connect;
        connect = new int[6];
        file >> connect[0] >> connect[1] >> connect[2] >> connect[3] >> 
        connect[4] >> connect[5];
            
        std::getline(file, line);

        connect[0] -= 1;connect[1] -= 1;connect[2] -= 1;connect[3] -= 1;
        connect[4] -= 1; connect[5] -= 1;

        Elements *el = new Elements(index++,connect,nodes_,0,fluidParameters,IsoPar_,1000); 
        
        elements_.push_back(el);
        //0 = Element finit type of element 
        //1000 = path number for FEM mesh
        for (int k = 0; k < 6; k++){
            nodes_[connect[k]] -> pushInverseIncidence(index);
        };
    };

    getline(file, line); getline(file, line); getline(file, line); 
    getline(file, line);
   
    // Read and print out number of FEM boundaries and FEM elements boundaries
    file >> numFemBoundaries >> numFemBoundElem;
    getline(file, line);

    mirrorData << "Number of Boundary groups (FEM MESH) = " << numFemBoundaries << std::endl;
    mirrorData << "Number of Boundary Elements (FEM MESH)= " << numFemBoundElem << std::endl;

    boundary_.reserve(numFemBoundElem);
   
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // +++++++++++++++++++++++++++++BONDARIES++++++++++++++++++++++++++++++++++
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //Read boundary information: connectivity and prescribed value
    index = 0;
    for (int ibound=0; ibound< numFemBoundaries; ibound++){
        
        getline(file,line);getline(file,line);getline(file,line);
        getline(file,line);getline(file,line);
        int constrain[dim+1];
        double value[dim+1];
        int numGroup;

        file >> numGroup >> constrain[0] >> constrain[1] >> constrain[2] 
        >> value[0] >> value[1] >> value[2];

        getline(file,line);getline(file,line);getline(file,line);
        getline(file,line);getline(file,line);
       
        for (int i = 0; i < numGroup; i++){        	
        	int *connectB;
        	connectB = new int[3];
            file >> connectB[0] >> connectB[1] >> connectB[2];
            getline(file,line);
            connectB[0] -= 1;
            connectB[1] -= 1;
            connectB[2] -= 1;
            Boundaries *bound = new Boundaries(connectB, index++,  
                                               constrain, value, ibound);
            boundary_.push_back(bound);
        };
    };     
     	

    //Sets boundary constrains
    for (int ibound = 0; ibound < numFemBoundElem; ibound++){
        int *connectB = boundary_[ibound] -> getBoundaryConnectivity();
        int no1 = connectB[0];
        int no2 = connectB[1];
        int no3 = connectB[2];

        for (int i = 0; i < dim; i++){
           
            if ((boundary_[ibound] -> getConstrain(i) == 1) || 
                (boundary_[ibound] -> getConstrain(i) == 3)){
                nodes_[no1] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no2] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no3] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
            };
        };
          
    };

    //Print nodal constrains
    for (int i=0; i< numFemNodes; i++){
        mirrorData<< "Constrains " << i
                  << " " << nodes_[i] -> getConstrains(0) 
                  << " " << nodes_[i] -> getConstrainValue(0) 
                  << " " << nodes_[i] -> getConstrains(1) 
                  << " " << nodes_[i] -> getConstrainValue(1) << std::endl;
    }; 

    
    //Sets fluid elements and sides on interface boundaries
    int boundElem[numFemBoundElem]; 
    for (int i = 0; i< numFemBoundElem; i++){
            if ((boundary_[i] -> getConstrain(0) > 0) ||
            (boundary_[i] -> getConstrain(1) > 0)) {
            int *connectB= boundary_[i] -> getBoundaryConnectivity();
            for (int j = 0; j < numFemElem; j++){               
                int *connect = elements_[j] -> getConnectivity();              
                int flag = 0;
                int side[3];
                for (int k = 0; k < 6; k++){
                    if ((connectB[0] == connect[k]) || 
                        (connectB[1] == connect[k]) ||
                        (connectB[2] == connect[k])){
                        side[flag] = k;
                        flag++;
                    };
                };
                if (flag == 3){
                    boundary_[i] -> setElement(elements_[j] -> getIndex());
                    //Sets element index and side
                    if ((side[0]==4) || (side[1]==4) || (side[2]==4)){
                        boundary_[i] -> setElementSide(0);
                        boundElem[i] = j;
                    };
                    if ((side[0]==5) || (side[1]==5) || (side[2]==5)){
                        boundary_[i] -> setElementSide(1);
                        boundElem[i] = j;
                    };
                    if ((side[0]==3) || (side[1]==3) || (side[2]==3)){
                        boundary_[i] -> setElementSide(2);
                        boundElem[i] = j;
                    };   
                };
            }; 
        }; 
    }; 

    std::vector<int> elemSide;
    std::vector<int> elemBound;

    for (int i =0; i <numFemElem; i++){
        for (int j =0; j<numFemBoundElem; j++){
            if (boundElem[j] == i){
                int side = boundary_[j] -> getElementSide();
                elemSide.push_back(side);
                elemBound.push_back(j);
            };
        };
        elements_[i] -> setElemSideInBoundary(elemSide,elemBound);
    };


    domainDecompositionMETIS_FEM(elements_,numFemElem,numFemNodes);

    //  //Closing the file
    //  file.close();
    //  if (deleteFiles)
    //     system((remove2 + inputFile).c_str());

    return;
};

template<>
void FluidData<2>::dataReading_ISO(const std::string& inputFile,const std::string& inputMeshIso,
                                  const std::string& mirror,const bool& deleteFiles){

    //Fluid data
    double pressInf;       			     //Undisturbed pressure 
    double rhoInf;         			     //Density
    double viscInf;        				 //Viscosity
    double fieldForces[3]; 				 //Field forces (constant) 
    double velocityInf[3];               //Infinity velocity
    int numTimeSteps;                    //Number of Time Steps
    double dTime;                        //Time Step
    int printFreq;         				 //Printing frequence of output files
    //Arlequin Problem
    double arlequinK1;                   //L2
    double arlequinK2;                   //H1
    double glueZoneThickness;			 //Thickness from gluing zone
    double arlequinEpsilon;				 //Constant > 0
    //Time integration parameters
    double integScheme;                  //Time Integration Scheme (0 - max. dissipation; 1 - no dissipation)
    double alpha_f;						 //Alphageneralized parameter - matrix except mass
    double alpha_m;						 //Alphageneralized parameter - mass
    double gamma;						 //Alphageneralized parameter 
    //Mesh data
    int numPatches;                      //Number of IGA patches
    int numCP;                           //Number of control points
    int numIsoElem;                      //Number of IGA elements 
    int numBoundariesIso;                //Number of fluid boundaries (IGA MESH)
    int numBoundElemIso;                 //Number of elements in fluid boundaries (IGA MESH)

    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);      
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    
    int dim = 2;
       
    if (rank == 0)std::cout << "Reading fluid data from \"" 
                            << inputFile << "\"" << std::endl;

   
    //Defines input and output files
    std::ifstream inputData(inputFile.c_str());
    std::ifstream file2(inputMeshIso);
    std::ofstream mirrorData(mirror.c_str());
    std::string line;
    std::string remove2 = "rm ";

    //Start to read the data problems  
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
    //Read and print out time steps and printing frequence
    inputData >> numTimeSteps >> printFreq;

    mirrorData << "Number of Time Steps   = " << numTimeSteps << std::endl;
    mirrorData << "Printing Frequence     = " << printFreq << std::endl;
    
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read and print out data from Arlequin problem
    inputData >> glueZoneThickness;
    mirrorData << "Thickness of glue zone   = " << glueZoneThickness << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    inputData >> arlequinEpsilon;
    mirrorData << "Epsilon from Arlequin formulation   = " << arlequinEpsilon << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    inputData >> arlequinK1 >> arlequinK2;
    mirrorData << "Arlequin k1   = " << arlequinK1 << std::endl;
    mirrorData << "Arlequin k2   = " << arlequinK2 << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
    //Read and print out undisturbed velocity and pressure components
    inputData >> velocityInf[0] >> velocityInf[1] >> velocityInf[2] >> pressInf;
    
    mirrorData << "Undisturbed Velocity x = " << velocityInf[0] << std::endl;
    mirrorData << "Undisturbed Velocity y = " << velocityInf[1] << std::endl;
    mirrorData << "Undisturbed Velocity z = " << velocityInf[2] << std::endl;
    mirrorData << "Undisturbed Pressure   = " << pressInf << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read and print out undisturbed density and viscosity
    inputData >> rhoInf >> viscInf;
    
    mirrorData << "Undisturbed Density    = " << rhoInf << std::endl;
    mirrorData << "Undisturbed Viscosity  = " << viscInf << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read and print out time step and time integration scheme
    inputData >> dTime >> integScheme;
    
    mirrorData << "Time Step              = " << dTime << std::endl;
    mirrorData << "Time Integration Scheme= " << integScheme << std::endl;

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);getline(inputData,line);

    //Read and print out field forces
    inputData >> fieldForces[0] >> fieldForces[1] >> fieldForces[2];
    
    mirrorData << "Field Forces x         = " << fieldForces[0] << std::endl;
    mirrorData << "Field Forces y         = " << fieldForces[1] << std::endl;
    mirrorData << "Field Forces z         = " << fieldForces[2] << std::endl 
               << std::endl;
    
    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    fluidParameters.setViscosity(viscInf);
    fluidParameters.setDensity(rhoInf);
    fluidParameters.setTimeStep(dTime);
    fluidParameters.setSpectralRadius(integScheme);
    fluidParameters.setFieldForce(fieldForces);
    fluidParameters.setArlequink1(arlequinK1);
    fluidParameters.setArlequink2(arlequinK2);

    //Read and print out Drag and lift coeficients
    inputData >> computeDragAndLift >> numberOfLines; 

    dragAndLiftBoundary.reserve(numberOfLines);
    for (int i = 0; i < numberOfLines; ++i)
    {
        int aux;
        inputData >> aux; 
        dragAndLiftBoundary.push_back(aux);
    }


    mirrorData << "Compute Drag and Lift  = " << computeDragAndLift<< std::endl;
    mirrorData << "Number of Lines  = " << numberOfLines << std::endl;
    for (int i = 0; i < numberOfLines; ++i)
    {
        mirrorData << "Lines  = " << dragAndLiftBoundary[i] << std::endl;
    }

    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    getline(inputData,line);getline(inputData,line);

    //Read and print out Printing results
    inputData >> printVelocity;              getline(inputData,line);
    inputData >> printRealVelocity;          getline(inputData,line);
    inputData >> printPressure;              getline(inputData,line);
    inputData >> printRealPressure;          getline(inputData,line);
    inputData >> printMeshVelocity;          getline(inputData,line);
    inputData >> printMeshDisplacement;      getline(inputData,line);
    inputData >> printDistFunction;          getline(inputData,line);
    inputData >> printGlueZone;              getline(inputData,line);
    inputData >> printEnergyWeightFunction;  getline(inputData,line);
    inputData >> printLagrangeMultipliers;   getline(inputData,line);
    inputData >> printNodalCorrespondence;   getline(inputData,line);
    inputData >> printProcess;          

    mirrorData << "PrintVelocity              = " << printVelocity << std::endl;
    mirrorData << "PrintRealVelocity          = " << printRealVelocity << std::endl;
    mirrorData << "PrintPressure              = " << printPressure << std::endl;
    mirrorData << "PrintRealPressure          = " << printRealPressure << std::endl;
    mirrorData << "PrintMeshVelocity          = " << printMeshVelocity << std::endl;
    mirrorData << "PrintMeshDisplacement      = " << printMeshDisplacement << std::endl;
    mirrorData << "printDistFunction          = " << printDistFunction << std::endl;
    mirrorData << "printGlueZone              = " << printGlueZone << std::endl;
    mirrorData << "printEnergyWeightFunction  = " << printEnergyWeightFunction << std::endl;
    mirrorData << "printLagrangeMultipliers   = " << printLagrangeMultipliers << std::endl;
    mirrorData << "printNodalCorrespondence   = " << printNodalCorrespondence << std::endl;
    mirrorData << "PrintProcess               = " << printProcess << std::endl << std::endl;

   
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++READING ISOGEOMETRIC MESH+++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    getline(file2, line); getline(file2, line); getline(file2, line);

    // Read and print out number of patches 
    file2 >> numPatches;

    IsoPar_.reserve(numPatches);

    getline(file2, line); getline(file2, line); getline(file2, line);
    getline(file2, line); getline(file2, line);

    // Read and print out number of control points and patches
    file2 >> numCP;    
    
    mirrorData << "Total number of Control Points " << numCP << std::endl;
    mirrorData << "Number of patches " << numPatches << std::endl;

    nodes_.reserve(numCP);

    getline(file2, line); getline(file2, line); getline(file2, line);
    getline(file2, line); getline(file2, line);
   
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++CONTROL POINTS AND KNOT VECTOR+++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   	
   	//Read and print out data from each patch:
   	//number of control points, control points/weight and data from knot vectors
    int index = 0;
    int numCPlocal; 
    double x3,wei;
    for (int ipatch = 0; ipatch < numPatches; ipatch ++){

        file2 >> numCPlocal;

        getline(file2, line); getline(file2, line); getline(file2, line);
        getline(file2, line); getline(file2, line);
            
        for (int i = 0; i < numCPlocal; i++){
            double x[dim];
            file2 >> x[0] >> x[1] >> x3 >> wei;
            std::getline(file2, line);
            Node *node = new Node(x,index++,wei);
            nodes_.push_back(node);
        };
      
        mirrorData << "Data in the " << ipatch << " PATCH" << std::endl;
        mirrorData << "Number of control points " << numCPlocal << std::endl;
        
        for (int i = (index - numCPlocal) ; i<index; i++){
           
            double *x = nodes_[i]->getCoordinates();     
            double x_[dim];
            for (int j = 0; j < dim; j++) x_[j] = x[j]; 
            for (int j = 0; j < dim; j++) mirrorData << x_[j] << " ";
            mirrorData << nodes_[i]-> getWeightPC() << std::endl;
            
            //Setting the initial nodal values variables
            double u[dim];
            for (int j = 0; j < dim; j++) u[j] = 0.;
            nodes_[i] -> setVelocity(velocityInf);
            nodes_[i] -> setPreviousVelocity(velocityInf);
            nodes_[i] -> setPreviousMeshVelocity(u);
            nodes_[i] -> setMeshVelocity(u);
            nodes_[i] -> setPreviousCoordinates(x_);
            
        };

        getline(file2, line); getline(file2, line); getline(file2, line);
        getline(file2, line); 

        int ncp[dim],deg[dim];
        int ncp2,deg2;

        file2 >> ncp[0] >> ncp[1] >> ncp2;

        getline(file2, line); getline(file2, line); getline(file2, line); 
        getline(file2, line); getline(file2, line);

        file2 >> deg[0] >> deg[1] >> deg2;

        IsoParameters *isopar = new IsoParameters(ipatch,deg,ncp);
        IsoPar_.push_back(isopar);

        int dim_u = ncp[0] + deg[0] + 1;
        int dim_v = ncp[1] + deg[1] + 1;

        double *uknot_;
        uknot_ = new double[dim_u];
        double *vknot_;
        vknot_ = new double[dim_v];

        getline(file2, line); getline(file2, line); getline(file2, line); 
        getline(file2, line); getline(file2, line);

        for (int i = 0; i < dim_u; i++) {
            file2 >> uknot_[i];
            getline(file2,line);
        }

        getline(file2, line); getline(file2, line); getline(file2, line); 
        getline(file2, line); 

        for (int j = 0; j < dim_v; j++) {
            file2 >> vknot_[j];
            getline(file2,line);
        }

        mirrorData << " Number of Control Points x: " << ncp[0] 
                   << " Function Degree: " << deg [0] << std::endl;
        mirrorData << " Number of Control Points y: " << ncp[1] 
                   << " Function Degree: " << deg [1] << std::endl;
        
        mirrorData << "uKnot " << std::endl;
        double uknotP[dim_u];
        for (int i=0;i <dim_u; i++) uknotP[i] = uknot_[i];
        for (int i = 0; i< dim_u; i++) mirrorData << uknotP[i] << std::endl;

        double vknotP[dim_v];
        for (int i=0;i <dim_v; i++) vknotP[i] = vknot_[i];
        mirrorData << "vKnot " << std::endl;
        for (int i = 0; i< dim_v; i++) mirrorData << vknotP[i] << std::endl;

        IsoPar_[ipatch] -> setuKnot(uknot_);
        IsoPar_[ipatch] -> setvKnot(vknot_);

        getline(file2, line); getline(file2, line); getline(file2, line); 
        getline(file2, line); getline(file2, line); getline(file2, line); 
        getline(file2, line); getline(file2, line); 


    }; //iPatch


    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++ELEMENTS++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //Read the element data from each patch

    int ngp = 0; // the index of control points
    numIsoElem = 0; // counts the number of elements
    for (int ipatch = 0; ipatch < numPatches; ipatch++){
        
        int ncp[dim],deg[dim]; 
        deg[0] = IsoPar_[ipatch] -> getDegree(0);
        deg[1] = IsoPar_[ipatch] -> getDegree(1);
        ncp[0] = IsoPar_[ipatch] -> getNcp(0);
        ncp[1] = IsoPar_[ipatch] -> getNcp(1);
        double *uknot= IsoPar_[ipatch] -> getuKnot();
        double *vknot= IsoPar_[ipatch] -> getvKnot();     
        
        for (int j = 0; j < ncp[1]; j++){   
            for (int i = 0; i < ncp[0]; i++ ) {
                int INC[2];
                INC[0] = i;
                INC[1] = j;
                nodes_[ngp++] -> setINC(INC);                      
                // finding the elements 
                if ((i >= deg[0]) && (j >= deg[1])) {
                    if (((uknot[i])!= (uknot[i+1])) && ((vknot[j])!= (vknot[j+1]))) {
                        numIsoElem ++;
                    };
                } ;
            };
        };
    };

    elements_.reserve(numIsoElem);

    //number of elements in each patch direction
    int numElemPatchx;
    int numElemPatchy;  

    //direction x
    for (int ipatch = 0; ipatch < numPatches; ipatch++){
        
        numElemPatchx = 0;
        int ncp,deg; 
        deg = IsoPar_[ipatch] -> getDegree(0);
        ncp = IsoPar_[ipatch] -> getNcp(0);
        double *uknot= IsoPar_[ipatch] -> getuKnot();
        
        for (int i = 0; i < ncp; i++ ) {
            //finding elements in each direction
            if (i >= deg) {
                if (uknot[i]!= uknot[i+1]) {
                    numElemPatchx++;
                };
            };
        };
        IsoPar_[ipatch] -> setNumElemPatchx(numElemPatchx);
    };

    //direction y
    for (int ipatch = 0; ipatch < numPatches; ipatch++){
        
        numElemPatchy = 0;
        int ncp,deg; 
        deg = IsoPar_[ipatch] -> getDegree(1);
        ncp = IsoPar_[ipatch] -> getNcp(1);
        double *vknot = IsoPar_[ipatch] -> getvKnot();
        
        for (int j = 0; j < ncp; j++ ) {
            //finding elements in each direction
            if (j >= deg) {
                if (vknot[j]!= vknot[j+1]) {
                    numElemPatchy++;
                };
            };
        };
        IsoPar_[ipatch] -> setNumElemPatchy(numElemPatchy);
    };


    int ElemType; // 0 - FEM elementos, 1 - IA elementos   
    ngp = 0;     // the index of control points
    index = 0;   // the index of elements

    for (int ipatch = 0; ipatch < numPatches; ipatch++){
        
        int ncp[dim],deg[dim]; 
        deg[0] = IsoPar_[ipatch] -> getDegree(0);
        deg[1] = IsoPar_[ipatch] -> getDegree(1);
        ncp[0] = IsoPar_[ipatch] -> getNcp(0);
        ncp[1] = IsoPar_[ipatch] -> getNcp(1);
        double *uknot = IsoPar_[ipatch] -> getuKnot();
        double *vknot = IsoPar_[ipatch] -> getvKnot();
        
        for (int j = 0; j < ncp[1]; j++){
            for (int i = 0; i < ncp[0]; i++ ) {
                ngp++;                      
                // finding the elements 
                if ((i >= deg[0]) && (j >= deg[1])) {
                    if (((uknot[i])!= (uknot[i+1])) && ((vknot[j])!= (vknot[j+1]))) {
                        int *connect;
                        connect = new int[9];
                        for (int jloc = 0; jloc <= deg[1] ; jloc++){
                            for (int iloc = 0; iloc <=  deg[0]; iloc++){
                                // global function number
                                int ngf = ngp - jloc*ncp[0] - iloc - 1; //
                                // local number function
                                int nlf = 8 - jloc*(deg[0] + 1) - iloc;
                                connect[nlf] = ngf;
                            }
                        }               
                        Elements *el = new Elements(index++,connect,nodes_,1,
                                                    fluidParameters,IsoPar_,ipatch); 
                        //1 - it's the element type to IA elements.
                        elements_.push_back(el);
                        for (int k = 0; k <9; k++){
                            nodes_[connect[k]] -> pushInverseIncidence(index);
                        }; 
                    };
                };    
            };
        };
    };


    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //+++++++++++++++++++++++++++++BONDARIES++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   
    
    //Read and print out number of boundary groups and number of elements by boundary
    file2 >> numBoundariesIso >> numBoundElemIso;                   

    mirrorData << "Number of Boundary groups (IA MESH) = " 
               << numBoundariesIso << std::endl;
    mirrorData << "Number of Boundary Elements (IA MESH)= " 
               << numBoundElemIso << std::endl;
    
    boundary_.reserve(numBoundElemIso);
  
    //Read each boundary information  
    index = 0;
    for (int ibound=0; ibound< numBoundariesIso; ibound++){

        getline(file2,line);getline(file2,line);getline(file2,line);
        getline(file2,line);getline(file2,line);getline(file2,line);

        int constrain[dim+1];
        double value[dim+1];
        int numGroup;
                
        file2 >> numGroup >> constrain[0] >> constrain[1] >> constrain[2] 
                  >> value[0] >> value[1] >> value[2];

        getline(file2,line);getline(file2,line);getline(file2,line);
        getline(file2,line);getline(file2,line);

        //Icpbc contain:
        //2 indexes of the control points with boundary conditions 
        //direction with boundary conditions (x = 0, y = 1)
        //belonging patch
        int Icpbc[numGroup][dim+2]; 
        for (int i = 0; i < numGroup; i++){
          file2 >> Icpbc[i][0] >> Icpbc[i][1] >> Icpbc[i][2] >> Icpbc[i][3];    
        };
        
        //any of the control points belong to the same patch
        int ncp[dim],deg[dim]; 
        deg[0] = IsoPar_[Icpbc[0][3]] -> getDegree(0);
        deg[1] = IsoPar_[Icpbc[0][3]] -> getDegree(1);
        ncp[0] = IsoPar_[Icpbc[0][3]] -> getNcp(0);
        ncp[1] = IsoPar_[Icpbc[0][3]] -> getNcp(1);
        double *uknot = IsoPar_[Icpbc[0][3]] -> getuKnot();
        double *vknot = IsoPar_[Icpbc[0][3]] -> getvKnot();

        int contcp = 0;
        int ncp0,ncp1;
        for (int aux = 0; aux < Icpbc[0][3]; aux++){
            ncp0 = IsoPar_[aux] -> getNcp(0);
            ncp1 = IsoPar_[aux] -> getNcp(1);
            contcp += ncp0 * ncp1;
        };

        //Finding the boundary connectivity
        for (int k = 0; k < numGroup; k++){
            int *connectB;
            connectB = new int[3];
            int i = Icpbc[k][0];
            int j = Icpbc[k][1];
            if (Icpbc[k][2] == 0){
                if ((uknot[i] != uknot[i+1]) && (i >= (Icpbc[0][0] + deg[0]))) { 
                    for (int l = 0; l <= deg[0]; l++){
                        connectB[l] = contcp + i + ncp[0]*j + l - deg[0];
                    };
                    Boundaries *bound = new Boundaries(connectB, index++,constrain, value, ibound);
                    boundary_.push_back(bound);
               };
            }; 
            if (Icpbc[k][2] == 1){
                if ((vknot[j] != vknot[j+1]) && (j >= (Icpbc[0][1] + deg[1]))) {
                    for (int l = 0; l <= deg[1]; l++){
                        connectB[l] = contcp + i + ncp[0]*j + l*ncp[0] - deg[0]*ncp[0];
                    };
                    Boundaries *bound = new Boundaries(connectB, index++,constrain, value, ibound);
                    boundary_.push_back(bound);
                };
            };                      
        };
    };  

    
    //Sets nodal boundary constrains
    numBoundElemIso = index;
    for (int ibound = 0; ibound < numBoundElemIso; ibound++){        
        int *connectB = boundary_[ibound] -> getBoundaryConnectivity();
        int no1 = connectB[0];
        int no2 = connectB[1];
        int no3 = connectB[2];
        
        for (int i = 0; i < dim; i++){
           
            if ((boundary_[ibound] -> getConstrain(i) == 1) || 
                (boundary_[ibound] -> getConstrain(i) == 3)){
                nodes_[no1] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no2] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no3] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
            };
        }; 
    };

    //Print nodal constrains
    for (int i=0; i< numCP; i++){
        mirrorData<< "Constrains " << i
                  << " " << nodes_[i] -> getConstrains(0) 
                  << " " << nodes_[i] -> getConstrainValue(0) 
                  << " " << nodes_[i] -> getConstrains(1) 
                  << " " << nodes_[i] -> getConstrainValue(1) << std::endl;
    }; 


    //Sets fluid elements and sides on interface boundaries
    int boundElem[numBoundElemIso]; 
    for (int i = 0; i< numBoundElemIso; i++){
        if ((boundary_[i] -> getConstrain(0) > 0) ||
            (boundary_[i] -> getConstrain(1) > 0)) {            
            int *connectB = boundary_[i] -> getBoundaryConnectivity();            
            for (int j = 0; j < numIsoElem; j++){
                int *connect = elements_[j] -> getConnectivity();              
                int flag = 0;
                int side[3];
                for (int k = 0; k < 9; k++){
                    if ((connectB[0] == connect[k]) || 
                        (connectB[1] == connect[k]) ||
                        (connectB[2] == connect[k])){
                        side[flag] = k;
                        flag++;
                    };
                };
                if (flag == 3){
                    boundary_[i] -> setElement(elements_[j] -> getIndex());
                    if ((side[0]==1) || (side[1]==1) || (side[2]==1)){
                        boundary_[i] -> setElementSide(0);
                        boundElem[i] = j;
                    };
                    if ((side[0]==5) || (side[1]==5) || (side[2]==5)){
                        boundary_[i] -> setElementSide(1);
                        boundElem[i] = j;
                    };
                    if ((side[0]==7) || (side[1]==7) || (side[2]==7)){
                        boundary_[i] -> setElementSide(2);
                        boundElem[i] = j;
                    }; 
                    if ((side[0]==3) || (side[1]==3) || (side[2]==3)){
                        boundary_[i] -> setElementSide(3);
                        boundElem[i] = j;
                    };  
                }; 
            }; 
        }; 
    }; 

    std::vector<int> elemSide;
    std::vector<int> elemBound;
    elemSide.clear();
    elemBound.clear();
    for (int i =0; i <numIsoElem; i++){
        for (int j =0; j<numBoundElemIso; j++){ 
            if (boundElem[j] == i){
                int side = boundary_[j] -> getElementSide();
                elemSide.push_back(side);
                elemBound.push_back(j);
            };
        };
        elements_[i] -> setElemSideInBoundary(elemSide,elemBound);
        elemSide.clear();
        elemBound.clear();
    };


    // Storing the matching control points of each control point
    double distol = 1.e-6;
    for (int i = 0; i < numCP; i++){
        std::vector<int> matchingCP_;
        for (int j = 0; j < numCP; j++){
            if (i != j){
                double *xi = nodes_[i]->getCoordinates(); 
                double *xj = nodes_[j]->getCoordinates();  
                double dist2 = ((xi[0]- xj[0]) * (xi[0]- xj[0])) + ((xi[1]- xj[1]) * (xi[1]- xj[1]));
                if (dist2 <= distol) {
                    matchingCP_.push_back(j);
                };
            };
        };
        int sizeMcp = matchingCP_.size();
        int *matchingCP;
        matchingCP = new int[sizeMcp];
        for (int i = 0; i<sizeMcp; i++) matchingCP[i] = matchingCP_[i];
        nodes_[i]->setMatchingCP(matchingCP);
        nodes_[i]->setSizeMcp(sizeMcp);
        matchingCP_.clear();   
    };

    // Creates a new numeration of control points (descarting the repeated ones)
    NCNumberNodes = 0; //non coincidente number nodes;
    int flag[numCP] = {};
    for (int inodes = 0; inodes < numCP; inodes++) {        
        int *matchingCP_ = nodes_[inodes]-> getMatchingCP();
        int sizemcp = nodes_[inodes]-> getSizeMcp();     
        if (sizemcp == 0) {
            nodes_[inodes]-> setnewcon(NCNumberNodes++);
        } else {
            if (flag[inodes] == 0){
                flag[inodes] = 1;
                nodes_[inodes]-> setnewcon(NCNumberNodes);
                for (int i = 0; i < sizemcp; i++) {
                    if (flag[matchingCP_[i]] == 0){
                        nodes_[matchingCP_[i]] -> setnewcon(NCNumberNodes);
                        flag[matchingCP_[i]] = 1;
                    }; 
                };
                NCNumberNodes++;
            };
        };
    };

    BezierConnectivity(numPatches);

    domainDecompositionMETIS_ISO(elements_,numIsoElem,numCP);

    //  //Closing the file
    //  file2.close();
    //  if (deleteFiles)
    //     system((remove2 + inputFile).c_str());

    return;
};

//----------------------------------------------------------------------------------
//----------------------------READS FLUID INPUT FILE--------------------------------
//----------------------------------------------------------------------------------
// template<>
// void FluidData<3>::dataReading_FEM(const std::string& inputFile,const std::string& inputMesh, 
//                            const std::string& mirror,const bool& deleteFiles){

 //    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);      
 //    MPI_Comm_size(PETSC_COMM_WORLD, &size);
    
 //    int dimension = 3;
       
 //    if (rank == 0)std::cout << "Reading fluid data from \"" 
 //                            << inputFile << "\"" << std::endl;

 //    //Defines input and output files
 //    std::ifstream inputData(inputFile.c_str());
 //    std::ifstream file(inputMesh);
 //    std::ofstream mirrorData(mirror.c_str());
 //    std::string line;
	// std::string remove2 = "rm ";
      
 //    //Start to read the data problems  
 //    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
 //    //Read and print out time steps and printing frequence
 //    inputData >> numTimeSteps >> printFreq;
    
 //    mirrorData << "Number of Time Steps   = " << numTimeSteps << std::endl;
 //    mirrorData << "Printing Frequence     = " << printFreq << std::endl;
    
 //    getline(inputData,line);getline(inputData,line);getline(inputData,line);
 //    getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
 //    //Read and print out undisturbed velocity and pressure components
 //    inputData >> velocityInf[0] >> velocityInf[1] >> velocityInf[2] >> pressInf;
    
 //    mirrorData << "Undisturbed Velocity x = " << velocityInf[0] << std::endl;
 //    mirrorData << "Undisturbed Velocity y = " << velocityInf[1] << std::endl;
 //    mirrorData << "Undisturbed Velocity z = " << velocityInf[2] << std::endl;
 //    mirrorData << "Undisturbed Pressure   = " << pressInf << std::endl;

 //    getline(inputData,line);getline(inputData,line);getline(inputData,line);
 //    getline(inputData,line);getline(inputData,line);

 //    //Read and print out undisturbed density and viscosity
 //    inputData >> rhoInf >> viscInf;
    
 //    mirrorData << "Undisturbed Density    = " << rhoInf << std::endl;
 //    mirrorData << "Undisturbed Viscosity  = " << viscInf << std::endl;

 //    getline(inputData,line);getline(inputData,line);getline(inputData,line);
 //    getline(inputData,line);getline(inputData,line);

 //    //Read and print out time step and time integration scheme
 //    inputData >> dTime >> integScheme;
    
 //    mirrorData << "Time Step              = " << dTime << std::endl;
 //    mirrorData << "Time Integration Scheme= " << integScheme << std::endl;

 //    alpha_f = 1. / (1. + integScheme);
 //    alpha_m = 0.5 * (3. - integScheme) / (1. + integScheme);
 //    gamma = 0.5 + alpha_m - alpha_f;

 //    // alpha_f = 1.;
 //    // alpha_m = 1.;
 //    // gamma = 1.;

 //    getline(inputData,line);getline(inputData,line);getline(inputData,line);
 //    getline(inputData,line);getline(inputData,line);getline(inputData,line);

 //    //Read and print out field forces
 //    inputData >> fieldForces[0] >> fieldForces[1] >> fieldForces[2];
    
 //    mirrorData << "Field Forces x         = " << fieldForces[0] << std::endl;
 //    mirrorData << "Field Forces y         = " << fieldForces[1] << std::endl;
 //    mirrorData << "Field Forces z         = " << fieldForces[2] << std::endl \
 //               << std::endl;
    
 //    getline(inputData,line);getline(inputData,line);getline(inputData,line);
 //    getline(inputData,line);getline(inputData,line);

 //    fluidParameters.setViscosity(viscInf);
 //    fluidParameters.setDensity(rhoInf);
 //    fluidParameters.setTimeStep(dTime);
 //    fluidParameters.setSpectralRadius(integScheme);
 //    fluidParameters.setFieldForce(fieldForces);

 //    //Read and print out Drag and lift coeficients
 //    inputData >> computeDragAndLift >> numberOfLines; 

 //    dragAndLiftBoundary.reserve(numberOfLines);
 //    for (int i = 0; i < numberOfLines; ++i)
 //    {
 //        int aux;
 //        inputData >> aux; 
 //        dragAndLiftBoundary.push_back(aux);
 //    }

 //    mirrorData << "Compute Drag and Lift  = " << computeDragAndLift<< std::endl;
 //    mirrorData << "Number of Lines  = " << numberOfLines << std::endl;
 //    for (int i = 0; i < numberOfLines; ++i)
 //    {
 //        mirrorData << "Lines  = " << dragAndLiftBoundary[i] << std::endl;
 //    }

 //    getline(inputData,line);getline(inputData,line);getline(inputData,line);
 //    getline(inputData,line);getline(inputData,line);

 //    //Read and print out Printing results
 //    inputData >> printVelocity;              getline(inputData,line);
 //    inputData >> printPressure;              getline(inputData,line);
 //    inputData >> printVorticity;             getline(inputData,line);
 //    inputData >> printMeshVelocity;          getline(inputData,line);
 //    inputData >> printMeshDisplacement;      getline(inputData,line);
 //    inputData >> printJacobian;              getline(inputData,line);
 //    inputData >> printProcess;              

 //    mirrorData << "PrintVelocity              = " << printVelocity << std::endl;
 //    mirrorData << "PrintPressure              = " << printPressure << std::endl;
 //    mirrorData << "PrintVorticity             = " << printVorticity 
 //               << std::endl;
 //    mirrorData << "PrintMeshVelocity          = " << printMeshVelocity
 //               << std::endl;
 //    mirrorData << "PrintMeshDisplacement      = " << printMeshDisplacement
 //               << std::endl;
 //    mirrorData << "PrintJacobian              = " << printJacobian << std::endl;
 //    mirrorData << "PrintProcess               = " << printProcess << std::endl 
 //               << std::endl;

   
 //    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 //    //++++++++++++++++++++++++++READIN FEM MESH+++++++++++++++++++++++++++++++
 //    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    

 //    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 //    //+++++++++++++++++++++++++++++++++NODES++++++++++++++++++++++++++++++++++
 //    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 //    //read and print out nodes coordinates  



 //    getline(file, line); getline(file, line); getline(file, line);

 //    // number of FEM elements and FEM nodes
 //    file >> numFemElem >> numFemNodes; 
 //    nodes_.reserve(numFemNodes);

 //    getline(file, line); getline(file, line); getline(file, line);
 //    getline(file, line); getline(file, line);
    
 //    int index = 0;
 //    int wei;
 //    for (int i = 0; i < numFemNodes; i++) {
 //        double x[3];
 //        file >> x[0] >> x[1] >> x[2] >> wei;
 //        std::getline(file, line);
 //        Node *node = new Node(x, index++, wei); 
 //        //1. its the weight for all nodes of FEM mesh
 //        nodes_.push_back(node);
 //    };
  
 //    mirrorData << "Nodal Coordinates " << std::endl;
 //    for (int i = 0 ; i<numFemNodes; i++){
       
 //        double *x;
 //        x = nodes_[i]->getCoordinates();   

 //        double x_[3];
 //        for (int i = 0; i<3; i++) x_[i] = x[i];

 //        for (int j=0; j<3; j++){
 //            mirrorData << x_[j] << " ";
 //        };
        
 //        mirrorData << std::endl;
        
 //        nodes_[i] -> setVelocity(velocityInf);
 //        nodes_[i] -> setPreviousVelocity(velocityInf);
        
 //        double u[3];
 //        u[0] = 0.; u[1] = 0.; u[2] = 0.;
 //        nodes_[i] -> setPreviousMeshVelocityComponent(0,0.);
 //        nodes_[i] -> setPreviousMeshVelocityComponent(1,0.);
 //        nodes_[i] -> setPreviousMeshVelocityComponent(2,0.);
 //        nodes_[i] -> setMeshVelocity(u);
            
 //        nodes_[i] -> setPreviousCoordinates(0,x_[0]);
 //        nodes_[i] -> setPreviousCoordinates(1,x_[1]);
 //        nodes_[i] -> setPreviousCoordinates(2,x_[2]);
 //    };


 //    std::vector<Elements *>   elementsAux_;
 //    elementsAux_.reserve(numFemElem);

 //    std::getline(file, line); std::getline(file, line); std::getline(file, line);
 //    std::getline(file, line);

 //    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 //    //++++++++++++++++++++++++++++++++ELEMENTS++++++++++++++++++++++++++++++++
 //    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 //    //read connectivity of the elements
 //    index = 0;
 //    for (int i = 0; i < numFemElem; i++){
        
 //        int *connect;
 //        connect = new int[10];
        
 //        file >> connect[0] >> connect[1] >> connect[2] >> connect[3] >> connect[4] >> connect[5] >>
 //                connect[6] >> connect[7] >> connect[8] >> connect[9];
            
 //        std::getline(file, line);
             
 //        connect[0] -= 1;
 //        connect[1] -= 1;
 //        connect[2] -= 1;
 //        connect[3] -= 1;
 //        connect[4] -= 1;
 //        connect[5] -= 1;
 //        connect[6] -= 1;
 //        connect[7] -= 1;
 //        connect[8] -= 1;
 //        connect[9] -= 1;
            
 //        Elements *el = new Elements(index++,connect,nodes_,0,fluidParameters,IsoPar_,1000); 
 //        elementsAux_.push_back(el);
 //        //0 - it's the element type to FEM mesh.
 //        //1000 - it's the patches number of FEM elements 

 //        for (int k = 0; k < 10; k++){
 //            nodes_[connect[k]] -> pushInverseIncidence(index);
 //        };

 //    };
    
 //    getline(file,line);getline(file,line);getline(file,line);
 //    getline(file,line);

 //    file >> numFemBoundaries >> numFemBoundElem;
 //    getline(file,line);

 //    mirrorData << "Number of Boundary groups (FEM MESH) = " << numFemBoundaries << std::endl;
 //    mirrorData << "Number of Boundary Elements (FEM MESH)= " << numFemBoundElem << std::endl;

	// boundary_.reserve(numFemBoundElem);

 //    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 //    //+++++++++++++++++++++++++++++BONDARIES++++++++++++++++++++++++++++++++++
 //    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

 //    //Read boundary information: connectivity and prescribed value
 //    index = 0;
 //    for (int ibound=0; ibound< numFemBoundaries; ibound++){
        
 //        getline(file,line);getline(file,line);getline(file,line);
 //        getline(file,line);getline(file,line);
 //        int constrain[3];
 //        double value[3];
 //        int numGroup;
                
 //        file >> numGroup >> constrain[0] >> constrain[1] >> constrain[2] \
 //                  >> value[0] >> value[1] >> value[2];

 //        getline(file,line);getline(file,line);getline(file,line);
 //        getline(file,line);getline(file,line);
       
 //        for (int i = 0; i < numGroup; i++){

 //        	int *connectB;
 //        	connectB = new int[6];

 //            file >> connectB[0] >> connectB[1] >> connectB[2] >> connectB[3] >> connectB[4] >> connectB[5];
 //            getline(file,line);
 //            connectB[0] -= 1;
 //            connectB[1] -= 1;
 //            connectB[2] -= 1;
 //            connectB[3] -= 1;
 //            connectB[4] -= 1;
 //            connectB[5] -= 1;
 //            Boundaries *bound = new Boundaries(connectB, index++,  \
 //                                               constrain, value, ibound);
 //            boundary_.push_back(bound);
 //        };
 //    };        

 //    //Decomposition of the mesh 
 //    domainDecompositionMETIS(elementsAux_);
    
 //    if (rank == 0){
 //        for (int i = 0; i < numFemElem; ++i) delete elementsAux_[i];
 //        elementsAux_.clear();
 //    };

 //    // Each rank reads and prints the information about its elements
 //    MPI_Barrier(PETSC_COMM_WORLD);

 //    std::string result;
 //    std::ostringstream convert;

 //    convert << rank+000;
 //    result = convert.str();
 //    std::string s = "mesh"+result+".dat";

 //    std::ifstream mesh(s.c_str(), std::ios_base::out);

 //    mesh >> numElem;

 //    elements_.reserve(numElem);
    
 //    for (int i = 0; i < numElem; i++){
        
 //    	int npatch;
 //    	int ind_;
 //    	int *connect;
 //    	int ElemType;

 //    	connect = new int[10];

 //    	mesh >>  ind_ >> ElemType >> npatch >> connect[0] >> connect[1] >> connect[2] >> connect[3] >> connect[4] 
 //                 >> connect[5] >> connect[6] >> connect[7] >> connect[8] >> connect[9];

 //        Elements *el = new Elements(ind_,connect,nodes_,ElemType,fluidParameters,IsoPar_,npatch);
 //    	elements_.push_back(el);
 //    };

 //    MPI_Barrier(PETSC_COMM_WORLD);
 //    MPI_Allreduce(&numElem,&numFemElem,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
   
 //    //Sets boundary constrains
 //    for (int ibound = 0; ibound < numFemBoundElem; ibound++){
      
 //        int *connectB;
 //        connectB = boundary_[ibound] -> getBoundaryConnectivity();

 //    	int no1 = connectB[0];
 //        int no2 = connectB[1];
 //        int no3 = connectB[2];
 //        int no4 = connectB[3];
 //        int no5 = connectB[4];
 //        int no6 = connectB[5];
      
 //        if ((boundary_[ibound] -> getConstrain(0) == 1) || (boundary_[ibound] -> getConstrain(0) == 3)){

 //            nodes_[no1] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
 //                                     boundary_[ibound] -> getConstrainValue(0));
 //            nodes_[no2] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
 //                                     boundary_[ibound] -> getConstrainValue(0));
 //            nodes_[no3] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
 //                                     boundary_[ibound] -> getConstrainValue(0));
 //            nodes_[no4] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
 //                                     boundary_[ibound] -> getConstrainValue(0));
 //            nodes_[no5] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
 //                                     boundary_[ibound] -> getConstrainValue(0));
 //            nodes_[no6] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
 //                                     boundary_[ibound] -> getConstrainValue(0));
 //        };

 //        if((boundary_[ibound] -> getConstrain(1) == 1) || (boundary_[ibound] -> getConstrain(1) == 3)){
            
 //            nodes_[no1] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
 //                                     boundary_[ibound] -> getConstrainValue(1));
 //            nodes_[no2] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
 //                                     boundary_[ibound] -> getConstrainValue(1));
 //            nodes_[no3] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
 //                                     boundary_[ibound] -> getConstrainValue(1));
 //            nodes_[no4] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
 //                                     boundary_[ibound] -> getConstrainValue(1));
 //            nodes_[no5] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
 //                                     boundary_[ibound] -> getConstrainValue(1));
 //            nodes_[no6] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
 //                                     boundary_[ibound] -> getConstrainValue(1));
 //        };   


 //        if((boundary_[ibound] -> getConstrain(2) == 1) || (boundary_[ibound] -> getConstrain(2) == 3)){
            
 //            nodes_[no1] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
 //                                     boundary_[ibound] -> getConstrainValue(2));
 //            nodes_[no2] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
 //                                     boundary_[ibound] -> getConstrainValue(2));
 //            nodes_[no3] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
 //                                     boundary_[ibound] -> getConstrainValue(2));
 //            nodes_[no4] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
 //                                     boundary_[ibound] -> getConstrainValue(2));
 //            nodes_[no5] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
 //                                     boundary_[ibound] -> getConstrainValue(2));
 //            nodes_[no6] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
 //                                     boundary_[ibound] -> getConstrainValue(2));
 //        };   

 //    };

 //    //Print nodal constrains
 //    for (int i=0; i< numFemNodes; i++){
 //        mirrorData<< "Constrains " << i
 //                  << " " << nodes_[i] -> getConstrains(0) \
 //                  << " " << nodes_[i] -> getConstrainValue(0) \
 //                  << " " << nodes_[i] -> getConstrains(1) \
 //                  << " " << nodes_[i] -> getConstrainValue(1) \
 //                  << " " << nodes_[i] -> getConstrains(2) \
 //                  << " " << nodes_[i] -> getConstrainValue(2) << std::endl;
 //    }; 

 //    //Sets fluid elements and sides on interface boundaries (just for drag and lift coeficients)
 //    //This part must be change in the overlap problem
 //    for (int i = 0; i < numFemBoundElem; i++){
 //        for (int k = 0; k < numberOfLines; k++){
 //            if (boundary_[i] -> getBoundaryGroup() == dragAndLiftBoundary[k]){
 //                int *connectB;
 //                connectB = boundary_[i] -> getBoundaryConnectivity();
 //                for (int j=0; j<numElem; j++){
                    
 //                    int *connect;
 //                    connect = elements_[j] -> getConnectivity();  
                    
 //                    int flag = 0;
 //                    int side[6];
 //                    for (int k=0; k<10; k++){
 //                        if ((connectB[0] == connect[k]) || 
 //                            (connectB[1] == connect[k]) ||
 //                            (connectB[2] == connect[k]) ||
 //                            (connectB[3] == connect[k]) ||
 //                            (connectB[4] == connect[k]) ||
 //                            (connectB[5] == connect[k])) {
 //                            side[flag] = k;
 //                            flag++;
 //                        };
 //                    };
 //                    if (flag == 6){   
 //                        boundary_[i] -> setElement(elements_[j] -> getIndex());
 //                        //Sets element index and side
 //                        if (((side[0]==1) || (side[1]==1) || (side[2]==1) || (side[3]==1) || (side[4]==1) || (side[5]==1)) &
 //                         ((side[0]==2) || (side[1]==2) || (side[2]==2) || (side[3]==2) || (side[4]==2) || (side[5]==2)) &
 //                         ((side[0]==3) || (side[1]==3) || (side[2]==3) || (side[3]==3) || (side[4]==3) || (side[5]==3))) {
 //                            boundary_[i] -> setElementSide(0);
 //                            elements_[j] -> setElemSideInBoundary(0);
 //                        };
 //                        if (((side[0]==0) || (side[1]==0) || (side[2]==0) || (side[3]==0) || (side[4]==0) || (side[5]==0)) &
 //                         ((side[0]==2) || (side[1]==2) || (side[2]==2) || (side[3]==2) || (side[4]==2) || (side[5]==2)) &
 //                         ((side[0]==3) || (side[1]==3) || (side[2]==3) || (side[3]==3) || (side[4]==3) || (side[5]==3))) {
 //                            boundary_[i] -> setElementSide(1);
 //                            elements_[j] -> setElemSideInBoundary(1);
 //                        };
 //                        if (((side[0]==0) || (side[1]==0) || (side[2]==0) || (side[3]==0) || (side[4]==0) || (side[5]==0)) &
 //                         ((side[0]==1) || (side[1]==1) || (side[2]==1) || (side[3]==1) || (side[4]==1) || (side[5]==1)) &
 //                         ((side[0]==3) || (side[1]==3) || (side[2]==3) || (side[3]==3) || (side[4]==3) || (side[5]==3))) {
 //                            boundary_[i] -> setElementSide(2);
 //                            elements_[j] -> setElemSideInBoundary(2);
 //                        };
 //                        if (((side[0]==0) || (side[1]==0) || (side[2]==0) || (side[3]==0) || (side[4]==0) || (side[5]==0)) &
 //                         ((side[0]==1) || (side[1]==1) || (side[2]==1) || (side[3]==1) || (side[4]==1) || (side[5]==1)) &
 //                         ((side[0]==2) || (side[1]==2) || (side[2]==2) || (side[3]==2) || (side[4]==2) || (side[5]==2))) {
 //                            boundary_[i] -> setElementSide(3);
 //                            elements_[j] -> setElemSideInBoundary(3);
 //                        };
 //                    };
 //                }; //elements
 //            }; //drag and lift boundary
 //        }; //number of lines          
 //    }; //number of boundary elements


 //     //Closing the file
 //     file.close();
 //     if (deleteFiles)
 //        system((remove2 + inputFile).c_str());

 //    return;
// };


template<>
void FluidData<3>::dataReading_ISO(const std::string& inputFile,const std::string& inputMeshIso,
                             const std::string& mirror,const bool& deleteFiles){

//     MPI_Comm_rank(PETSC_COMM_WORLD, &rank);      
//     MPI_Comm_size(PETSC_COMM_WORLD, &size);
    
//     int dimension = 3;
       
//     if (rank == 0)std::cout << "Reading fluid data from \"" 
//                             << inputFile << "\"" << std::endl;

//     //Defines input and output files
//     std::ifstream inputData(inputFile.c_str());
//     std::ifstream file2(inputMeshIso);
//     std::ofstream mirrorData(mirror.c_str());
//     std::string line;
	// std::string remove2 = "rm ";
      
//     //Start to read the data problems  
//     getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
//     //Read and print out time steps and printing frequence
//     inputData >> numTimeSteps >> printFreq;
    
//     mirrorData << "Number of Time Steps   = " << numTimeSteps << std::endl;
//     mirrorData << "Printing Frequence     = " << printFreq << std::endl;
    
//     getline(inputData,line);getline(inputData,line);getline(inputData,line);
//     getline(inputData,line);getline(inputData,line);getline(inputData,line);
    
//     //Read and print out undisturbed velocity and pressure components
//     inputData >> velocityInf[0] >> velocityInf[1] >> velocityInf[2] >> pressInf;
    
//     mirrorData << "Undisturbed Velocity x = " << velocityInf[0] << std::endl;
//     mirrorData << "Undisturbed Velocity y = " << velocityInf[1] << std::endl;
//     mirrorData << "Undisturbed Velocity z = " << velocityInf[2] << std::endl;
//     mirrorData << "Undisturbed Pressure   = " << pressInf << std::endl;

//     getline(inputData,line);getline(inputData,line);getline(inputData,line);
//     getline(inputData,line);getline(inputData,line);

//     //Read and print out undisturbed density and viscosity
//     inputData >> rhoInf >> viscInf;
    
//     mirrorData << "Undisturbed Density    = " << rhoInf << std::endl;
//     mirrorData << "Undisturbed Viscosity  = " << viscInf << std::endl;

//     getline(inputData,line);getline(inputData,line);getline(inputData,line);
//     getline(inputData,line);getline(inputData,line);

//     //Read and print out time step and time integration scheme
//     inputData >> dTime >> integScheme;
    
//     mirrorData << "Time Step              = " << dTime << std::endl;
//     mirrorData << "Time Integration Scheme= " << integScheme << std::endl;

//     alpha_f = 1. / (1. + integScheme);
//     alpha_m = 0.5 * (3. - integScheme) / (1. + integScheme);
//     gamma = 0.5 + alpha_m - alpha_f;

//     // alpha_f = 1.;
//     // alpha_m = 1.;
//     // gamma = 1.;

//     getline(inputData,line);getline(inputData,line);getline(inputData,line);
//     getline(inputData,line);getline(inputData,line);getline(inputData,line);

//     //Read and print out field forces
//     inputData >> fieldForces[0] >> fieldForces[1] >> fieldForces[2];
    
//     mirrorData << "Field Forces x         = " << fieldForces[0] << std::endl;
//     mirrorData << "Field Forces y         = " << fieldForces[1] << std::endl;
//     mirrorData << "Field Forces z         = " << fieldForces[2] << std::endl \
//                << std::endl;
    
//     getline(inputData,line);getline(inputData,line);getline(inputData,line);
//     getline(inputData,line);getline(inputData,line);

//     fluidParameters.setViscosity(viscInf);
//     fluidParameters.setDensity(rhoInf);
//     fluidParameters.setTimeStep(dTime);
//     fluidParameters.setSpectralRadius(integScheme);
//     fluidParameters.setFieldForce(fieldForces);

//     //Read and print out Drag and lift coeficients
//     inputData >> computeDragAndLift >> numberOfLines; 

//     dragAndLiftBoundary.reserve(numberOfLines);
//     for (int i = 0; i < numberOfLines; ++i)
//     {
//         int aux;
//         inputData >> aux; 
//         dragAndLiftBoundary.push_back(aux);
//     }

//     mirrorData << "Compute Drag and Lift  = " << computeDragAndLift<< std::endl;
//     mirrorData << "Number of Lines  = " << numberOfLines << std::endl;
//     for (int i = 0; i < numberOfLines; ++i)
//     {
//         mirrorData << "Lines  = " << dragAndLiftBoundary[i] << std::endl;
//     }

//     getline(inputData,line);getline(inputData,line);getline(inputData,line);
//     getline(inputData,line);getline(inputData,line);

//     //Read and print out Printing results
//     inputData >> printVelocity;              getline(inputData,line);
//     inputData >> printPressure;              getline(inputData,line);
//     inputData >> printVorticity;             getline(inputData,line);
//     inputData >> printMeshVelocity;          getline(inputData,line);
//     inputData >> printMeshDisplacement;      getline(inputData,line);
//     inputData >> printJacobian;              getline(inputData,line);
//     inputData >> printProcess;              

//     mirrorData << "PrintVelocity              = " << printVelocity << std::endl;
//     mirrorData << "PrintPressure              = " << printPressure << std::endl;
//     mirrorData << "PrintVorticity             = " << printVorticity 
//                << std::endl;
//     mirrorData << "PrintMeshVelocity          = " << printMeshVelocity
//                << std::endl;
//     mirrorData << "PrintMeshDisplacement      = " << printMeshDisplacement
//                << std::endl;
//     mirrorData << "PrintJacobian              = " << printJacobian << std::endl;
//     mirrorData << "PrintProcess               = " << printProcess << std::endl 
//                << std::endl;

   
//     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     //++++++++++++++++++++++++READING ISOGEOMETRIC MESH+++++++++++++++++++++++
//     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//     getline(file2, line); getline(file2, line); getline(file2, line);

//     // number of patches 
//     file2 >> numPatches;

//     IsoPar_.reserve(numPatches);

//     getline(file2, line); getline(file2, line); getline(file2, line);
//     getline(file2, line); getline(file2, line);

//     // number of control points 
//     file2 >> numCP;    
    
//     mirrorData << "Number of Control Points " << numCP << std::endl;
//     mirrorData << "Number of patches " << numPatches << std::endl;

//     nodes_.reserve(numCP);

//     getline(file2, line); getline(file2, line); getline(file2, line);
//     getline(file2, line); getline(file2, line);
   
//     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     //+++++++++++++++++++CONTROL POINTS AND KNOT VECTOR+++++++++++++++++++++++
//     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   
//     int index = 0;
//     int numCPlocal; //number of control points in the patch
//     double wei;

//     for (int ipatch = 0; ipatch < numPatches; ipatch ++){

//         file2 >> numCPlocal;

//         getline(file2, line); getline(file2, line); getline(file2, line);
//         getline(file2, line); getline(file2, line);
            
//         for (int i = 0; i < numCPlocal; i++){
//             double x[3];
//             file2 >> x[0] >> x[1] >> x[2] >> wei;
//             std::getline(file2, line);
//             Node *node = new Node(x,index++,wei);
//             nodes_.push_back(node);
//         };
      
//         mirrorData << "Data in the " << ipatch << " PATCH" << std::endl;
//         mirrorData << "Number of control points " << numCPlocal << std::endl;
        
//         for (int i = (index - numCPlocal) ; i<index; i++){
//             double *x;
//             x = nodes_[i]->getCoordinates();  

//             double x_[3];
//             for (int i=0; i<3; i++) x_[i] = x[i];

//             for (int j=0; j<3; j++){
//                 mirrorData << x_[j] << " ";
//             };
//             double w = nodes_[i]-> getWeightPC(); 
//             mirrorData << w << std::endl;
            
//             nodes_[i] -> setVelocity(velocityInf);
//             nodes_[i] -> setPreviousVelocity(velocityInf);
//             double u[3];
//             u[0] = 0.; u[1] = 0.; u[2] = 0.;
//             nodes_[i] -> setPreviousMeshVelocityComponent(0,0.);
//             nodes_[i] -> setPreviousMeshVelocityComponent(1,0.);
//             nodes_[i] -> setPreviousMeshVelocityComponent(2,0.);
//             nodes_[i] -> setMeshVelocity(u);
                
//             nodes_[i] -> setPreviousCoordinates(0,x_[0]);
//             nodes_[i] -> setPreviousCoordinates(1,x_[1]);
//             nodes_[i] -> setPreviousCoordinates(2,x_[2]);
//         };

//         getline(file2, line); getline(file2, line); getline(file2, line);
//         getline(file2, line); 

//         int ncp[3],deg[3];

//         file2 >> ncp[0] >> ncp[1] >> ncp[2];

//         getline(file2, line); getline(file2, line); getline(file2, line); 
//         getline(file2, line); getline(file2, line);

//         file2 >> deg[0] >> deg[1] >> deg[2];

//         IsoParameters *isopar = new IsoParameters(ipatch,deg,ncp);
//         IsoPar_.push_back(isopar);

//         int dim_u = ncp[0] + deg[0] + 1;
//         int dim_v = ncp[1] + deg[1] + 1;
//         int dim_t = ncp[2] + deg[2] + 1;

//         double *uknot_;
//         uknot_ = new double[dim_u];

//         double *vknot_;
//         vknot_ = new double[dim_v];

//         double *tknot_;
//         tknot_ = new double[dim_t];

//         getline(file2, line); getline(file2, line); getline(file2, line); 
//         getline(file2, line); getline(file2, line);

//         for (int i = 0; i < dim_u; i++) {
//             file2 >> uknot_[i];
//             getline(file2,line);
//         }

//         getline(file2, line); getline(file2, line); getline(file2, line); 
//         getline(file2, line); 

//         for (int j = 0; j < dim_v; j++) {
//             file2 >> vknot_[j];
//             getline(file2,line);
//         }

//         getline(file2, line); getline(file2, line); getline(file2, line); 
//         getline(file2, line); 

//         for (int k = 0; k < dim_t; k++) {
//             file2 >> tknot_[k];
//             getline(file2,line);
//         }

//         IsoPar_[ipatch] -> setuKnot(uknot_);
//         IsoPar_[ipatch] -> setvKnot(vknot_);
//         IsoPar_[ipatch] -> settKnot(tknot_);

//         mirrorData << "Number of Control Points x: " << ncp[0] << " Function Degree: " << deg [0] << std::endl;
//         mirrorData << "Number of Control Points y: " << ncp[1] << " Function Degree: " << deg [1] << std::endl;
//         mirrorData << "Number of Control Points w: " << ncp[2] << " Function Degree: " << deg [2] << std::endl;
        
//         mirrorData << "uKnot " << std::endl;
//         double uknot[dim_u];
//         for (int i = 0; i < dim_u; i++) uknot[i] = uknot_[i];
//         for (int i = 0; i< dim_u; i++) {
//             mirrorData << uknot[i] << std::endl;
//         }

//         mirrorData << "vKnot " << std::endl;
//         double  vknot[dim_v];
//         for (int i = 0; i < dim_v; i++) vknot[i] = vknot_[i];
//         for (int i = 0; i< dim_v; i++) {
//             mirrorData << vknot[i] << std::endl;
//         }

//         mirrorData << "tKnot " << std::endl;
//         double tknot[dim_t];
//         for (int i = 0; i < dim_t; i++) tknot[i] = tknot_[i];
//         for (int i = 0; i< dim_t; i++) {
//             mirrorData << tknot[i] << std::endl;
//         }

//         getline(file2, line); getline(file2, line); 
//         getline(file2, line); getline(file2, line); 


//     }; //iPatch

    

//     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     //++++++++++++++++++++++++++++++++ELEMENTS++++++++++++++++++++++++++++++++
//     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//     int ngp = 0; // the index of control points
//     numIsoElem = 0; // counts the number of elements  
    
//     for (int ipatch = 0; ipatch < numPatches; ipatch++){

//         int ncp[3],deg[3]; 
//         deg[0] = IsoPar_[ipatch] -> getDegree(0);
//         deg[1] = IsoPar_[ipatch] -> getDegree(1);
//         deg[2] = IsoPar_[ipatch] -> getDegree(2);
//         ncp[0] = IsoPar_[ipatch] -> getNcp(0);
//         ncp[1] = IsoPar_[ipatch] -> getNcp(1);
//         ncp[2] = IsoPar_[ipatch] -> getNcp(2);
//         double *uknot;
//         double *vknot;
//         double *tknot;
//         uknot = IsoPar_[ipatch] -> getuKnot();
//         vknot = IsoPar_[ipatch] -> getvKnot();
//         tknot = IsoPar_[ipatch] -> gettKnot();

//         for (int k = 0; k < ncp[2]; k++){ 
//             for (int j = 0; j < ncp[1]; j++){   
//                 for (int i = 0; i < ncp[0]; i++ ) {
//                     int INC[3];
//                     INC[0] = i;
//                     INC[1] = j;
//                     INC[2] = k;
//                     nodes_[ngp++] -> setINC(INC);                      
//                     // finding the elements 
//                     if ((i >= deg[0]) && (j >= deg[1]) && (k >= deg[2])) {
//                         if ( ((uknot[i])!= (uknot[i+1])) && ((vknot[j])!= (vknot[j+1]))  && ((tknot[k])!= (tknot[k+1]))) {
//                             numIsoElem ++;
//                         }
//                     } 
//                 }
//             }
//         }
//     }

//     //number of elements in each physical direction
//     int numElemPatchx;
//     int numElemPatchy;  
//     int numElemPatchz;

//     //direction x
//     for (int ipatch = 0; ipatch < numPatches; ipatch++){
//         numElemPatchx = 0;
//         int ncp,deg; 
//         deg = IsoPar_[ipatch] -> getDegree(0);
//         ncp = IsoPar_[ipatch] -> getNcp(0);
//         double *uknot;
//         uknot = IsoPar_[ipatch] -> getuKnot();

//         for (int i = 0; i < ncp; i++ ) {
//             //finding elements in each direction
//             if (i >= deg) {
//                 if (uknot[i]!= uknot[i+1]) {
//                     numElemPatchx++;
//                 }
//             }
//         }
//         IsoPar_[ipatch] -> setNumElemPatchx(numElemPatchx);
//     }

//     //direction y
//     for (int ipatch = 0; ipatch < numPatches; ipatch++){
//         numElemPatchy = 0;
//         int ncp,deg; 
//         deg = IsoPar_[ipatch] -> getDegree(1);
//         ncp = IsoPar_[ipatch] -> getNcp(1);
//         double *vknot;
//         vknot = IsoPar_[ipatch] -> getvKnot();

//         for (int j = 0; j < ncp; j++ ) {
//             //finding elements in each direction
//             if (j >= deg) {
//                 if (vknot[j]!= vknot[j+1]) {
//                     numElemPatchy++;
//                 }
//             }
//         }
//         IsoPar_[ipatch] -> setNumElemPatchy(numElemPatchy);
//     }

//     //direction z
//     for (int ipatch = 0; ipatch < numPatches; ipatch++){
//         numElemPatchz = 0;
//         int ncp,deg; 
//         deg = IsoPar_[ipatch] -> getDegree(2);
//         ncp = IsoPar_[ipatch] -> getNcp(2);
//         double *tknot;
//         tknot = IsoPar_[ipatch] -> gettKnot();

//         for (int j = 0; j < ncp; j++ ) {
//             //finding elements in each direction
//             if (j >= deg) {
//                 if (tknot[j]!= tknot[j+1]) {
//                     numElemPatchz++;
//                 }
//             }
//         }
//         IsoPar_[ipatch] -> setNumElemPatchz(numElemPatchz);
//     }


//     std::vector<Elements *>   elementsAux_;
//     elementsAux_.reserve(numIsoElem);
    
//     int ElemType; // 0 - FEM elementos, 1 - IGA elementos   
//     ngp = 0;      // the index of control pointst
//     index = 0;   // the index of elements

//     for (int ipatch = 0; ipatch < numPatches; ipatch++){

//         int ncp[3],deg[3]; 
//         deg[0] = IsoPar_[ipatch] -> getDegree(0);
//         deg[1] = IsoPar_[ipatch] -> getDegree(1);
//         deg[2] = IsoPar_[ipatch] -> getDegree(2);
//         ncp[0] = IsoPar_[ipatch] -> getNcp(0);
//         ncp[1] = IsoPar_[ipatch] -> getNcp(1);
//         ncp[2] = IsoPar_[ipatch] -> getNcp(2);
//         double *uknot;
//         double *vknot;
//         double *tknot;
//         uknot = IsoPar_[ipatch] -> getuKnot();
//         vknot = IsoPar_[ipatch] -> getvKnot();
//         tknot = IsoPar_[ipatch] -> gettKnot();
//         int dim_u = ncp[0] + deg[0] + 1;
//         int dim_v = ncp[1] + deg[1] + 1;
//         int dim_t = ncp[2] + deg[2] + 1;

//         for (int k = 0; k < ncp[2]; k++){
//             for (int j = 0; j < ncp[1]; j++){               
//                 for (int i = 0; i < ncp[0]; i++ ) {                 
//                     ngp++;                      
//                     // finding the elements 
//                     if ((i >= deg[0]) && (j >= deg[1]) && (k >= deg[2])) {
//                         if (((uknot[i])!= (uknot[i+1])) && ((vknot[j])!= (vknot[j+1])) && ((tknot[k])!= (tknot[k+1]))) {
//                             int *connect;
//                             connect = new int[27];
//                             for (int kloc = 0; kloc <= deg[2] ; kloc++){                                
//                                 for (int jloc = 0; jloc <= deg[1] ; jloc++){
//                                     for (int iloc = 0; iloc <=  deg[0]; iloc++){
//                                         // global function number
//                                         int ngf = ngp - jloc*ncp[0] - kloc*ncp[1]*ncp[0] - iloc - 1; 
//                                         // local number function
//                                         int nlf = 26 - jloc*(deg[0] + 1) - kloc*(deg[1] + 1)*(deg[0] + 1) - iloc;
//                                         connect[nlf] = ngf;
//                                     }
//                                 } 
//                             }
//                             Elements *el = new Elements(index++,connect,nodes_,1,fluidParameters,IsoPar_,ipatch); 
//                             elementsAux_.push_back(el);
//                             for (int k = 0; k <27; k++){
//                                 nodes_[connect[k]] -> pushInverseIncidence(index);
//                             }; 
//                         }
//                     }      
//                 }
//             }
//         }
//     }

//     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//     //+++++++++++++++++++++++++++++BONDARIES++++++++++++++++++++++++++++++++++
//     //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   
//     //read number of boundary groups and number of elements that belongs to boundaries 
//     file2 >> numBoundariesIso >> numBoundElemIso;

//     //numBoundElemIso = number of control points by boundary - 
//                         //degree of the fuction in the parametric direction
//                         //But can be less than this if it has repeated knots                        

//     mirrorData << "Number of Boundary groups (IA MESH) = " << numBoundariesIso << std::endl;
//     mirrorData << "Number of Boundary Elements (IA MESH)= " << numBoundElemIso << std::endl;
    

//     boundary_.reserve(numBoundElemIso);
  
//     //read boundary information  
//     index = 0;
//     for (int ibound=0; ibound< numBoundariesIso; ibound++){

//         getline(file2,line);getline(file2,line);getline(file2,line);
//         getline(file2,line);getline(file2,line);getline(file2,line);

//         int constrain[3];
//         double value[3];
//         int numGroup;
                
//         file2 >> numGroup >> constrain[0] >> constrain[1] >> constrain[2] \
//                   >> value[0] >> value[1] >> value[2];

//         getline(file2,line);getline(file2,line);getline(file2,line);
//         getline(file2,line);getline(file2,line);

//         //Icpbc contain:
//         //Indexes of the control points with boundary conditions
//         //plan with boundary conditions (xy = 0, xz = 1, yz = 2)
//         //belonging patch
//         int Icpbc[numGroup][dimension+2]; 
        
//         for (int i = 0; i < numGroup; i++){
//           file2 >> Icpbc[i][0] >> Icpbc[i][1] >> Icpbc[i][2] >> Icpbc[i][3] >> Icpbc[i][4];   
//         };
       
//         int ncp[3],deg[3]; 

//         //any of the control points belong to the same patch
//         deg[0] = IsoPar_[Icpbc[0][4]] -> getDegree(0);
//         deg[1] = IsoPar_[Icpbc[0][4]] -> getDegree(1);
//         deg[2] = IsoPar_[Icpbc[0][4]] -> getDegree(2);
//         ncp[0] = IsoPar_[Icpbc[0][4]] -> getNcp(0);
//         ncp[1] = IsoPar_[Icpbc[0][4]] -> getNcp(1);
//         ncp[2] = IsoPar_[Icpbc[0][4]] -> getNcp(2);
//         double *uknot;
//         double *vknot;
//         double *tknot;
//         uknot = IsoPar_[Icpbc[0][4]] -> getuKnot();
//         vknot = IsoPar_[Icpbc[0][4]] -> getvKnot();
//         tknot = IsoPar_[Icpbc[0][4]] -> gettKnot();

//         int contcp = 0;
//         int ncp0,ncp1,ncp2;
//         for (int aux = 0; aux < Icpbc[0][4]; aux++){
//             ncp0 = IsoPar_[aux] -> getNcp(0);
//             ncp1 = IsoPar_[aux] -> getNcp(1);
//             ncp2 = IsoPar_[aux] -> getNcp(2);
//             contcp += ncp0 * ncp1 * ncp2;
//         }

//        //Finding the boundary connectivity
//         for (int k = 0; k < numGroup; k++){
            
//             int *connectB;
//             connectB = new int[9];

//             int i = Icpbc[k][0];
//             int j = Icpbc[k][1];
//             int w = Icpbc[k][2];
//             // face xy
//             if (Icpbc[k][3] == 0){
//                 if ((i >= (Icpbc[0][0] + deg[0])) &&  (j >= (Icpbc[0][1] + deg[1])) ) {
//                     if ((uknot[i] != uknot[i+1]) && (vknot[j] != vknot[j+1]) ) { 
//                         int contcon = 8;
//                         for (int l = 0; l <= deg[1]; l++){
//                             for (int z = 0; z <= deg[0]; z++){  
//                                 connectB[contcon--] = contcp + ncp[0]*ncp[1]*w + ncp[0]*j + i - z -l*ncp[0];
//                             }
//                         }
//                         Boundaries *bound = new Boundaries(connectB, index++,constrain, value, ibound);
//                         boundary_.push_back(bound);
//                     }
//                 }
//             } 
//             // face xz
//             if (Icpbc[k][3] == 1){
//                 if ((i >= (Icpbc[0][0] + deg[0])) &&  (w >= (Icpbc[0][2] + deg[2])) ) {
//                     if ((uknot[i] != uknot[i+1]) && (tknot[w] != tknot[w+1]) ) { 
//                         int contcon = 8;
//                         for (int l = 0; l <= deg[2]; l++){
//                             for (int z = 0; z <= deg[0]; z++){  
//                                 connectB[contcon--] = contcp + ncp[0]*ncp[1]*w + ncp[0]*j + i - z - l*ncp[0]*ncp[1];
//                             }
//                         }
//                         Boundaries *bound = new Boundaries(connectB, index++,constrain, value, ibound);
//                         boundary_.push_back(bound);
//                     }
//                 }
//             }

//             // face yz
//             if (Icpbc[k][3] == 2){
//                 if ((j >= (Icpbc[0][1] + deg[1])) &&  (w >= (Icpbc[0][2] + deg[2])) ) {
//                     if ((vknot[j] != vknot[j+1]) && (tknot[w] != tknot[w+1]) ) { 
//                         int contcon = 8;
//                         for (int l = 0; l <= deg[1]; l++){
//                             for (int z = 0; z <= deg[2]; z++){  
//                                 connectB[contcon--] = contcp + ncp[0]*ncp[1]*w + ncp[0]*j + i - l*ncp[0] - z*ncp[0]*ncp[1];
//                             }
//                         }

//                         Boundaries *bound = new Boundaries(connectB, index++,constrain, value, ibound);
//                         boundary_.push_back(bound);

//                     }
//                 }
//             }                    
//         }//for
//     };  

//     numBoundElemIso = index;

    
//     domainDecompositionMETISIso(elementsAux_);
    
//     if (rank == 0){
//         for (int i = 0; i < numIsoElem; ++i) delete elementsAux_[i];
//         elementsAux_.clear();
//     };

//     // Each rank reads and prints the information about its elements
//     MPI_Barrier(PETSC_COMM_WORLD);

//     std::string result;
//     std::ostringstream convert;

//     convert << rank+000;
//     result = convert.str();
//     std::string s = "mesh"+result+".dat";

//     std::ifstream mesh(s.c_str(), std::ios_base::out);

//     mesh >> numElem;

//     elements_.reserve(numElem);
    
//     for (int i = 0; i < numElem; i++){
        
//         int npatch;
//         int ind_;
//         int ElemType;
//         int *connect;
//         connect = new int[27];

//         mesh >>  ind_ >> ElemType >> npatch >> connect[0] >> connect[1] >> connect[2] >> connect[3] 
//                  >> connect[4]  >> connect[5] >> connect[6] >> connect[7] >> connect[8] >> connect[9]
//                  >> connect[10] >> connect[11] >> connect[12] >> connect[13] >> connect[14]
//                  >> connect[15] >> connect[16] >> connect[17] >> connect[18] >> connect[19]
//                  >> connect[20] >> connect[21] >> connect[22] >> connect[23] >> connect[24]
//                  >> connect[25] >> connect[26];

//         Elements *el = new Elements(ind_,connect,nodes_,ElemType,fluidParameters,IsoPar_,npatch);
//         elements_.push_back(el);

//     };

//     MPI_Barrier(PETSC_COMM_WORLD);
//     MPI_Allreduce(&numElem,&numIsoElem,1,MPI_INT,MPI_SUM,PETSC_COMM_WORLD);
   
//     //Sets boundary constrains
//     for (int ibound = 0; ibound < numBoundElemIso; ibound++){
      
//         int *connectB;
//         connectB = boundary_[ibound] -> getBoundaryConnectivity();

//         int no1 = connectB[0];
//         int no2 = connectB[1];
//         int no3 = connectB[2];
//         int no4 = connectB[3];
//         int no5 = connectB[4];
//         int no6 = connectB[5];
//         int no7 = connectB[6];
//         int no8 = connectB[7];
//         int no9 = connectB[8];
      
//         if ((boundary_[ibound] -> getConstrain(0) == 1) || (boundary_[ibound] -> getConstrain(0) == 3)){

//             nodes_[no1] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
//                                      boundary_[ibound] -> getConstrainValue(0));
//             nodes_[no2] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
//                                      boundary_[ibound] -> getConstrainValue(0));
//             nodes_[no3] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
//                                      boundary_[ibound] -> getConstrainValue(0));
//             nodes_[no4] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
//                                      boundary_[ibound] -> getConstrainValue(0));
//             nodes_[no5] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
//                                      boundary_[ibound] -> getConstrainValue(0));
//             nodes_[no6] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
//                                      boundary_[ibound] -> getConstrainValue(0));
//             nodes_[no7] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
//                                      boundary_[ibound] -> getConstrainValue(0));
//             nodes_[no8] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
//                                      boundary_[ibound] -> getConstrainValue(0));
//             nodes_[no9] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
//                                      boundary_[ibound] -> getConstrainValue(0));
//         };

//         if((boundary_[ibound] -> getConstrain(1) == 1) || (boundary_[ibound] -> getConstrain(1) == 3)){
            
//             nodes_[no1] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
//                                      boundary_[ibound] -> getConstrainValue(1));
//             nodes_[no2] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
//                                      boundary_[ibound] -> getConstrainValue(1));
//             nodes_[no3] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
//                                      boundary_[ibound] -> getConstrainValue(1));
//             nodes_[no4] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
//                                      boundary_[ibound] -> getConstrainValue(1));
//             nodes_[no5] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
//                                      boundary_[ibound] -> getConstrainValue(1));
//             nodes_[no6] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
//                                      boundary_[ibound] -> getConstrainValue(1));
//             nodes_[no7] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
//                                      boundary_[ibound] -> getConstrainValue(1));
//             nodes_[no8] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
//                                      boundary_[ibound] -> getConstrainValue(1));
//             nodes_[no9] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
//                                      boundary_[ibound] -> getConstrainValue(1));
//         };   


//         if((boundary_[ibound] -> getConstrain(2) == 1) || (boundary_[ibound] -> getConstrain(2) == 3)){
            
//             nodes_[no1] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
//                                      boundary_[ibound] -> getConstrainValue(2));
//             nodes_[no2] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
//                                      boundary_[ibound] -> getConstrainValue(2));
//             nodes_[no3] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
//                                      boundary_[ibound] -> getConstrainValue(2));
//             nodes_[no4] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
//                                      boundary_[ibound] -> getConstrainValue(2));
//             nodes_[no5] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
//                                      boundary_[ibound] -> getConstrainValue(2));
//             nodes_[no6] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
//                                      boundary_[ibound] -> getConstrainValue(2));
//             nodes_[no7] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
//                                      boundary_[ibound] -> getConstrainValue(2));
//             nodes_[no8] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
//                                      boundary_[ibound] -> getConstrainValue(2));
//             nodes_[no9] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
//                                      boundary_[ibound] -> getConstrainValue(2));
//         };    
        
//     };

//     //Print nodal constrains
//     for (int i=0; i< numCP; i++){
//         mirrorData<< "Constrains " << i
//                   << " " << nodes_[i] -> getConstrains(0) \
//                   << " " << nodes_[i] -> getConstrainValue(0) \
//                   << " " << nodes_[i] -> getConstrains(1) \
//                   << " " << nodes_[i] -> getConstrainValue(1) \
//                   << " " << nodes_[i] -> getConstrains(2) \
//                   << " " << nodes_[i] -> getConstrainValue(2) << std::endl;
//     }; 

//     //Sets fluid elements and sides on interface boundaries (just for drag and lift coeficients)
//     //This part must be change in the overlap problem
//     for (int i = 0; i < numBoundElemIso; i++){
//         for (int k = 0; k < numberOfLines; k++){
//             if (boundary_[i] -> getBoundaryGroup() == dragAndLiftBoundary[k]){
//                 int *connectB;
//                 connectB = boundary_[i] -> getBoundaryConnectivity();
//                 for (int j=0; j<numElem; j++){
//                     int *connect;
//                     connect = elements_[j] -> getConnectivity();  
//                     int flag = 0;
//                     int side[9];
//                     for (int k=0; k<27; k++){
//                         if ((connectB[0] == connect[k]) || 
//                             (connectB[1] == connect[k]) ||
//                             (connectB[2] == connect[k]) ||
//                             (connectB[3] == connect[k]) ||
//                             (connectB[4] == connect[k]) ||
//                             (connectB[5] == connect[k]) ||
//                             (connectB[6] == connect[k]) ||
//                             (connectB[7] == connect[k]) ||
//                             (connectB[8] == connect[k])) {
//                             side[flag] = k;
//                             flag++;
//                         };
//                     };
//                     if (flag == 9){
//                         // sides:
//                         // 0 - front
//                         // 1 - right
//                         // 2 - back
//                         // 3 - left
//                         // 4 - bottom
//                         // 5 - upper

//                         //Element numeration (front to back - bottom to upper - left to right)
//                         // firs layer of points 0-1-2-3-4-5-6-7-8
//                         // second layer of points 9-10-11-12-13-14-15-16-17
//                         // last lauer of points 18-19-20-21-22-23-24-25-26
//                         boundary_[i] -> setElement(elements_[j] -> getIndex());
//                         //Sets element index and side
//                         if ((side[0]==4) || (side[1]==4) || (side[2]==4) || (side[3]==4) || (side[4]==4) || (side[5]==4) || (side[6]==4) || (side[7]==4) || (side[8]==4) ) {
//                             boundary_[i] -> setElementSide(0);
//                             elements_[j] -> setElemSideInBoundary(0);
//                         };
//                         if ((side[0]==14) || (side[1]==14) || (side[2]==14) || (side[3]==14) || (side[4]==14) || (side[5]==14) || (side[6]==14) || (side[7]==14) || (side[8]==14)) {
//                             boundary_[i] -> setElementSide(1);
//                             elements_[j] -> setElemSideInBoundary(1);
//                         };
//                         if ((side[0]==22) || (side[1]==22) || (side[2]==22) || (side[3]==22) || (side[4]==22) || (side[5]==22) || (side[6]==22) || (side[7]==22) || (side[8]==22)) {
//                             boundary_[i] -> setElementSide(2);
//                             elements_[j] -> setElemSideInBoundary(2);
//                         };
//                         if ((side[0]==12) || (side[1]==12) || (side[2]==12) || (side[3]==12) || (side[4]==12) || (side[5]==12) || (side[6]==12) || (side[7]==12) || (side[8]==12)) { 
//                             boundary_[i] -> setElementSide(3);
//                             elements_[j] -> setElemSideInBoundary(3);
//                         };
//                         if ((side[0]==10) || (side[1]==10) || (side[2]==10) || (side[3]==10) || (side[4]==10) || (side[5]==10) || (side[6]==10) || (side[7]==10) || (side[8]==10)) {
//                             boundary_[i] -> setElementSide(4);
//                             elements_[j] -> setElemSideInBoundary(4);
//                         };
//                         if ((side[0]==16) || (side[1]==16) || (side[2]==16) || (side[3]==16) || (side[4]==16) || (side[5]==16) || (side[6]==16) || (side[7]==16) || (side[8]==16)) {
//                             boundary_[i] -> setElementSide(5);
//                             elements_[j] -> setElemSideInBoundary(5);
//                         };
//                     }; 
//                 }; //elements
//             }; //drag and lift boundary
//         }; //number of lines          
//     }; //number of boundary elements

 
//     // stores the matching control points of each control point CORRETO
//     double distol = 1.e-6;
//     for (int i = 0; i < numCP; i++){
//         std::vector<int> matchingCP_;
//         for (int j = 0; j < numCP; j++){
//             if (i != j){
//                 double *xi,*xj;
//                 xi = nodes_[i]->getCoordinates(); 
//                 xj = nodes_[j]->getCoordinates();  
//                 double dist2 = ((xi[0]- xj[0]) * (xi[0]- xj[0])) + ((xi[1]- xj[1]) * (xi[1]- xj[1])) + ((xi[2]- xj[2]) * (xi[2]- xj[2]));
//                 if (dist2 <= distol) {
//                     matchingCP_.push_back(j);
//                 }
//             }
//         }   

//         int sizeMcp = matchingCP_.size();
//         int *matchingCP;
//         matchingCP = new int[sizeMcp];
//         for (int i = 0; i<sizeMcp; i++) matchingCP[i] = matchingCP_[i];
//         nodes_[i]->setMatchingCP(matchingCP);
//         nodes_[i]->setSizeMcp(sizeMcp);

//         matchingCP_.clear();   
//     }
    
//     // creates a new numeration of control points (descarting the repeated ones)
//     NCNumberNodes = 0; //non coincidente number nodes;
//     int flag[numCP] = {};

//     for (int inodes = 0; inodes < numCP; inodes++) {
        
//         int *matchingCP_;
//         matchingCP_ = nodes_[inodes]-> getMatchingCP();
        
//         int sizemcp = nodes_[inodes]-> getSizeMcp();
        
//         if (sizemcp == 0) {
//             nodes_[inodes]-> setnewcon(NCNumberNodes++);
//         } else {
//             if (flag[inodes] == 0){
//                 flag[inodes] = 1;
//                 nodes_[inodes]-> setnewcon(NCNumberNodes);
//                 for (int i = 0; i < sizemcp; i++) {
//                     if (flag[matchingCP_[i]] == 0){
//                         nodes_[matchingCP_[i]] -> setnewcon(NCNumberNodes);
//                         flag[matchingCP_[i]] = 1;
//                     } 
//                 }
//                 NCNumberNodes++;
//             }
//         }
//     }

//     BezierConnectivity();

            
//      //Closing the file
//      file2.close();
//      if (deleteFiles)
//         system((remove2 + inputFile).c_str());

    // return;
};