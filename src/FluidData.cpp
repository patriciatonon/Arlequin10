#include "FluidData.h"

//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//---------------------SUBDIVIDES THE FINITE ELEMENT DOMAIN---------------------
//------------------------------------------------------------------------------
template<int DIM>
void FluidData<DIM>::domainDecompositionMETIS_FEM(std::vector<Elements *> &elem_, int &numFemElem, int &numFemNodes) {
    
    int LNN = 4*DIM - 2;
    
    std::string mirror2;
    mirror2 = "domain_decompositionFEM.txt";
    std::ofstream mirrorData(mirror2.c_str());
    int size;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    idx_t objval;
    idx_t numEl = numFemElem;
    idx_t numNd = numFemNodes;
    idx_t dd = DIM;
    idx_t ssize = size;
    idx_t three = 3;
    idx_t one = 1;
    idx_t elem_start[numEl+1], elem_connec[LNN*numEl];
   

    MPI_Bcast(&numEl,1,MPI_INT,0,PETSC_COMM_WORLD);
    MPI_Bcast(&numNd,1,MPI_INT,0,PETSC_COMM_WORLD);

    part_elem = new idx_t[numEl];
    part_nodes = new idx_t[numNd];

    if (rank == 0){

        for (idx_t i = 0; i < numEl+1; i++){
            elem_start[i]=LNN*i;
        };
        for (idx_t jel = 0; jel < numEl; jel++){
            int *connec;
            connec = elem_[jel] -> getConnectivity();      
            for (idx_t i = 0; i < LNN; i++){
                elem_connec[LNN*jel+i] = connec[i];
            };
        };
    

        //Performs the domain decomposition
        if (size == 1){

            for(int i = 0; i < numFemElem; i++)
                part_elem[i] = 0;


            for(int i = 0; i < numFemNodes; i++)
                part_nodes[i] = 0;
        
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

template<int DIM>
void FluidData<DIM>::domainDecompositionMETIS_ISO(std::vector<Elements *> &elem_, int &numIsoElem, int &numCP) {
    
    int LNN = 18*DIM - 27;

    std::string mirror2;
    mirror2 = "domain_decompositionISO.txt";
    std::ofstream mirrorData(mirror2.c_str());
    
    int size;

    MPI_Comm_size(PETSC_COMM_WORLD, &size);

    idx_t objval;
    idx_t numEl = numIsoElem;
    idx_t numNd = numCP;
    idx_t dd = DIM;
    idx_t ssize = size;
    idx_t three = 3;
    idx_t one = 1;
    idx_t elem_start[numEl+1], elem_connec[LNN*numEl];
    

    MPI_Bcast(&numEl,1,MPI_INT,0,PETSC_COMM_WORLD);
    MPI_Bcast(&numNd,1,MPI_INT,0,PETSC_COMM_WORLD);

    part_elem = new idx_t[numEl];
    part_nodes = new idx_t[numNd];

    if (rank == 0){

        for (idx_t i = 0; i < numEl+1; i++){
            elem_start[i]=LNN*i;
        };
        for (idx_t jel = 0; jel < numEl; jel++){
            int *connec;
            connec = elem_[jel] -> getConnectivity();      
            for (idx_t i = 0; i < LNN; i++){
                elem_connec[LNN*jel+i] = connec[i];
            };
        };
    

        //Performs the domain decomposition

        if (size == 1){

            for(int i = 0; i < numIsoElem; i++)
                 part_elem[i] = 0;
 
            for(int i = 0; i < numCP; i++)
                part_nodes[i] = 0;

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
    };

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
                };
                if (j != 0) {
                    if (rest == 0) {
                        int kp = (j/2);
                        if (i == 0) {
                            Beconnec[kp*numElemx+i][0] = numNodes; 
                            Beconnec[(kp-1)*numElemx+i][6] = numNodes++;  
                        };
                        Beconnec[kp*numElemx+i][1] = numNodes;
                        Beconnec[(kp-1)*numElemx+i][7] = numNodes++;
                        Beconnec[kp*numElemx+i][2] = numNodes;
                        Beconnec[(kp-1)*numElemx+i][8] = numNodes;  
                        if (i != numElemx-1) {
                            Beconnec[kp*numElemx+i+1][0] = numNodes; 
                            Beconnec[(kp-1)*numElemx+i+1][6] = numNodes;       
                        };     
                        numNodes ++;  
                    } else {
                        int kp = ((j-1)/2);
                        if (i == 0) Beconnec[kp*numElemx+i][3] = numNodes++;
                        Beconnec[kp*numElemx+i][4] = numNodes++;
                        Beconnec[kp*numElemx+i][5] = numNodes;
                        if (i != numElemx -1) Beconnec[kp*numElemx+i+1][3] = numNodes;
                        numNodes ++;
                    };
                }; //j!=0
            }; //elemx
        };//elemy
        for (int i = 0; i < numElemx; i++) {
            if (i == 0) Beconnec[(numElemy-1)*numElemx+i][6] = numNodes++;  
            Beconnec[(numElemy-1)*numElemx+i][7] = numNodes++;
            Beconnec[(numElemy-1)*numElemx+i][8] = numNodes;
            if (i != numElemx-1) Beconnec[(numElemy-1)*numElemx+i+1][6] = numNodes;
            numNodes ++;
        };
        int contel = 0;
        for (int aux = 0; aux < ipatch; aux++){
             int numElemxx = IsoPar_[aux] -> getNumElemPatchx();
             int numElemyy = IsoPar_[aux] -> getNumElemPatchy();
             contel +=  numElemxx*numElemyy;
        };  

        for (int iElem = 0; iElem < numElemx*numElemy; iElem++){
            int *Becon_;
            Becon_ = new int[9];
            for (int j=0; j<9; j++){
                Becon_[j] = Beconnec[iElem][j];
            };
            elements_[iElem+contel] -> setBezierConnectivity(Becon_);
        };
    };//patch
    NumBezierNodes = numNodes;
};

template<>
void FluidData<3>::BezierConnectivity(int &numPatches) {

    int numNodes = 0;
    int kne = 0;
   
    for (int ipatch = 0; ipatch < numPatches; ipatch++){
        
        int numElemx = IsoPar_[ipatch] -> getNumElemPatchx();
        int numElemy = IsoPar_[ipatch] -> getNumElemPatchy();
        int numElemz = IsoPar_[ipatch] -> getNumElemPatchz();

        int  Beconnec[numElemx*numElemy*numElemz][27];
        
        for (int k = 0; k <= 2*numElemz; k++){
        	if (k == 0){
		        for (int j = 0; j< 2*numElemy; j++){    
		            for (int i = 0; i< numElemx; i++){
		                int rest = j%2;                
		                if (j == 0) {
		                    if (i == 0) Beconnec[j*numElemx+i][0] = numNodes++;  
		                    Beconnec[j*numElemx+i][1] = numNodes++;
		                    Beconnec[j*numElemx+i][2] = numNodes;
		                    if (i != numElemx-1) Beconnec[j*numElemx+i+1][0] = numNodes;
		                    numNodes ++;
		                };
		                if (j != 0) {
		                    if (rest == 0) {
		                        int kp = (j/2);
		                        if (i == 0) {
		                            Beconnec[kp*numElemx+i][0] = numNodes; 
		                            Beconnec[(kp-1)*numElemx+i][6] = numNodes++;  
		                        };
		                        Beconnec[kp*numElemx+i][1] = numNodes;
		                        Beconnec[(kp-1)*numElemx+i][7] = numNodes++;
		                        Beconnec[kp*numElemx+i][2] = numNodes;
		                        Beconnec[(kp-1)*numElemx+i][8] = numNodes;   
		                        if (i != numElemx-1) {
		                            Beconnec[kp*numElemx+i+1][0] = numNodes; 
		                            Beconnec[(kp-1)*numElemx+i+1][6] = numNodes;       
		                        };       
		                        numNodes ++;
		                    } else {
		                        int kp = ((j-1)/2);
		                        if (i == 0) Beconnec[kp*numElemx+i][3] = numNodes++;
		                        Beconnec[kp*numElemx+i][4] = numNodes++;
		                        Beconnec[kp*numElemx+i][5] = numNodes;
		                        if (i != numElemx -1) Beconnec[kp*numElemx+i+1][3] = numNodes;
		                        numNodes ++;
		                    };
		                }; //j!=0
		            }; //elemx
		        };//elemy
		        for (int i = 0; i < numElemx; i++) {
		            if (i == 0) Beconnec[(numElemy-1)*numElemx+i][6] = numNodes++;  
		            Beconnec[(numElemy-1)*numElemx+i][7] = numNodes++;
		            Beconnec[(numElemy-1)*numElemx+i][8] = numNodes;
		            if (i != numElemx-1) Beconnec[(numElemy-1)*numElemx+i+1][6] = numNodes;
		            numNodes ++;
		        };
        	};

        	if ((k != 0) && (k != 2*numElemz)){
        		int restk = k%2;
        		if (restk != 0){
        			int ind = (k-1)/2;
		    		for (int j = 0; j< 2*numElemy; j++){    
			            for (int i = 0; i< numElemx; i++){
			                int rest = j%2;               
			                if (j == 0) {
			                    if (i == 0) Beconnec[ind*numElemx*numElemy + j*numElemx + i][9] = numNodes++;  
			                    Beconnec[ind*numElemx*numElemy + j*numElemx+i][10] = numNodes++;
			                    Beconnec[ind*numElemx*numElemy + j*numElemx+i][11] = numNodes;
			                    if (i != numElemx-1) Beconnec[ind*numElemx*numElemy + j*numElemx+i+1][9] = numNodes;
			                    numNodes ++;
			                };
			                if (j != 0) {
			                    if (rest == 0) {
			                        int kp = (j/2);
			                        if (i == 0) {
			                            Beconnec[ind*numElemx*numElemy + kp*numElemx+i][9] = numNodes; 
			                            Beconnec[ind*numElemx*numElemy + (kp-1)*numElemx+i][15] = numNodes++;  
			                        };
			                        Beconnec[ind*numElemx*numElemy + kp*numElemx+i][10] = numNodes;
			                        Beconnec[ind*numElemx*numElemy + (kp-1)*numElemx+i][16] = numNodes++;
			                        Beconnec[ind*numElemx*numElemy + kp*numElemx+i][11] = numNodes;
			                        Beconnec[ind*numElemx*numElemy + (kp-1)*numElemx+i][17] = numNodes;			                            
			                        if (i != numElemx-1) {
			                            Beconnec[ind*numElemx*numElemy + kp*numElemx+i+1][9] = numNodes; 
			                            Beconnec[ind*numElemx*numElemy + (kp-1)*numElemx+i+1][15] = numNodes;       
			                        };       
			                        numNodes ++;			                        
			                    } else {
			                        int kp = ((j-1)/2);
			                        if (i == 0) Beconnec[ind*numElemx*numElemy + kp*numElemx+i][12] = numNodes++;
			                        Beconnec[ind*numElemx*numElemy + kp*numElemx+i][13] = numNodes++;
			                        Beconnec[ind*numElemx*numElemy + kp*numElemx+i][14] = numNodes;
			                        if (i != numElemx -1) Beconnec[ind*numElemx*numElemy + kp*numElemx+i+1][12] = numNodes;
			                        numNodes ++;
			                    };			                
			                }; //j!=0
			            }; //elemx
			        };//elemy			      
			        for (int i = 0; i < numElemx; i++) {
			            if (i == 0) Beconnec[ind*numElemx*numElemy + (numElemy-1)*numElemx+i][15] = numNodes++;  
			            Beconnec[ind*numElemx*numElemy + (numElemy-1)*numElemx+i][16] = numNodes++;
			            Beconnec[ind*numElemx*numElemy + (numElemy-1)*numElemx+i][17] = numNodes;
			            if (i != numElemx-1) Beconnec[ind*numElemx*numElemy + (numElemy-1)*numElemx+i+1][15] = numNodes;
			            numNodes ++;
			        };
        		}; // if restz !=0
        		
        		if (restk == 0){
        			int ind = (k/2);
		    		for (int j = 0; j< 2*numElemy; j++){    
			            for (int i = 0; i< numElemx; i++){
			                int rest = j%2;			                                
			                if (j == 0) {
			                    if (i == 0) {
			                    	Beconnec[ind*numElemx*numElemy + j*numElemx+i][0] = numNodes;  
			                    	Beconnec[(ind-1)*numElemx*numElemy + j*numElemx+i][18] = numNodes++; 
			                    };
			                    Beconnec[ind*numElemx*numElemy + j*numElemx+i][1] = numNodes;
			                    Beconnec[(ind-1)*numElemx*numElemy + j*numElemx+i][19] = numNodes++;
			                    Beconnec[ind*numElemx*numElemy + j*numElemx+i][2] = numNodes;
			                    Beconnec[(ind-1)*numElemx*numElemy + j*numElemx+i][20] = numNodes;
			                    if (i != numElemx-1) {
			                    	Beconnec[ind*numElemx*numElemy + j*numElemx+i+1][0] = numNodes;
			                    	Beconnec[(ind-1)*numElemx*numElemy + j*numElemx+i+1][18] = numNodes;
			                    };
			                    numNodes ++;
			                };
			                if (j != 0){
			                    if (rest == 0) {
			                        int kp = (j/2);
			                        if (i == 0) {
			                            Beconnec[ind*numElemx*numElemy + kp*numElemx+i][0] = numNodes; 
			                            Beconnec[(ind-1)*numElemx*numElemy + kp*numElemx+i][18] = numNodes;
			                            Beconnec[ind*numElemx*numElemy +(kp-1)*numElemx+i][6] = numNodes; 
			                            Beconnec[(ind-1)*numElemx*numElemy + (kp-1)*numElemx+i][24] = numNodes++;  
			                        };
			                        Beconnec[ind*numElemx*numElemy + kp*numElemx+i][1] = numNodes;
			                        Beconnec[(ind-1)*numElemx*numElemy + kp*numElemx+i][19] = numNodes;
			                        Beconnec[ind*numElemx*numElemy +(kp-1)*numElemx+i][7] = numNodes;
			                        Beconnec[(ind-1)*numElemx*numElemy +(kp-1)*numElemx+i][25] = numNodes++;
			                        Beconnec[ind*numElemx*numElemy +kp*numElemx+i][2] = numNodes;
			                        Beconnec[(ind-1)*numElemx*numElemy + kp*numElemx+i][20] = numNodes;
			                        Beconnec[ind*numElemx*numElemy +(kp-1)*numElemx+i][8] = numNodes;
			                        Beconnec[(ind-1)*numElemx*numElemy +(kp-1)*numElemx+i][26] = numNodes;			                            
			                        if (i != numElemx-1) {
			                            Beconnec[ind*numElemx*numElemy +kp*numElemx+i+1][0] = numNodes; 
			                            Beconnec[(ind-1)*numElemx*numElemy + kp*numElemx+i+1][18] = numNodes; 
			                            Beconnec[ind*numElemx*numElemy +(kp-1)*numElemx+i+1][6] = numNodes; 
			                            Beconnec[(ind-1)*numElemx*numElemy + (kp-1)*numElemx+i+1][24] = numNodes;       
			                        };      
			                        numNodes ++;			                        
			                    } else {
			                        int kp = ((j-1)/2);
			                        if (i == 0) {
										Beconnec[ind*numElemx*numElemy +kp*numElemx+i][3] = numNodes;
										Beconnec[(ind-1)*numElemx*numElemy +kp*numElemx+i][21] = numNodes++;
			                        };
			                        Beconnec[ind*numElemx*numElemy + kp*numElemx+i][4] = numNodes;
			                        Beconnec[(ind-1)*numElemx*numElemy + kp*numElemx+i][22] = numNodes++;
			                        Beconnec[ind*numElemx*numElemy + kp*numElemx+i][5] = numNodes;
			                        Beconnec[(ind-1)*numElemx*numElemy + kp*numElemx+i][23]= numNodes;
			                        if (i != numElemx -1) {
			                        	Beconnec[ind*numElemx*numElemy + kp*numElemx+i+1][3] = numNodes;
			                        	Beconnec[(ind-1)*numElemx*numElemy + kp*numElemx+i+1][21] = numNodes;
			                        };
			                        numNodes ++;
			                    }; //else
			                }; //j!=0
			            }; //elemx
			        };//elemy
			        for (int i = 0; i < numElemx; i++) {
			            if (i == 0) {
			            	Beconnec[ind*numElemx*numElemy + (numElemy-1)*numElemx+i][6] = numNodes;  
			            	Beconnec[(ind-1)*numElemx*numElemy +(numElemy-1)*numElemx+i][24] = numNodes++;  
			            };
			            Beconnec[ind*numElemx*numElemy + (numElemy-1)*numElemx+i][7] = numNodes;
			            Beconnec[(ind-1)*numElemx*numElemy + (numElemy-1)*numElemx+i][25] = numNodes++;
			            Beconnec[ind*numElemx*numElemy + (numElemy-1)*numElemx+i][8] = numNodes;
			            Beconnec[(ind-1)*numElemx*numElemy + (numElemy-1)*numElemx+i][26] = numNodes;
			            if (i != numElemx-1) {
			            	Beconnec[ind*numElemx*numElemy + (numElemy-1)*numElemx+i+1][6] = numNodes;
			            	Beconnec[(ind-1)*numElemx*numElemy + (numElemy-1)*numElemx+i+1][24] = numNodes;
			            };
			            numNodes ++;
			        };    
			    }; // if restz == 0	
	       	}; // k != 0 or k != 2*numElemz
        	
        	if (k == 2*numElemz){
        		int ind = numElemz - 1;
        		for (int j = 0; j< 2*numElemy; j++){    
		            for (int i = 0; i< numElemx; i++){
		                int rest = j%2;               
		                if (j == 0) {
		                    if (i == 0) Beconnec[ind*numElemx*numElemy +j*numElemx+i][18] = numNodes++;  
		                    Beconnec[ind*numElemx*numElemy +j*numElemx+i][19] = numNodes++;
		                    Beconnec[ind*numElemx*numElemy +j*numElemx+i][20] = numNodes;
		                    if (i != numElemx-1) Beconnec[ind*numElemx*numElemy +j*numElemx+i+1][18] = numNodes;
		                    numNodes ++;
		                };
		                if (j != 0) {
		                    if (rest == 0) {
		                        int kp = (j/2);
		                        if (i == 0) {
		                            Beconnec[ind*numElemx*numElemy +kp*numElemx+i][18] = numNodes; 
		                            Beconnec[ind*numElemx*numElemy +(kp-1)*numElemx+i][24] = numNodes++;  
		                        };
		                        Beconnec[ind*numElemx*numElemy +kp*numElemx+i][19] = numNodes;
		                        Beconnec[ind*numElemx*numElemy +(kp-1)*numElemx+i][25] = numNodes++;
		                        Beconnec[ind*numElemx*numElemy +kp*numElemx+i][20] = numNodes;
		                        Beconnec[ind*numElemx*numElemy +(kp-1)*numElemx+i][26] = numNodes;
		                        if (i != numElemx-1) {
		                            Beconnec[ind*numElemx*numElemy +kp*numElemx+i+1][18] = numNodes; 
		                            Beconnec[ind*numElemx*numElemy +(kp-1)*numElemx+i+1][24] = numNodes;       
		                        };      
		                        numNodes ++;    
		                    } else {
		                        int kp = ((j-1)/2);
		                        if (i == 0) Beconnec[ind*numElemx*numElemy +kp*numElemx+i][21] = numNodes++;
		                        Beconnec[ind*numElemx*numElemy +kp*numElemx+i][22] = numNodes++;
		                        Beconnec[ind*numElemx*numElemy +kp*numElemx+i][23] = numNodes;
		                        if (i != numElemx -1) Beconnec[ind*numElemx*numElemy +kp*numElemx+i+1][21] = numNodes;
		                        numNodes ++;
		                    };
		                };//j!=0
		            };//elemx
		        }//elemy
		        for (int i = 0; i < numElemx; i++) {
		            if (i == 0) Beconnec[ind*numElemx*numElemy +(numElemy-1)*numElemx+i][24] = numNodes++;  
		            Beconnec[ind*numElemx*numElemy +(numElemy-1)*numElemx+i][25] = numNodes++;
		            Beconnec[ind*numElemx*numElemy +(numElemy-1)*numElemx+i][26] = numNodes;
		            if (i != numElemx-1) Beconnec[ind*numElemx*numElemy +(numElemy-1)*numElemx+i+1][24] = numNodes;
		            numNodes ++;
		        };
        	};//if k = 2*numElemz

        }; // loop in z direction
		
		int contel = 0;
        for (int aux = 0; aux < ipatch; aux++){
             int numElemxx = IsoPar_[ipatch] -> getNumElemPatchx();
             int numElemyy = IsoPar_[ipatch] -> getNumElemPatchy();
             int numElemzz = IsoPar_[ipatch] -> getNumElemPatchz();
             contel +=  numElemxx*numElemyy*numElemzz;
        };
      

        for (int iElem = 0; iElem < numElemx*numElemy*numElemz; ++iElem){
            
            int *Becon_;
            Becon_ = new int[27];

            for (int j=0; j<27; j++) Becon_[j] = Beconnec[iElem][j];

            elements_[kne] -> setBezierConnectivity(Becon_);
            
            kne ++;
        };
    
    };//patch

    NumBezierNodes = numNodes;
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
    double dTime;                        //Time Step
    //Arlequin Problem
    double arlequinK1;                   //L2
    double arlequinK2;                   //H1
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
     	

    // //Sets boundary constrains
    // for (int ibound = 0; ibound < numFemBoundElem; ibound++){
    //     int *connectB = boundary_[ibound] -> getBoundaryConnectivity();
    //     int no1 = connectB[0];
    //     int no2 = connectB[1];
    //     int no3 = connectB[2];

    //     for (int i = 0; i < dim; i++){
           
    //         if ((boundary_[ibound] -> getConstrain(i) == 1) || 
    //             (boundary_[ibound] -> getConstrain(i) == 3)){
    //             nodes_[no1] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
    //                                          boundary_[ibound] -> getConstrainValue(i));
    //             nodes_[no2] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
    //                                          boundary_[ibound] -> getConstrainValue(i));
    //             nodes_[no3] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
    //                                          boundary_[ibound] -> getConstrainValue(i));
    //         };
    //     };
          
    // };


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


        if (boundary_[ibound] -> getConstrain(0) == 4){

            double constrainvalue1 = sin(2*pi*nodes_[no1]->getCoordinateValue(0));
            double constrainvalue2 = sin(2*pi*nodes_[no2]->getCoordinateValue(0));
            double constrainvalue3 = sin(2*pi*nodes_[no3]->getCoordinateValue(0));

            double zero = 0.;
            
            nodes_[no1] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                        zero);
            nodes_[no1] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                        constrainvalue1);
            nodes_[no2] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                        zero);
            nodes_[no2] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                         constrainvalue2);
            nodes_[no3] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                         zero);
            nodes_[no3] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                         constrainvalue3);



        };

        if (boundary_[ibound] -> getConstrain(0) == 5){

            double constrainvalue1 = -cos(2*pi*nodes_[no1]->getCoordinateValue(1));
            double constrainvalue2 = -cos(2*pi*nodes_[no2]->getCoordinateValue(1));
            double constrainvalue3 = -cos(2*pi*nodes_[no3]->getCoordinateValue(1));
            double zero = 0.;
            
            nodes_[no1] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                         zero);
            nodes_[no1] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                         constrainvalue1);
            nodes_[no2] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                         zero);
            nodes_[no2] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                         constrainvalue2);
            nodes_[no3] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                         zero);
            nodes_[no3] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                         constrainvalue3);

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
void FluidData<3>::dataReading_FEM(const std::string& inputFile,const std::string& inputMeshFem, 
                                   const std::string& mirror,const bool& deleteFiles){
    
    int DIM = 3;
    //variables
    //Fluid data
    double pressInf;       			     //Undisturbed pressure 
    double rhoInf;         			     //Density
    double viscInf;        				 //Viscosity
    double fieldForces[DIM]; 			 //Field forces (constant) 
    double velocityInf[DIM];             //Infinity velocity
    double dTime;                        //Time Step
    //Arlequin Problem
    double arlequinK1;                   //L2
    double arlequinK2;                   //H1
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
        
        double wei,x[DIM];
        
        file >> x[0] >> x[1] >> x[2] >> wei;
        std::getline(file, line);
        
        Node *node = new Node(x, index++, wei); 
        nodes_.push_back(node);
    };

    mirrorData << "Nodal Coordinates " << std::endl;
    
    for (int i = 0 ; i<numFemNodes; i++){       
        
        double *x = nodes_[i]->getCoordinates();       
        // double x_[DIM];
        // for (int j = 0; j < dim; j++) x_[j] = x[j];
        for (int j = 0; j < DIM; j++) mirrorData << x[j] << " ";
        mirrorData << std::endl;
        
        //Setting the initial nodal values variables
        double u[DIM] = {};
        nodes_[i] -> setVelocity(velocityInf);
        nodes_[i] -> setPreviousVelocity(velocityInf);
        nodes_[i] -> setPreviousMeshVelocity(u);
        nodes_[i] -> setMeshVelocity(u);
        nodes_[i] -> setPreviousCoordinates(x);

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
        connect = new int[10];
        
        file >> connect[0] >> connect[1] >> connect[2] >> connect[3] >> connect[4] >> connect[5] >>
                connect[6] >> connect[7] >> connect[8] >> connect[9];
            
        std::getline(file, line);
             
        connect[0] -= 1; connect[1] -= 1; connect[2] -= 1;
        connect[3] -= 1; connect[4] -= 1; connect[5] -= 1;
        connect[6] -= 1; connect[7] -= 1; connect[8] -= 1;
        connect[9] -= 1;

        Elements *el = new Elements(index++,connect,nodes_,0,fluidParameters,IsoPar_,1000); 
        
        elements_.push_back(el);
        //0 = Element finit type of element 
        //1000 = path number for FEM mesh
        for (int k = 0; k < 10; k++){
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
        int constrain[DIM];
        double value[DIM];
        int numGroup;

        file >> numGroup >> constrain[0] >> constrain[1] >> constrain[2] 
                         >> value[0] >> value[1] >> value[2];

        getline(file,line);getline(file,line);getline(file,line);
        getline(file,line);getline(file,line);
       
        for (int i = 0; i < numGroup; i++){        	
        	
            int *connectB;
        	connectB = new int[6];
            
            file >> connectB[0] >> connectB[1] >> connectB[2] >> connectB[3] >> connectB[4] >> connectB[5];
            
            getline(file,line);
            
            connectB[0] -= 1; connectB[1] -= 1; connectB[2] -= 1;
            connectB[3] -= 1; connectB[4] -= 1; connectB[5] -= 1;

            Boundaries *bound = new Boundaries(connectB, index++,  
                                               constrain, value, ibound);
            boundary_.push_back(bound);
        };
    };     
     	

    // //Sets boundary constrains
    // for (int ibound = 0; ibound < numFemBoundElem; ibound++){
        
    //     int *connectB = boundary_[ibound] -> getBoundaryConnectivity();
    //     int no1 = connectB[0];
    //     int no2 = connectB[1];
    //     int no3 = connectB[2];
    //     int no4 = connectB[3];
    //     int no5 = connectB[4];
    //     int no6 = connectB[5];

    //     for (int i = 0; i < DIM; i++){
           
    //         if ((boundary_[ibound] -> getConstrain(i) == 1) || 
    //             (boundary_[ibound] -> getConstrain(i) == 3)){
    //             nodes_[no1] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
    //                                          boundary_[ibound] -> getConstrainValue(i));
    //             nodes_[no2] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
    //                                          boundary_[ibound] -> getConstrainValue(i));
    //             nodes_[no3] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
    //                                          boundary_[ibound] -> getConstrainValue(i));
    //             nodes_[no4] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
    //                                          boundary_[ibound] -> getConstrainValue(i));
    //             nodes_[no5] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
    //                                          boundary_[ibound] -> getConstrainValue(i));
    //             nodes_[no6] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
    //                                          boundary_[ibound] -> getConstrainValue(i));
    //         };
    //     };
          
    // };


    for (int ibound = 0; ibound < numFemBoundElem; ibound++){
        
        int *connectB = boundary_[ibound] -> getBoundaryConnectivity();
        int no1 = connectB[0];
        int no2 = connectB[1];
        int no3 = connectB[2];
        int no4 = connectB[3];
        int no5 = connectB[4];
        int no6 = connectB[5];


        for (int i = 0; i < DIM; i++){
           
            if ((boundary_[ibound] -> getConstrain(i) == 1) || 
                (boundary_[ibound] -> getConstrain(i) == 3)){
                nodes_[no1] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no2] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no3] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no4] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no5] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no6] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
            };
        };
        
        
        // for (int i = 0; i < dim; i++){
           
            if (boundary_[ibound] -> getConstrain(0) == 4){

                double constrainvalue1 = sin(2*pi*nodes_[no1]->getCoordinateValue(0));
                double constrainvalue2 = sin(2*pi*nodes_[no2]->getCoordinateValue(0));
                double constrainvalue3 = sin(2*pi*nodes_[no3]->getCoordinateValue(0));
                double constrainvalue4 = sin(2*pi*nodes_[no4]->getCoordinateValue(0));
                double constrainvalue5 = sin(2*pi*nodes_[no5]->getCoordinateValue(0));
                double constrainvalue6 = sin(2*pi*nodes_[no6]->getCoordinateValue(0));
                double zero = 0.;
                
                nodes_[no1] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                             zero);
                nodes_[no1] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                             constrainvalue1);
                nodes_[no1] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
                                             zero);
                
                nodes_[no2] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                             zero);
                nodes_[no2] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                             constrainvalue2);
                nodes_[no2] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
                                             zero);
                
                nodes_[no3] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                             zero);
                nodes_[no3] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                             constrainvalue3);
                nodes_[no3] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
                                             zero);
                
                nodes_[no4] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                             zero);
                nodes_[no4] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                             constrainvalue4);
                nodes_[no4] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
                                             zero);

                nodes_[no5] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                             zero);
                nodes_[no5] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                             constrainvalue5);
                nodes_[no5] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
                                             zero);

                nodes_[no6] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                             zero);
                nodes_[no6] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                             constrainvalue6);
                nodes_[no6] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
                                             zero);



            };

            if (boundary_[ibound] -> getConstrain(0) == 5){

                double constrainvalue1 = -cos(2*pi*nodes_[no1]->getCoordinateValue(1));
                double constrainvalue2 = -cos(2*pi*nodes_[no2]->getCoordinateValue(1));
                double constrainvalue3 = -cos(2*pi*nodes_[no3]->getCoordinateValue(1));
                double constrainvalue4 = -cos(2*pi*nodes_[no4]->getCoordinateValue(1));
                double constrainvalue5 = -cos(2*pi*nodes_[no5]->getCoordinateValue(1));
                double constrainvalue6 = -cos(2*pi*nodes_[no6]->getCoordinateValue(1));
                double zero = 0.;
                
                nodes_[no1] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                             zero);
                nodes_[no1] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                             constrainvalue1);
                nodes_[no1] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
                                             zero);
                
                nodes_[no2] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                             zero);
                nodes_[no2] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                             constrainvalue2);
                nodes_[no2] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
                                             zero);
                
                nodes_[no3] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                             zero);
                nodes_[no3] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                             constrainvalue3);
                nodes_[no3] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
                                             zero);
                
                nodes_[no4] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                             zero);
                nodes_[no4] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                             constrainvalue4);
                nodes_[no4] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
                                             zero);

                nodes_[no5] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                             zero);
                nodes_[no5] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                             constrainvalue5);
                nodes_[no5] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
                                             zero);

                nodes_[no6] -> setConstrains(0,boundary_[ibound] -> getConstrain(0),
                                             zero);
                nodes_[no6] -> setConstrains(1,boundary_[ibound] -> getConstrain(1),
                                             constrainvalue6);
                nodes_[no6] -> setConstrains(2,boundary_[ibound] -> getConstrain(2),
                                             zero);
               
            };


        // };
          
    };

    //Print nodal constrains
    for (int i = 0; i < numFemNodes; i++){
        mirrorData<< "Constrains " << i
                  << " " << nodes_[i] -> getConstrains(0) 
                  << " " << nodes_[i] -> getConstrainValue(0) 
                  << " " << nodes_[i] -> getConstrains(1) 
                  << " " << nodes_[i] -> getConstrainValue(1) 
                  << " " << nodes_[i] -> getConstrains(2) 
                  << " " << nodes_[i] -> getConstrainValue(2) << std::endl;
    }; 

    
    //Sets fluid elements and sides on interface boundaries
    int boundElem[numFemBoundElem]; 
    for (int i = 0; i< numFemBoundElem; i++){
        int *connectB= boundary_[i] -> getBoundaryConnectivity();
        for (int j = 0; j < numFemElem; j++){               
            int *connect = elements_[j] -> getConnectivity();              
            int flag = 0;
            int side[6];
            for (int k = 0; k < 10; k++){
                if ((connectB[0] == connect[k]) || 
                    (connectB[1] == connect[k]) ||
                    (connectB[2] == connect[k]) ||
                    (connectB[3] == connect[k]) ||
                    (connectB[4] == connect[k]) ||
                    (connectB[5] == connect[k])) {
                    side[flag] = k;
                    flag++;
                };
            };
            if (flag == 6){
                boundary_[i] -> setElement(elements_[j] -> getIndex());
                //Sets element index and side
                if (((side[0]==1) || (side[1]==1) || (side[2]==1) || (side[3]==1) || (side[4]==1) || (side[5]==1)) &
                    ((side[0]==2) || (side[1]==2) || (side[2]==2) || (side[3]==2) || (side[4]==2) || (side[5]==2)) &
                    ((side[0]==3) || (side[1]==3) || (side[2]==3) || (side[3]==3) || (side[4]==3) || (side[5]==3))) {
                    boundary_[i] -> setElementSide(0);
                    boundElem[i] = j;
                };
                if (((side[0]==0) || (side[1]==0) || (side[2]==0) || (side[3]==0) || (side[4]==0) || (side[5]==0)) &
                    ((side[0]==2) || (side[1]==2) || (side[2]==2) || (side[3]==2) || (side[4]==2) || (side[5]==2)) &
                    ((side[0]==3) || (side[1]==3) || (side[2]==3) || (side[3]==3) || (side[4]==3) || (side[5]==3))) {
                    boundary_[i] -> setElementSide(1);
                    boundElem[i] = j;
                };
                if (((side[0]==0) || (side[1]==0) || (side[2]==0) || (side[3]==0) || (side[4]==0) || (side[5]==0)) &
                    ((side[0]==1) || (side[1]==1) || (side[2]==1) || (side[3]==1) || (side[4]==1) || (side[5]==1)) &
                    ((side[0]==3) || (side[1]==3) || (side[2]==3) || (side[3]==3) || (side[4]==3) || (side[5]==3))) {
                    boundary_[i] -> setElementSide(2);
                    boundElem[i] = j;
                };
                if (((side[0]==0) || (side[1]==0) || (side[2]==0) || (side[3]==0) || (side[4]==0) || (side[5]==0)) &
                    ((side[0]==1) || (side[1]==1) || (side[2]==1) || (side[3]==1) || (side[4]==1) || (side[5]==1)) &
                    ((side[0]==2) || (side[1]==2) || (side[2]==2) || (side[3]==2) || (side[4]==2) || (side[5]==2))) {
                    boundary_[i] -> setElementSide(3);
                    boundElem[i] = j;
                };  
            };
        };
    }; 

    std::vector<int> elemSide;
    std::vector<int> elemBound;

    for (int i = 0; i < numFemElem; i++){
        for (int j = 0; j < numFemBoundElem; j++){
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
    double dTime;                        //Time Step
    //Arlequin Problem
    double arlequinK1;                   //L2
    double arlequinK2;                   //H1
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
    fluidParameters.setVelocityInf(velocityInf);

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

template<>
void FluidData<3>::dataReading_ISO(const std::string& inputFile,const std::string& inputMeshIso,
                                  const std::string& mirror,const bool& deleteFiles){
    
    int DIM = 3;                                
    
    //Fluid data
    double pressInf;       			     //Undisturbed pressure 
    double rhoInf;         			     //Density
    double viscInf;        				 //Viscosity
    double fieldForces[DIM]; 			 //Field forces (constant) 
    double velocityInf[DIM];             //Infinity velocity
    double dTime;                        //Time Step
    //Arlequin Problem
    double arlequinK1;                   //L2
    double arlequinK2;                   //H1
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
    fluidParameters.setVelocityInf(velocityInf);

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
    double wei;
    for (int ipatch = 0; ipatch < numPatches; ipatch ++){

        file2 >> numCPlocal;

        getline(file2, line); getline(file2, line); getline(file2, line);
        getline(file2, line); getline(file2, line);
            
        for (int i = 0; i < numCPlocal; i++){
            double x[DIM];
            file2 >> x[0] >> x[1] >> x[2] >> wei;
            std::getline(file2, line);
            Node *node = new Node(x,index++,wei);
            nodes_.push_back(node);
        };
      
        mirrorData << "Data in the " << ipatch << " PATCH" << std::endl;
        mirrorData << "Number of control points " << numCPlocal << std::endl;
        
        for (int i = (index - numCPlocal) ; i < index; i++){
           
            double *x = nodes_[i]->getCoordinates();     
            for (int j = 0; j < DIM; j++) mirrorData << x[j] << " ";
            mirrorData << nodes_[i]-> getWeightPC() << std::endl;
            
            //Setting the initial nodal values variables
            double u[DIM] = {};
            nodes_[i] -> setVelocity(velocityInf);
            nodes_[i] -> setPreviousVelocity(velocityInf);
            nodes_[i] -> setPreviousMeshVelocity(u);
            nodes_[i] -> setMeshVelocity(u);
            nodes_[i] -> setPreviousCoordinates(x);
            
        };

        getline(file2, line); getline(file2, line); getline(file2, line);
        getline(file2, line); 

        int ncp[DIM],deg[DIM];

        file2 >> ncp[0] >> ncp[1] >> ncp[2];

        getline(file2, line); getline(file2, line); getline(file2, line); 
        getline(file2, line); getline(file2, line);

        file2 >> deg[0] >> deg[1] >> deg[2];

        IsoParameters *isopar = new IsoParameters(ipatch,deg,ncp);
        IsoPar_.push_back(isopar);

        int dim_u = ncp[0] + deg[0] + 1;
        int dim_v = ncp[1] + deg[1] + 1;
        int dim_t = ncp[2] + deg[2] + 1;

        double *uknot_;
        uknot_ = new double[dim_u];
        double *vknot_;
        vknot_ = new double[dim_v];
        double *tknot_;
        tknot_ = new double[dim_t];

        getline(file2, line); getline(file2, line);
        getline(file2, line); getline(file2, line);
         getline(file2, line); 

        for (int i = 0; i < dim_u; i++) {
            file2 >> uknot_[i];
            getline(file2,line);
        };

        getline(file2, line); getline(file2, line); 
        getline(file2, line); getline(file2, line); 

        for (int j = 0; j < dim_v; j++) {
            file2 >> vknot_[j];
            getline(file2,line);
        };

        getline(file2, line); getline(file2, line); 
        getline(file2, line); getline(file2, line); 

        for (int k = 0; k < dim_t; k++) {
            file2 >> tknot_[k];
            getline(file2,line);
        };

        mirrorData << " Number of Control Points x: " << ncp[0] 
                   << " Function Degree: " << deg [0] << std::endl;
        mirrorData << " Number of Control Points y: " << ncp[1] 
                   << " Function Degree: " << deg [1] << std::endl;
        mirrorData << " Number of Control Points z: " << ncp[2] 
                   << " Function Degree: " << deg [2] << std::endl;
        
        mirrorData << "uKnot " << std::endl;
        double uknotP[dim_u];
        for (int i=0;i <dim_u; i++) uknotP[i] = uknot_[i];
        for (int i = 0; i< dim_u; i++) mirrorData << uknotP[i] << std::endl;

        double vknotP[dim_v];
        for (int i=0;i <dim_v; i++) vknotP[i] = vknot_[i];
        mirrorData << "vKnot " << std::endl;
        for (int i = 0; i< dim_v; i++) mirrorData << vknotP[i] << std::endl;

        double tknotP[dim_t];
        for (int i=0;i <dim_t; i++) tknotP[i] = tknot_[i];
        mirrorData << "tKnot " << std::endl;
        for (int i = 0; i< dim_t; i++) mirrorData << tknotP[i] << std::endl;

        IsoPar_[ipatch] -> setuKnot(uknot_);
        IsoPar_[ipatch] -> setvKnot(vknot_);
        IsoPar_[ipatch] -> settKnot(tknot_);

        getline(file2, line); getline(file2, line); 
        getline(file2, line); getline(file2, line);  


    }; //iPatch


    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++ELEMENTS++++++++++++++++++++++++++++++++
    //++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

    //Read the element data from each patch

    int ngp = 0; // the index of control points
    numIsoElem = 0; // counts the number of elements
    for (int ipatch = 0; ipatch < numPatches; ipatch++){
        
        int ncp[DIM],deg[DIM];
        for (int i = 0; i < DIM; i++){
            deg[i] = IsoPar_[ipatch] -> getDegree(i);
            ncp[i] = IsoPar_[ipatch] -> getNcp(i);
        }; 
                
        double *uknot= IsoPar_[ipatch] -> getuKnot();
        double *vknot= IsoPar_[ipatch] -> getvKnot();  
        double *tknot= IsoPar_[ipatch] -> gettKnot();     
        
        for (int k = 0; k < ncp[2]; k++){ 
            for (int j = 0; j < ncp[1]; j++){   
                for (int i = 0; i < ncp[0]; i++ ) {
                    int INC[DIM];
                    INC[0] = i;
                    INC[1] = j;
                    INC[2] = k;
                    nodes_[ngp++] -> setINC(INC);                      
                    // finding the elements 
                    if ((i >= deg[0]) && 
                        (j >= deg[1]) && 
                        (k >= deg[2])) {

                        if ( ((uknot[i])!= (uknot[i+1])) && 
                             ((vknot[j])!= (vknot[j+1])) && 
                             ((tknot[k])!= (tknot[k+1]))) {
                            numIsoElem ++;
                        };
                    };
                };
            };
        };
    };

    elements_.reserve(numIsoElem);

    //number of elements in each patch direction
    int numElemPatchx;
    int numElemPatchy;  
    int numElemPatchz;

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

    //direction z
    for (int ipatch = 0; ipatch < numPatches; ipatch++){
        numElemPatchz = 0;
        int ncp,deg; 
        deg = IsoPar_[ipatch] -> getDegree(2);
        ncp = IsoPar_[ipatch] -> getNcp(2);
        double *tknot;
        tknot = IsoPar_[ipatch] -> gettKnot();

        for (int j = 0; j < ncp; j++ ) {
            //finding elements in each direction
            if (j >= deg){
                if (tknot[j]!= tknot[j+1]) {
                    numElemPatchz++;
                };
            };
        };
        IsoPar_[ipatch] -> setNumElemPatchz(numElemPatchz);
    };


    int ElemType; // 0 - FEM elementos, 1 - IA elementos   
    ngp = 0;     // the index of control points
    index = 0;   // the index of elements

    for (int ipatch = 0; ipatch < numPatches; ipatch++){
        
        int ncp[DIM],deg[DIM];
        for (int i = 0; i < DIM; i++){

            deg[i] = IsoPar_[ipatch] -> getDegree(i);
            ncp[i] = IsoPar_[ipatch] -> getNcp(i);
        }; 
                
        double *uknot= IsoPar_[ipatch] -> getuKnot();
        double *vknot= IsoPar_[ipatch] -> getvKnot();  
        double *tknot= IsoPar_[ipatch] -> gettKnot();    
        
        for (int k = 0; k < ncp[2]; k++){
            for (int j = 0; j < ncp[1]; j++){               
                for (int i = 0; i < ncp[0]; i++ ) {                 
                    ngp++;                      
                    // finding the elements 
                    if ((i >= deg[0]) && 
                        (j >= deg[1]) && 
                        (k >= deg[2])) {
                        if (((uknot[i])!= (uknot[i+1])) && 
                            ((vknot[j])!= (vknot[j+1])) && 
                            ((tknot[k])!= (tknot[k+1]))) {
                            
                            int *connect;
                            connect = new int[27];
                            
                            for (int kloc = 0; kloc <= deg[2] ; kloc++){                                
                                for (int jloc = 0; jloc <= deg[1] ; jloc++){
                                    for (int iloc = 0; iloc <=  deg[0]; iloc++){
                                        // global function number
                                        int ngf = ngp - jloc*ncp[0] - kloc*ncp[1]*ncp[0] - iloc - 1; 
                                        // local number function
                                        int nlf = 26 - jloc*(deg[0] + 1) - kloc*(deg[1] + 1)*(deg[0] + 1) - iloc;
                                        connect[nlf] = ngf;
                                    };
                                };
                            };
                            
                            Elements *el = new Elements(index++,connect,nodes_,1,fluidParameters,IsoPar_,ipatch); 
                            elements_.push_back(el);
                            
                            for (int k = 0; k <27; k++){
                                nodes_[connect[k]] -> pushInverseIncidence(index);
                            }; 
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

        int constrain[DIM];
        double value[DIM];
        int numGroup;
                
        file2 >> numGroup >> constrain[0] >> constrain[1] >> constrain[2] 
                          >> value[0] >> value[1] >> value[2];

        getline(file2,line);getline(file2,line);getline(file2,line);
        getline(file2,line);getline(file2,line);

        //Icpbc contain:
        //Indexes of the control points with boundary conditions
        //plan with boundary conditions (xy = 0, xz = 1, yz = 2)
        //belonging patch
        int Icpbc[numGroup][DIM+2]; 
        for (int i = 0; i < numGroup; i++){
            file2 >> Icpbc[i][0] >> Icpbc[i][1] >> Icpbc[i][2] >> Icpbc[i][3] >> Icpbc[i][4];       
        };
        
        //any of the control points belong to the same patch
        int ncp[DIM],deg[DIM]; 
        for (int i = 0; i < DIM; i++){
            deg[i] = IsoPar_[Icpbc[0][4]] -> getDegree(i);
            ncp[i] = IsoPar_[Icpbc[0][4]] -> getNcp(i);  
        };
        
        double *uknot = IsoPar_[Icpbc[0][4]] -> getuKnot();
        double *vknot = IsoPar_[Icpbc[0][4]] -> getvKnot();
        double *tknot = IsoPar_[Icpbc[0][4]] -> gettKnot();

        int contcp = 0;
        int ncp0,ncp1,ncp2;
        for (int aux = 0; aux < Icpbc[0][4]; aux++){
            ncp0 = IsoPar_[aux] -> getNcp(0);
            ncp1 = IsoPar_[aux] -> getNcp(1);
            ncp2 = IsoPar_[aux] -> getNcp(2);
            contcp += ncp0 * ncp1 * ncp2;
        };


        //Finding the boundary connectivity
        for (int k = 0; k < numGroup; k++){

            int *connectB;
            connectB = new int[9];

            int i = Icpbc[k][0];
            int j = Icpbc[k][1];
            int w = Icpbc[k][2];
            
            // face xy
            if (Icpbc[k][3] == 0){
                if ((i >= (Icpbc[0][0] + deg[0])) &&  
                    (j >= (Icpbc[0][1] + deg[1]))) {
                    if ((uknot[i] != uknot[i+1]) && 
                        (vknot[j] != vknot[j+1]) ) { 
                        int contcon = 8;
                        for (int l = 0; l <= deg[1]; l++){
                            for (int z = 0; z <= deg[0]; z++){  
                                connectB[contcon--] = contcp + ncp[0]*ncp[1]*w + ncp[0]*j + i - z -l*ncp[0];
                            };
                        };
                                                
                        Boundaries *bound = new Boundaries(connectB, index++,constrain, value, ibound);
                        boundary_.push_back(bound);
                    };
                };
            };
            // face xz
            if (Icpbc[k][3] == 1){
                if ((i >= (Icpbc[0][0] + deg[0])) &&  (w >= (Icpbc[0][2] + deg[2])) ) {
                    if ((uknot[i] != uknot[i+1]) && (tknot[w] != tknot[w+1]) ) { 
                        int contcon = 8;
                        for (int l = 0; l <= deg[2]; l++){
                            for (int z = 0; z <= deg[0]; z++){  
                                connectB[contcon--] = contcp + ncp[0]*ncp[1]*w + ncp[0]*j + i - z - l*ncp[0]*ncp[1];
                            };
                        };
                        Boundaries *bound = new Boundaries(connectB, index++,constrain, value, ibound);
                        boundary_.push_back(bound);
                    };
                };
            };
            // face yz
            if (Icpbc[k][3] == 2){
                if ((j >= (Icpbc[0][1] + deg[1])) &&  (w >= (Icpbc[0][2] + deg[2])) ) {
                    if ((vknot[j] != vknot[j+1]) && (tknot[w] != tknot[w+1]) ) { 
                        int contcon = 8;
                        for (int l = 0; l <= deg[1]; l++){
                            for (int z = 0; z <= deg[2]; z++){  
                                connectB[contcon--] = contcp + ncp[0]*ncp[1]*w + ncp[0]*j + i - l*ncp[0] - z*ncp[0]*ncp[1];
                            };
                        };
                        Boundaries *bound = new Boundaries(connectB, index++,constrain, value, ibound);
                        boundary_.push_back(bound);
                    };
                };
            };  
             
        };
    };  

    numBoundElemIso = index;
   

    //Sets nodal boundary constrains
    
    for (int ibound = 0; ibound < numBoundElemIso; ibound++){        
        
        int *connectB = boundary_[ibound] -> getBoundaryConnectivity();

        int no1 = connectB[0];
        int no2 = connectB[1];
        int no3 = connectB[2];
        int no4 = connectB[3];
        int no5 = connectB[4];
        int no6 = connectB[5];
        int no7 = connectB[6];
        int no8 = connectB[7];
        int no9 = connectB[8];

        for (int i = 0; i < DIM; i++){
           
            if ((boundary_[ibound] -> getConstrain(i) == 1) || 
                (boundary_[ibound] -> getConstrain(i) == 3)){
                nodes_[no1] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no2] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no3] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no4] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no5] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no6] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no7] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no8] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
                nodes_[no9] -> setConstrains(i,boundary_[ibound] -> getConstrain(i),
                                             boundary_[ibound] -> getConstrainValue(i));
            };
        }; 
    };

    //Print nodal constrains
    for (int i = 0; i < numCP; i++){
        mirrorData<< "Constrains " << i
                  << " " << nodes_[i] -> getConstrains(0) 
                  << " " << nodes_[i] -> getConstrainValue(0) 
                  << " " << nodes_[i] -> getConstrains(1) 
                  << " " << nodes_[i] -> getConstrainValue(1)
                  << " " << nodes_[i] -> getConstrains(2) 
                  << " " << nodes_[i] -> getConstrainValue(2) << std::endl;
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
                int side[9];
                for (int k = 0; k < 27; k++){
                    if ((connectB[0] == connect[k]) || 
                        (connectB[1] == connect[k]) ||
                        (connectB[2] == connect[k]) ||
                        (connectB[3] == connect[k]) ||
                        (connectB[4] == connect[k]) ||
                        (connectB[5] == connect[k]) ||
                        (connectB[6] == connect[k]) ||
                        (connectB[7] == connect[k]) ||
                        (connectB[8] == connect[k])) {
                        side[flag] = k;
                        flag++;
                    };
                };
                if (flag == 9){
                    boundary_[i] -> setElement(elements_[j] -> getIndex());
                    // sides:
                    // 0 - front, 1 - right, 2 - back, 3 - left, 4 - bottom, 5 - upper
                    // Element numeration 
                    // first layer of points 0-1-2-3-4-5-6-7-8
                    // second layer of points 9-10-11-12-13-14-15-16-17
                    // last layer of points 18-19-20-21-22-23-24-25-26
                    
                    boundary_[i] -> setElement(elements_[j] -> getIndex());
                    //Sets element index and side
                    if ((side[0]==4) || (side[1]==4) || (side[2]==4) || 
                        (side[3]==4) || (side[4]==4) || (side[5]==4) || 
                        (side[6]==4) || (side[7]==4) || (side[8]==4) ) {
                        boundary_[i] -> setElementSide(0);
                        boundElem[i] = j;
                    };
                    if ((side[0]==14) || (side[1]==14) || (side[2]==14) || 
                        (side[3]==14) || (side[4]==14) || (side[5]==14) || 
                        (side[6]==14) || (side[7]==14) || (side[8]==14)) {
                        boundary_[i] -> setElementSide(1);
                        boundElem[i] = j;
                    };
                    if ((side[0]==22) || (side[1]==22) || (side[2]==22) || 
                        (side[3]==22) || (side[4]==22) || (side[5]==22) || 
                        (side[6]==22) || (side[7]==22) || (side[8]==22)) {
                        boundary_[i] -> setElementSide(2);
                        boundElem[i] = j;
                    };
                    if ((side[0]==12) || (side[1]==12) || (side[2]==12) || 
                        (side[3]==12) || (side[4]==12) || (side[5]==12) || 
                        (side[6]==12) || (side[7]==12) || (side[8]==12)) { 
                        boundary_[i] -> setElementSide(3);
                        boundElem[i] = j;
                    };
                    if ((side[0]==10) || (side[1]==10) || (side[2]==10) || 
                        (side[3]==10) || (side[4]==10) || (side[5]==10) || 
                        (side[6]==10) || (side[7]==10) || (side[8]==10)) {
                        boundary_[i] -> setElementSide(4);
                        boundElem[i] = j;
                    };
                    if ((side[0]==16) || (side[1]==16) || (side[2]==16) || 
                        (side[3]==16) || (side[4]==16) || (side[5]==16) || 
                        (side[6]==16) || (side[7]==16) || (side[8]==16)) {
                        boundary_[i] -> setElementSide(5);
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
                double dist2 = 0;
                for (int k = 0; k < DIM; k++) dist2 += ((xi[k]- xj[k]) * (xi[k]- xj[k]));
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


template class FluidData<2>;
template class FluidData<3>;