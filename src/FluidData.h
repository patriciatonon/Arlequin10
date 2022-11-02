//------------------------------------------------------------------------------
// 
//        Jeferson W D Fernandes,Rodolfo A K Sanches and Patrícia Tonon
//                             University of Sao Paulo
//                           (C) 2019 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//---------------------------------FLUID DATA-----------------------------------
//------------------------------------------------------------------------------

#ifndef FLUIDDATA_H
#define FLUIDDATA_H

#include "Element.h"
#include "Boundary.h"

#include<cstdlib>
#include<fstream>
#include<iostream>
#include <math.h>       /* atan */

// PETSc libraries
#include <metis.h>
#include <petscksp.h> 

#include <boost/timer.hpp> 
#include <boost/thread.hpp>

// Mounts the incompressible flow problems
template<int DIM>
class FluidData{
public:
	//Defines the classes locally:
    // Class of element
    typedef Element<DIM> Elements;
    // Class of Node 
    typedef typename Elements::Nodes  Node;
    // Class of Boundary
    typedef Boundary<DIM> Boundaries;
    // Class of Fluid Parameters
    typedef FluidParameters<DIM> Parameters;
    // Class of Isogeometric Parameters locally
    typedef IsogeometricParameters<DIM> IsoParameters;

private:
    //MPI variables
    idx_t* part_elem;      				 //Fluid Domain Decomposition - Elements
    idx_t* part_nodes;     				 //Fluid Domain Decomposition - Nodes
    int rank;			   			     //Number of thread 
    int size;			   				 //Parameter to the division of mesh between the threads
    double pi = M_PI;					 //Paralelization parameter 


public:

	// Defines the classes objetcs:
    std::vector<Node *>          nodes_;
    std::vector<IsoParameters *> IsoPar_;
    std::vector<Elements *>      elements_;
    std::vector<Boundaries *>    boundary_;
    Parameters  fluidParameters;    

    //mesh data
    int NumBezierNodes;                   //Number of Bezier Nodes 
    int NCNumberNodes;                    //Number of non coincidents control points in the mesh (for solve the matrixes)
    
    //Print data
    bool printVelocity;
    bool printPressure;
    bool printRealVelocity;
    bool printRealPressure;
    bool printMeshVelocity;
    bool printMeshDisplacement;
    bool printProcess;
    bool printDistFunction;
    bool printGlueZone;
    bool printLagrangeMultipliers;
    bool printEnergyWeightFunction;
    bool printNodalCorrespondence;
    
    //Drag and Lift coefficients
    bool computeDragAndLift;			 //Check if is necessary compute Drag and Lift parameters
    int numberOfLines;                   //Number of the boundaries in the problem to compute Drag and Lift parameters
    std::vector<int> dragAndLiftBoundary;//Store the number of boundary that match with numberOfLines 

    
public:
    //Reads the input file 
    void dataReading_FEM(const std::string& inputFile, const std::string& inputMesh, const std::string& mirror, const bool& deleteFiles);
    void dataReading_ISO(const std::string& inputFile, const std::string& inputMeshIso, const std::string& mirror, const bool& deleteFiles);
    
    //Performs the domain decomposition for parallel processing
    void domainDecompositionMETIS_FEM(std::vector<Elements *> &elem_, int &numFemElem, int &numFemNodes); 
    void domainDecompositionMETIS_ISO(std::vector<Elements *> &elem_, int &numIsoElem, int &numCP); 

    //Export the domain decomposition 
    std::pair<idx_t*,idx_t*> getDomainDecomposition(){
    return std::make_pair(part_elem,part_nodes);};

    //Compute the Bézier connectivity of the IGA elements
    void BezierConnectivity(int &numPatches);

};

#endif
