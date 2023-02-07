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

#include "FluidData.h"
#include "Glue.h"

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
    typedef typename FluidMesh::Parameters      Parameters;
    // Class glueParameters
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
    int NCnumNodesGlueZoneFine;             // Number of NC fine nodes in the gluing zone       
    int NCNumberNodesC;						// NUmber of non coincident control points coarse mesh
    int NCNumberNodesF;			            // NUmber of non coincident control points fine mesh
    std::vector<int>elementsGlueZoneFine_;  // Vector with fine elements in the gluing zone
    std::vector<int>nodesGlueZoneFine_;     // Vector with fine nodes in the gluing zone
    std::vector<int>NCnodesGlueZoneFine_;   // Vector with non coincidents fine nodes in the gluing zone
    std::vector<int>elementsGlueZoneCoarse_;// Vector with coarse elements in the gluing zone
    std::vector<int>nodesGlueZoneCoarse_;   // Vector with coarse nodes in the gluing zone
    
    // Integration parameters 
    double integScheme;            // Integration scheme (0 - max. dissipation; 1 - without dissipation)
         
    //MPI parameters
    int rank;                   				// Name of the process
    std::pair<idx_t*,idx_t*> domDecompCoarse;   // Coarse Model Domain Decomposition
    std::pair<idx_t*,idx_t*> domDecompFine;     // Fine Model Domain Decomposition
    double pi = M_PI;
    Mat               A,F;
    Vec               b, u, All;
    PetscErrorCode    ierr;
    PetscInt          Istart, Iend, Ii, Ij, Ione, iterations, *dof;
    KSP               ksp;
    PC                pc;
    VecScatter        ctx;
    PetscScalar       val;
    PetscViewer       viewer;
    
    int NumBezierNodesLagrange;   

public:
    // Sets the coarse and fine mesh models
    void setFluidModels_FEM_ISO(FluidMesh& coarse, FluidMesh& fine);
    
    // Computes and stores the element boxes for coarse mesh for improving the correspondence searching process
    void setElementBoxes_ISO();
    
    // Computes and store the signaled distance function from Nodes to a defined FEM boundary
    void setSignaledDistance_FEM_ISO();
    void searchMinDist_FEM(int &inode, double *pointCoord, int &ibound, double &mindist);

    // Defines the Gluing zone
    void setGluingZone_FEM_ISO();

    // Computes the energy weight function
    void setWeightFunction_FEM_ISO();

    //Sets the nodal and integration points correspondence in the coarse mesh
    void setCorrespondenceFine_FEM_ISO();

   	//Sets the Dirichelet Constrains in the domain
   	void setDirichletConstrain_FEM_ISO(std::vector<int> &dofTemp);

    //Sets the Dirichelet Constrains in the domain
   	void setDirichletConstrainLaplace_FEM_ISO(std::vector<int> &dofTemp);
   
    // Searchs point correspondence in the coarse mesh
    void searchPointCorrespondence_ISO(double *x,std::vector<Nodes *> nodes,
                                      std::vector<Element *> elements, 
                                      std::vector<IsoParameters* > isopar, int numElem,
                                      double *xsiC, int &elemC, int elSearch);
    
    //Solves the Arlequin Problem
    int solveArlequinProblem_FEM_ISO(int iterNumber,double tolerance);

    int solveArlequinProblemLaplace_FEM_ISO(int iterNumber,double tolerance);

    //Assemble system
    void setMatVecValuesCoarse_ISO();
    void setMatVecValuesFine_FEM();
    void setMatVecValuesLagrangeFine_FEM_ISO();
    void setMatVecValuesLagrangeCoarse_FEM_ISO();

    void setMatVecValuesCoarseLaplace_ISO();
    void setMatVecValuesFineLaplace_FEM();
    void setMatVecValuesLagrangeFineLaplace_FEM_ISO();
    void setMatVecValuesLagrangeCoarseLaplace_FEM_ISO();

    // Compute and print drag and lift coefficients
    void dragAndLiftCoefficients_FEM(std::ofstream& dragLift, int &iTimeStep);


    //Prints the results for Paraview post-processing
    void printResults_FEM_ISO(int step);

    void printResultsLaplace_FEM_ISO(int step);

    //Print the results of the gluing in the integration points for Paraview post-processing
    void printResultsIP_FEM_ISO(int step);

};

#endif
