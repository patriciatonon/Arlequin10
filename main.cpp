//------------------------------------------------------------------------------
//--------------------------Universidade de Sao Paulo---------------------------
//----------------------Escola de Engenharia de Sao Carlos----------------------
//----------------Departamento de Engenharia de Estruturas - SET----------------
//------------------------------Sao Carlos - 2018-------------------------------
//------------------------------------------------------------------------------
 
//------------------------------------------------------------------------------
//----Software developed for analysis of 2D/3D incompressible flow problems-----
//----in the Arbitrary Lagrangian-Eulerian (ALE) description. The governing-----
//--equations are approximated by the Finite Element Method or by Isogeometric--
//-----Analysis with a mixed formulation, stabilized elements and quadratic-----
//----------------approximation for both velocity and pressure.-----------------
//------------------------------------------------------------------------------
  
//------------------------------------------------------------------------------
//---------------------------------Developed by---------------------------------
//-------Jeferson Wilian Dossa Fernandes, Rodolfo Andre Kuche Sanches and-------
//--------------------------------Patricia Tonon--------------------------------
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------STARTS MAIN PROGRAM-----------------------------
//------------------------------------------------------------------------------

static char help[] = "Solves the Incompressible flow problem";

// C++ standard libraries
#include <fstream> 
  
// Developed Header Files
#include "src/FluidData.h"
#include "src/Arlequin.h" 


int main(int argc, char **args) {
 
    // Starts main program invoking PETSc
    PetscInitialize(&argc, &args, (char*)0, help);

    int rank, size; 
 
    MPI_Comm_rank(PETSC_COMM_WORLD, &rank);
    MPI_Comm_size(PETSC_COMM_WORLD, &size);
  
    if (rank == 0){
        std::cout << "2D/3D Incompressible Flows Numerical Analysis" << std::endl;
        std::cout << "Starting.." << std::endl;
    };

    // Defines the problem dimension
    const int dimension = 3;

    //Type definition
    typedef FluidData<dimension>     FluidModel;
    typedef Arlequin<dimension>      Arlequin;
 
    //Create problem variables 
    FluidModel coarseModel,fineModel; 
    Arlequin   ArlequinProblem; 

//===============================================================================
//=================================PROBLEM MESH===================================
//================================================================================

    MPI_Barrier(PETSC_COMM_WORLD);   
    //data reading functions need two files from coarse or fine mesh:  
    //1- Fluid Flow Data
    //2- FEM mesh (function dataReading) or IGA mesh (function dataReadingIso)

    //FEM/ISO
    //3D LAPLACE - FEM/ISO
    // coarseModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/Laplace_iso3d.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/Laplace_fem3d.msh","mirror_fine.txt",0);

    // //2D Laplace - FEM/ISO
    // coarseModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/coarse_iso2.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/fineJ_fem.msh","mirror_fine.txt",0);

    // ArlequinProblem.setFluidModels_FEM_ISO(coarseModel,fineModel);
    // ArlequinProblem.solveArlequinProblemLaplace_FEM_ISO(10, 1.e-6);

    //ISO/ISO

    //2d
    // coarseModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/coarse_iso2.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/fine_iso.msh","mirror_fine.txt",0);


    //3d
    // coarseModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/Laplace_iso3d.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/LaplaceFine_iso3d.msh","mirror_fine.txt",0);

    // ArlequinProblem.setFluidModels_ISO_ISO(coarseModel,fineModel);
    // ArlequinProblem.solveArlequinProblemLaplace_ISO_ISO(10, 1.e-6);


    //FEM/FEM
    //2D
    // coarseModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/coarseJ_fem.msh","mirror_fine.txt",0);
    // fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/fineJ_fem.msh","mirror_fine.txt",0);

    //3D
    coarseModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/coarseR_fem.msh","mirror_coarse.txt",0);
    fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/fineR_fem.msh","mirror_fine.txt",0);

    ArlequinProblem.setFluidModels_FEM_FEM(coarseModel,fineModel);
    ArlequinProblem.solveArlequinProblemLaplace_FEM_FEM(10, 1.e-6);


    //Finalize main program   
    PetscFinalize();
 
    return 0; 
};
 
  




 
