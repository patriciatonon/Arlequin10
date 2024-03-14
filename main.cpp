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
    const int dimension = 2;

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
    //2D - Airfoil 
    coarseModel.dataReading_ISO("../../mesh/aerofemfem/meshfine_data.txt","../../mesh/aero1/coarse_iso.msh","mirror_coarse.txt",
                                "../../mesh/aero1/COARSECPOutput100032.dat",0);
    fineModel.dataReading_FEM("../../mesh/aerofemfem/meshfine_data.txt","../../mesh/aerofemfem/fine_fem.msh","mirror_fine.txt", "nada.txt", 0);

    
    //teste que nao deu certo cd e cl
    // coarseModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/aero1/coarse1_iso.msh","mirror_coarse.txt",
    //                             "../../mesh/aero1/COARSECPOutput100032.dat",0);
    // fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/aero1/fineRef_fem.msh","mirror_fine.txt",
    //                             "../../mesh/aero1/FineNodeOutput100032.dat",0);
   
    //Cylinder
    //coarseModel.dataReading_ISO("../../mesh/cylinder/meshcoarse_data.txt","../../mesh/cylinder/coarse_iso_cylinder.msh","mirror_coarse.txt", "nada",0);
    // fineModel.dataReading_FEM("../../mesh/cylinder/meshfine_data.txt","../../mesh/cylinder/fine_fem_cylinder_jefe_nre.msh","mirror_fine.txt","nada", 0);

    //2D - cavity
    // coarseModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/coarse_iso2.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/fineJ_fem.msh","mirror_fine.txt",0);

    //3D - cavity
    // coarseModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/Laplace_iso3d.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/Laplace_fem3d1.msh","mirror_fine.txt",0);

    ArlequinProblem.setFluidModels_FEM_ISO(coarseModel,fineModel);
    // ArlequinProblem.solveArlequinProblemLaplace_FEM_ISO(10, 1.e-6);
    ArlequinProblem.solveArlequinProblem_FEM_ISO(3, 1.e-6);

    //ISO/ISO
    //2d
    // coarseModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/coarse_iso2.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/fine_iso.msh","mirror_fine.txt",0);

    //3d
    // coarseModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/Laplace_iso3d.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/LaplaceFine_iso3d.msh","mirror_fine.txt",0);

    // ArlequinProblem.setFluidModels_ISO_ISO(coarseModel,fineModel);
    // ArlequinProblem.solveArlequinProblemLaplace_ISO_ISO(10, 1.e-6);
    // ArlequinProblem.solveArlequinProblem_ISO_ISO(10, 1.e-6);


    //FEM/FEM
    //2D

    // AIRFOIL
    // coarseModel.dataReading_FEM("../../mesh/aerofemfem/meshfine_data.txt","../../mesh/aerofemfem/coarse_fem.msh","mirror_coarse.txt", "nada.txt",0);
    // fineModel.dataReading_FEM("../../mesh/aerofemfem/meshfine_data.txt","../../mesh/aerofemfem/fine_fem.msh","mirror_fine.txt", "nada.txt", 0);

    // coarseModel.dataReading_FEM("../../mesh/meshfinecyl_data.txt","../../mesh/cylinder/coarse_fem_cylinder_jefe.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("../../mesh/meshfinecyl_data.txt","../../mesh/cylinder/fine_fem_cylinder_jefe.msh","mirror_fine.txt",0);

    // cavidade faixa ao redor
    // coarseModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/coarseT_fem.msh","mirror_fine.txt",0);
    // fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/fineT_fem.msh","mirror_fine.txt",0);

    //cavidade faixa embaixo
    // coarseModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/coarseJ_fem.msh","mirror_fine.txt", ".", 0);
    // fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/fineJ_fem.msh","mirror_fine.txt", ".", 0);
    
    // Stokes problem with analytic resolution
    // coarseModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/coarseJT_fem.msh","mirror_fine.txt",0);
    // fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/fineJT_fem.msh","mirror_fine.txt",0);

    
    //3D
    // fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/Laplace_fem3d1.msh","mirror_fine.txt",0);
    // coarseModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/LaplaceCoarse_fem3d1.msh","mirror_coarse.txt",0);

    //Stokes problem with analytic resolution
    // fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/LaplaceT_fem3d1.msh","mirror_fine.txt",0);
    // coarseModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/LaplaceCoarseT_fem3d1.msh","mirror_coarse.txt",0);

    
    //ArlequinProblem.setFluidModels_FEM_FEM(coarseModel,fineModel);
    // ArlequinProblem.solveArlequinProblemLaplace_FEM_FEM(10, 1.e-6);
    //ArlequinProblem.solveArlequinProblem_FEM_FEM(3, 1.e-6);


    //Finalize main program   
    PetscFinalize();
 
    return 0; 
};
 
  




 
