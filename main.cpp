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

    //EXEMPLOS FEM/ISO (FINE/COARSE) 
    //Exemplo tipo Jeferson - Lembrar de colocar os valores de normais no Set Signaled Distance
    // coarseModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/coarse_iso2.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/fineT_fem.msh","mirror_fine.txt",0);

    //Exemplo inicial, com uma faixa na regiao inferior da cavidade de elementos fine
    // coarseModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/coarse_iso2.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/fine_fem.msh","mirror_fine.txt",0);

    //EXEMPLO CILINDRO
    // coarseModel.dataReading_ISO("../../mesh/meshfinecyl_data.txt","../../mesh/coarse_iso_cylinder.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("../../mesh/meshfinecyl_data.txt","../../mesh/fine_fem_cylinder_jefe_nre.msh","mirror_fine.txt",0);


    //3D problem cavity 
    // coarseModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/NEWcavity_iso3d.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/NEWcavity_fem3dEST.msh","mirror_fine.txt",0);

    //3D LAPLACE
    coarseModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/Laplace_iso3d.msh","mirror_coarse.txt",0);
    fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/Laplace_fem3d.msh","mirror_fine.txt",0);


    //ArlequinProblem needs two objects from fluid:
    //1 - coarseModel;
    //2 - fineModel;
    ArlequinProblem.setFluidModels_FEM_ISO(coarseModel,fineModel);
    
    //solveArlequinProblem function needs two parameters:  
    //1- The maximum number of iterations in the Newton-Raphson process
    //2- The maximum relative error in the Newton-Raphson process (DU) 
    //ArlequinProblem.solveArlequinProblem_FEM_ISO(10, 1.e-6);

    ArlequinProblem.solveArlequinProblemLaplace_FEM_ISO(10, 1.e-6);
    
    //Finalize main program   
    PetscFinalize();
 
    return 0; 
};
 
  




 
