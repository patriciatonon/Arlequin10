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
#include "src/FluidData.hpp"
#include "src/Arlequin.hpp" 


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
    typedef Arlequin<dimension>  Arlequin;
 
    //Create problem variables 
    FluidModel coarseModel,fineModel; 
    Arlequin   ArlequinProblem; 

//================================================================================
//=================================PROBLEM MESH===================================
//================================================================================

    MPI_Barrier(PETSC_COMM_WORLD);


    // Recomendações que envio a mim mesma do futuro:

    // Temos vetores normais colocados manualmente quando a malha local é a fine_iso3;
    // Se você esquecer de tirar esses vetores quando for rodar outras malhas,
    // provavelmente você terá um erro no seu código;
    
    
    //data reading functions need two files from coarse or fine mesh:  
    //1- Fluid Flow Data
    //2- FEM mesh (function dataReading) or IGA mesh (function dataReadingIso)
    coarseModel.dataReading_ISO("meshcoarse_data.txt","coarse_iso.msh","mirror_coarse.txt",0);
    //fineModel.dataReading_ISO("meshfine_data.txt","fine_iso.msh","mirror_fine.txt",0);
    fineModel.dataReading_FEM("meshfine_data.txt","fine_fem.msh","mirror_fine.txt",0);
    
    //ArlequinProblem needs two objects from fluid:
    //1 - coarseModel;
    //2 - fineModel;
    ArlequinProblem.setFluidModels(coarseModel,fineModel);
    
    //solveArlequinProblem function needs two parameters:  
    //1- The maximum number of iterations in the Newton-Raphson process
    //2- The maximum relative error in the Newton-Raphson process (DU) 
    ArlequinProblem.solveArlequinProblem(10, 1.e-6);
    
    //Finalize main program   
    PetscFinalize();
 
    return 0; 
};
 
  




 
