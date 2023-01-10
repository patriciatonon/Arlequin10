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


    // Recomendações que envio a mim mesma do futuro:

    // Temos vetores normais colocados manualmente quando a malha local é a fine_iso3;
    // Se você esquecer de tirar esses vetores quando for rodar outras malhas,
    // provavelmente você terá um erro no seu código;
    
    
    //data reading functions need two files from coarse or fine mesh:  
    //1- Fluid Flow Data
    //2- FEM mesh (function dataReading) or IGA mesh (function dataReadingIso)




    //EXEMPLOS ISO/ISO
    //Exemplo inicial, com uma faixa na regiao inferior da cavidade de elementos fine
    // coarseModel.dataReading_ISO("meshcoarse_data.txt","coarse_iso2.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_ISO("meshfine_data.txt","fine_iso.msh","mirror_fine.txt",0);

    //Exemplo tipo do Jeferson - Lembrar de colocar a normal no Set Signaled Distance
    // coarseModel.dataReading_ISO("meshcoarse_data.txt","coarse_iso2.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_ISO("meshfine_data.txt","fine_iso3.msh","mirror_fine.txt",0);


    //EXEMPLOS FEM/FEM 
    //Exemplo tipo do Jeferson da Cavidade 2D - Lembrar de colocar os valores de normais no Set Signaled Distance
    // coarseModel.dataReading_FEM("meshcoarse_data.txt","coarseT_fem.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("meshfine_data.txt","fineT_fem.msh","mirror_fine.txt",0);

    //Exemplo tipo do Jeferson com malha coarse não estruturada - Lembrar de colocar os valores de normais no Set Signaled Distance
    // coarseModel.dataReading_FEM("meshcoarse_data.txt","coarse_fem.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("meshfine_data.txt","fineT_fem.msh","mirror_fine.txt",0);


    //EXEMPLO CILINDRO
    // coarseModel.dataReading_FEM("meshcoarse_data.txt","coarse_fem_cylinder_jefe.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("meshfine_data.txt","fine_fem_cylinder_jefe.msh","mirror_fine.txt",0);

    //Exemplo inicial, com uma faixa na regiao inferior da cavidade de elementos fine
    // coarseModel.dataReading_FEM("meshcoarse_data.txt","coarseIni_fem.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("meshfine_data.txt","fineIni_fem.msh","mirror_fine.txt",0);
    // mais um
    // coarseModel.dataReading_FEM("meshcoarse_data.txt","coarseIni_fem.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("meshfine_data.txt","fine_fem.msh","mirror_fine.txt",0);
    
    // *****problema quando mais de 1 elemento coarse sobre um elemento fine - nao converge - FUNCIONA NA OUTRA VERSAO DO JEFERSON******
    // coarseModel.dataReading_FEM("meshcoarse_data.txt","coarse_fem.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("meshfine_data.txt","fine_fem.msh","mirror_fine.txt",0);



    //EXEMPLOS FEM/ISO (FINE/COARSE) 
    //Exemplo tipo Jeferson - Lembrar de colocar os valores de normais no Set Signaled Distance
    coarseModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/coarse_iso2.msh","mirror_coarse.txt",0);
    fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/fineT_fem.msh","mirror_fine.txt",0);

    //Exemplo inicial, com uma faixa na regiao inferior da cavidade de elementos fine
    // coarseModel.dataReading_ISO("../../mesh/meshfine_data.txt","../../mesh/coarse_iso2.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("../../mesh/meshfine_data.txt","../../mesh/fine_fem.msh","mirror_fine.txt",0);

    //EXEMPLO CILINDRO
    // coarseModel.dataReading_ISO("../../mesh/meshfinecyl_data.txt","../../mesh/coarse_iso_cylinder.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_FEM("../../mesh/meshfinecyl_data.txt","../../mesh/fine_fem_cylinder_jefe_nre.msh","mirror_fine.txt",0);


    //EXEMPLOS ISO/FEM
    //PROBLEMA 1 - Cavidade com uma faixa de elementos mais refinidos na parte inferior
    // coarseModel.dataReading_FEM("meshcoarse_data.txt","coarseIni_fem.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_ISO("meshfine_data.txt","fine_iso.msh","mirror_fine.txt",0);

    //PROBLEMA 1 - opção 2 - malha coarse não estruturada
	// coarseModel.dataReading_FEM("meshcoarse_data.txt","coarse_fem.msh","mirror_coarse.txt",0);
 //    fineModel.dataReading_ISO("meshfine_data.txt","fine_iso.msh","mirror_fine.txt",0);


    //PROBLEMA 2  - Cavidade com faixa ao redor de elementos refinados - Lembrar de colocar valores normais no Set Signaled Distance
    // coarseModel.dataReading_FEM("meshcoarse_data.txt","coarseIni_fem.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_ISO("meshfine_data.txt","fine_iso3.msh","mirror_fine.txt",0);
    //Problema 2 - opção 3 - malha não estruturada coarse
    // coarseModel.dataReading_FEM("meshcoarse_data.txt","coarse_fem.msh","mirror_coarse.txt",0);
    // fineModel.dataReading_ISO("meshfine_data.txt","fine_iso3.msh","mirror_fine.txt",0);

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
 
  




 
