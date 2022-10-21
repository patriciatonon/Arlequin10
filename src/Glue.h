//------------------------------------------------------------------------------
// 
//         Jeferson W D Fernandes,Rodolfo A K Sanches and Patricia Tonon
//                             University of Sao Paulo
//                           (C) 2021 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------
 
//------------------------------------------------------------------------------
//----------------------------------GLUE ZONE-----------------------------------
//------------------------------------------------------------------------------

#ifndef GLUE_H
#define GLUE_H

#include "Node.h"

template<int DIM>
class Glue{
private:
    int *connect_;          //New connectivity from the gluing zone 
    int *newconnect_;          //New connectivity from the gluing zone          
    int index_;             //Element Number in the gluing zone
    int elemCorrespondent_; //Correspondent number in the fluid mesh
 

public:
    // Gluing Zone Element constructor
    Glue(int index, int elemCorrespondent){
        index_ = index;
        elemCorrespondent_ = elemCorrespondent;
    };

    // Returns the fluid element correspondence
    int getElemCorrespondent(){return elemCorrespondent_;};

    // Sets the element connectivity
    void setConnectivity(int *connect){connect_ = connect;};

    // Sets the new element connectivity
    void setNewConnectivity(int *connect){newconnect_ = connect;};

    // Gets the element connectivity
    int *getConnectivity(){return connect_;};

    // Gets the new element connectivity
    int *getNewConnectivity(){return newconnect_;};
    
};


//------------------------------------------------------------------------------
//--------------------------------IMPLEMENTATION--------------------------------
//------------------------------------------------------------------------------


#endif

