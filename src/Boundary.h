//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------------FLUID BOUNDARY--------------------------------
//------------------------------------------------------------------------------

#ifndef BOUNDARY_H
#define BOUNDARY_H

// Defines the fluid boundary object and its properties

template<int DIM>
class Boundary{

private:
    int          *connectB_;         //Boundary element connectivity
    int          index_;            //Boundary element index
    int          constrainType[3];  //Element type of constrain
    double       constrainValue[3]; //Element constrain value
    int          element_;          //Fluid Element
    int          elementSide_;      //Fluid Element Side
    int          group_;            //Element boundary group

public: 

    // Boundary element constructor
    Boundary(int *connec, int index,  \
             int *constrain, double *values, int gr){
        
        connectB_ = connec;
        index_ = index;
        group_ = gr;

        constrainType[0] = constrain[0];
        constrainType[1] = constrain[1];
        constrainType[2] = constrain[2];
        
        constrainValue[0] = values[0];
        constrainValue[1] = values[1];
        constrainValue[2] = values[2];

        element_ = 10000000;
        elementSide_ = 1000;
    };     

    // Returns the boundary element constrain component type
    int getConstrain(int dir){return constrainType[dir];}

    // Returns the boundary element constrain component value
    double getConstrainValue(int dir){return constrainValue[dir];}

    // Returns the boundary element connectivity
    int* getBoundaryConnectivity(){return connectB_;}

    // Sets the boundary element group
    void setBoundaryGroup(int gr){group_ = gr;}

    // Gets the boundary element group
    int getBoundaryGroup(){return group_;};

    // Sets the fluid element correspondence
    void setElement(int el){element_ = el;}

    // Sets the fluid element correspondence side at the boundary
    void setElementSide(int el){elementSide_ = el;}

    // Gets the fluid element correspondence
    int getElement(){return element_;}

    // Gets the fluid element correspondence side at the boundary
    int getElementSide(){return elementSide_;}

};






















#endif
