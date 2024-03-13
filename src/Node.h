//------------------------------------------------------------------------------
// 
//        Jeferson W D Fernandes,Rodolfo A K Sanches and Patrícia Tonon
//                             University of Sao Paulo
//                           (C) 2019 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//------------------------------------NODES-------------------------------------
//------------------------------------------------------------------------------

#ifndef NODE_H
#define NODE_H

#include "DataTypes.h"
#include<cstdlib>
#include<fstream>
#include<iostream>
#include <math.h>       /* atan */
#include <vector>

// Defines the node object and stores all nodal variables information

template<int DIM>

class Node{
private:
    //Main variables
    double           coord_[DIM];                   //Nodal coordinate vector
    double           coordUpdated_[DIM];            //Updated nodal coordinate vector
    double           initialCoord_[DIM];            //Initial nodal coordinate vector
    double           previousCoord_[DIM];           //Previous nodal coordinate vector
    double           velocity_[DIM];                //Nodal velocity
    double           previousVelocity_[DIM];        //Previous time step velocity
    double           acceleration_[DIM];            //Nodal acceleration
    double           previousAcceleration_[DIM];    //Previous time step acceleration
    double           pressure_;                     //Nodal pressure 
    double           meshVelocity_[DIM];            //Nodal mesh velocity
    double           previousMeshVelocity_[DIM];    //Previous time step mesh velocity  
    double           Force_[DIM];                   //Nodal force in Nodes
             
    int              constrainType[3];              //Constrain direction
    double           constrainValue[3];             //Nodal prescribed velocity
    std::vector<int> invIncidence;					//Elements that have this node

    //IGA mesh variables
    int              INC_[DIM];                     //Base Function IA index
    int              sizeMcp_;                      //Size vector matchingCP
    double           weightPC_;                     //Weight from Control Points    
    int              *matchingCP_;                  //Vector that stores the repeated control points in patches
    int              newcon_;                       //Connectivity between a global number with
                                                    //a new global numeration that discarts repeated CP  
    //Arlequin
    double           weightFunction_;               //Nodal Energy Weight Function
    double           previousWeightFunction_;       //Previous nodal Energy Weight Function
    double           nNodal_[DIM];                  //Nodal normal vector
    double           distGlueZone;                  //Nodes signaled distance to a defined boundary in gluing zone
    bool             glueZone;                      //Defines if the node is in the gluing zone (true = gluing zone)  
    int              elemCorresp;                   //Element correspondence
    double           xsiCorresp[DIM];               //Adimensional coord correspond
    double           lagMultiplier_[DIM];           //Nodal lagrange multiplier
    double           presArlequin_;                 //Glue zone pressure    
    double           velArlequin_[DIM];             //Glue zone velocity
    double           globalVelocity_[DIM];          //Global velocity for Arlequin stabilization problem

    double           analsol_[DIM];
    double           analpress_;

public:
    // Constructor - Defines a node with index,coordinates,weight and beloging patch
    Node(double *coor, int index, double wei){
        
        for (int i = 0; i < DIM; i++) {
        	
            coord_[i] = coor[i];
            coordUpdated_[i] = coor[i];
        	initialCoord_[i] = coor[i];
        	previousCoord_[i] = coor[i];
        	velocity_[i] = 0.0;
        	previousVelocity_[i] = 0.0;
        	acceleration_[i] = 0.0;
        	previousAcceleration_[i] = 0.0;
            meshVelocity_[i] = 0.0;
            previousMeshVelocity_[i] = 0.0;

            lagMultiplier_[i]= 0.;
            velArlequin_[i] = 0.;

            nNodal_[i] = 0.;
            xsiCorresp[i] = 0.;
        	
        	constrainType[i] = 0;
        	constrainValue[i] = 0.0;
            Force_[i] = 0.0;
        
        };

        pressure_ = 0.;
        invIncidence.clear();
        weightPC_ = wei;

        constrainType[DIM+1] = 0;
        constrainValue[DIM+1] = 0.;
        
        distGlueZone = 0.;
        glueZone = false;     
        elemCorresp = 0;
        weightFunction_ = 0.;
        previousWeightFunction_ = 0.;
        presArlequin_=0.;
    };

    //............................Coordinates functions............................
    void setCoordinates(double *u) {for (int i = 0; i <DIM; i++) coord_[i] = u[i];};

    // Returns the node initial coordinate vector
    double*  getInitialCoordinates() {return initialCoord_;};
    double   getInitialCoordinateValue(int dir) {return initialCoord_[dir];};

    // Returns the node coordinate vector
    double* getCoordinates() {return coord_;};
    double getCoordinateValue(int dir) {return coord_[dir];};

    // Increment the coordinate vector
    void incrementCoordinate(int dir, double u){coord_[dir] += u;};
   
    // Sets the previous coordinate vector
    void setPreviousCoordinates(double *u){for (int i = 0; i <DIM; i++) previousCoord_[i] = u[i];};

    // Returns the node coordinate vector at the previous time step
    double* getPreviousCoordinates() {return previousCoord_;}
    double getPreviousCoordinateValue(int dir) {return previousCoord_[dir];};

    // Sets the updated coordinate vector
    void setUpdatedCoordinates(int dir, double &u) {coordUpdated_[dir] = u;};
   
    // Returns the node updated coordinate vector
    double getUpdatedCoordinates(int dir) {return coordUpdated_[dir];};

    //............................IGA functions............................
    // Returns the weight of de control points
    double getWeightPC() {return weightPC_;};
    
    // Sets the base function indexes
    void setINC(int *inc){for (int i=0; i < DIM; i++) INC_[i] = inc[i];};

    // Returns the base function index
    int* getINC() {return INC_;};

    // Sets the new connectivity to solve the system
    void setnewcon(int newcon) {newcon_ = newcon;}

    // Returns the new connectivity to solve the system
    int getnewcon() {return newcon_;}

    // Sets to each control point the control points that have the same coordinates;
    void setMatchingCP(int *matchingCP) {matchingCP_ = matchingCP;}
    void setSizeMcp(int sizeMcp) {sizeMcp_ = sizeMcp;};
    
    // Returns to each control point the control points that have the same coordinates;
    int* getMatchingCP(){return matchingCP_;};
    
    //Returns the size of the vector matchingCP
    int  getSizeMcp(){return sizeMcp_;};


    //.................Elements respective to the Node......................
    // Pushs back a term of the inverse incidence, i.e., an element which
    // contains the node
    void pushInverseIncidence(int el) {invIncidence.push_back(el);}

    // Gets the number of elements which contains the node
    int getNumberOfElements(){return invIncidence.size();}

    // Gets an specific member of the inverse incidence
    int getInverseIncidenceElement(int i){return invIncidence[i];}
    
    
    //............................Velocity functions............................
    // Sets the velocity vector
    void setAnalSol(double *analsol, double &analpress) {
        for (int i=0; i < DIM; i++) analsol_[i] = analsol[i];
        analpress_ = analpress;
    };

    double getAnalSol(int dir) {return analsol_[dir];};
    double getAnalPress() {return analpress_;};

    void setVelocity(double *u) {for (int i=0; i < DIM; i++) velocity_[i] = u[i];};
    void setVelocityComponent(int dir, double val){velocity_[dir] = val;};

    // Returns the node velocity vector
    double getVelocity(int dir) {return velocity_[dir];}

    // Increment the velocity vector
    void incrementVelocity(int dir, double u) {velocity_[dir] += u;};

    // Sets the previous velocity vector
    void setPreviousVelocity(double *u) {for (int i=0; i < DIM; i++) previousVelocity_[i] = u[i];};
	void setPreviousVelocityComponent(int dir, double val){previousVelocity_[dir]= val;};
    
    // Returns the node previous time step velocity vector
    double getPreviousVelocity(int dir) {return previousVelocity_[dir];}


    //..........................Acceleration functions..........................
    // Sets the acceleration vector
    void setAcceleration(double *u){for (int i=0; i < DIM; i++) acceleration_[i] = u[i];};
    void setAccelerationComponent(int dir, double val){acceleration_[dir] = val;};

    // Gets the acceleration vectore
    double getAcceleration(int dir) {return acceleration_[dir];}

    // Increment the acceleration vector
    void incrementAcceleration(int dir, double u) {acceleration_[dir] += u;};  

    // Sets the previous time step acceleration vector
    void setPreviousAcceleration(double *u){for (int i=0; i < DIM; i++) previousAcceleration_[i] = u[i];};

    // Gets the previous time step acceleration vector
    double getPreviousAcceleration(int dir) {return previousAcceleration_[dir];}

    
    //............................Pressure functions............................
    // Sets the nodal pressure
    void setPressure(double p){pressure_ = p;};

    // Increments the nodal pressure
    void incrementPressure(double p){pressure_ += p;};

    // Gets the nodal pressure value
    double getPressure() {return pressure_;};


    //.........................Mesh Velocity functions..........................
    // Sets the node mesh velocity
    void setMeshVelocity(double *u){for (int i=0; i < DIM; i++) meshVelocity_[i] = u[i];}
    void setMeshVelocityValue(int dir, double &u){meshVelocity_[dir] = u;}

    // Gets the node mesh velocity
    double getMeshVelocity(int dir) {return meshVelocity_[dir];}

    // Gets the previous time step mesh velocity
    void setPreviousMeshVelocity(double *u){for (int i=0; i < DIM; i++) previousMeshVelocity_[i] = u[i];}

    // Gets the previous time step mesh velocity
    double getPreviousMeshVelocity(int dir) {return previousMeshVelocity_[dir];}

    
    //...........................Constrains functions...........................
    // Sets all node constrains     
    void setConstrains(int dir, int type, double value){
        constrainType[dir] = type;
        constrainValue[dir] = value;
        velocity_[dir] = value;
    };

    // Gets node constrain type
    int getConstrains(int dir) {return constrainType[dir];}

    void setConstrainValue(int dir, double value){
        constrainValue[dir] = value;
        velocity_[dir] = value;
    }

    // Gets node constrain value
    double getConstrainValue(int dir) {return constrainValue[dir];}
  

    //.......................Signaled distance functions........................
    // Sets nodal normal vector
    void setInnerNormal(double *n){ for (int i =0; i < DIM; i++) nNodal_[i] = n[i];};

    // Gets the nodal normal vector
    double* getInnerNormal() {return nNodal_;};

    // Clears nodal normal vector
    void clearInnerNormal(){for (int i = 0; i < DIM; i++) nNodal_[i] = 0.;};

    // Sets the Signaled distance function value
    void setDistFunction(double dist) {distGlueZone = dist;};

    // Gets the Signaled distance function value
    double getDistFunction(){return distGlueZone;};
   

    //............................Gluing Zone....................................
    // Sets if the Node is part of the blend zone
    void setGlueZone() {glueZone = true;};

    // Gets if the Node is part of the blend zone
    bool getGlueZone(){return glueZone;};

    // Sets nodal correspondence of overlapped mesh
    void setNodalCorrespondence(double elem, double* xsi){
    elemCorresp = elem;
    for (int i=0; i<DIM; i++) xsiCorresp[i] = xsi[i];};

    // Gets nodal correspondence of overlapped mesh - element
    int getNodalElemCorrespondence() {return elemCorresp;}

    // Gets nodal correspondence of overlapped mesh - adim. coordinate
    double* getNodalXsiCorrespondence() {return xsiCorresp;}


    //............................Arlequin functions............................
    // Sets the nodal energy weight function value
    void setWeightFunction(double val) {previousWeightFunction_ = weightFunction_; weightFunction_ = val;};

    // Gets the nodal energy weight function value
    double getWeightFunction() {return weightFunction_;};

    // Gets the nodal previous energy weight function value
    double getPreviousWeightFunction() {return previousWeightFunction_;};

    // Sets Lagrange Multiplier value
    void setLagrangeMultiplier(int dir, double lMult){
        lagMultiplier_[dir] = lMult;};

    // Increment the Lagrange Multiplier vector
    void incrementLagrangeMultiplier(int dir, double lMult){
        lagMultiplier_[dir] += lMult;};

    // Gets Lagrange Multiplier component value
    double getLagrangeMultiplier(int dir) {return lagMultiplier_[dir];};

    // Sets the interpolated Arlequin pressure 
    void setPressureArlequin(double p) {presArlequin_ = p;};

    // Gets interpolated Arlequin pressure
    double getPressureArlequin() {return presArlequin_;};

    // Sets the interpolated Arlequin velocity
    void setVelocityArlequin(int dir, double v) {velArlequin_[dir] = v;};

    // Gets the interpolated Arlequin velocity component
    double getVelocityArlequin(int dir) {return velArlequin_[dir];};

    // Sets global velocity for the Arlequin stabilization paremeter computation
    void setGlobalVelocity(double *u) {for (int i=0; i < DIM; i++) globalVelocity_[i] = u[i];};

    // Returns global velocity for the Arlequin stabilization paremeter computation
    double getGlobalVelocity(int dir) {return globalVelocity_[dir];}

    //Nodal force
    void setNodalForce(double *force) {for (int i=0; i < DIM; i++) Force_[i] = force[i];};
    
    double getNodalForce(int dir) {return Force_[dir];};

};


#endif

