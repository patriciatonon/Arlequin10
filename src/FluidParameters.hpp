//------------------------------------------------------------------------------
// 
//                   Jeferson W D Fernandes and Rodolfo A K Sanches
//                             University of Sao Paulo
//                           (C) 2017 All Rights Reserved
//
// <LicenseText>
//------------------------------------------------------------------------------

//------------------------------------------------------------------------------
//--------------------------FLUID CONSTANT PARAMETERS---------------------------
//------------------------------------------------------------------------------

#ifndef FLUID_PARAMETERS_H
#define FLUID_PARAMETERS_H

using namespace boost::numeric;

/// Defines the fluid boundary shape functions

template<int DIM>
class FluidParameters {
public:
    
    // Sets the element viscosity
    void setViscosity(double& visc){viscosity = visc;}

    // Sets the element density
    void setDensity(double& dens){density = dens;}

    // Sets the time step size
    void setTimeStep(double& dt){timeStepSize = dt;}

    // Sets the time integration scheme
    void setSpectralRadius(double& b){
        spectralRadius = b;
        alpha_f = 1. / (1. + spectralRadius);
        alpha_m = 0.5 * (3. - spectralRadius) / (1. + spectralRadius);
        gamma = 0.5 + alpha_m - alpha_f;
    }

    // Sets body forces
    void setFieldForce(double *ff);

    // Sets arlequin parameters k1 and k2
    void setArlequink1(double &arlk1) {arlequinK1 = arlk1;};
    void setArlequink2(double &arlk2) {arlequinK2 = arlk2;};

    double& getTimeStep() {return timeStepSize;}
    double& getDensity() {return density;}
    double& getViscosity() {return viscosity;}
    double& getAlphaM() {return alpha_m;}
    double& getAlphaF() {return alpha_f;}
    double& getGamma() {return gamma;}
    double& getFieldForce(int dir) {return fieldForce[dir];}
    double& getArlequinK1() {return arlequinK1;};
    double& getArlequinK2() {return arlequinK2;};

private:
    double viscosity;
    double density;
    double timeStepSize;
    double spectralRadius;
    double alpha_m;
    double alpha_f;
    double gamma;
    double fieldForce [3];
    double arlequinK1;
    double arlequinK2;
    
};

template<>
void FluidParameters<2>::setFieldForce(double *ff){
    fieldForce[0] =  ff[0];
    fieldForce[1] =  ff[1];
}

template<>
void FluidParameters<3>::setFieldForce(double *ff){
    fieldForce[0] =  ff[0];
    fieldForce[1] =  ff[1];
    fieldForce[2] =  ff[2];
}





#endif
