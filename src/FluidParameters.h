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
    void setFieldForce(double *ff){
        for (int i = 0; i < DIM; i++){
            fieldForce[i] = ff[i];
        }
    };

    // Sets arlequin parameters 
    void setArlequink1(double &arlk1) {arlequinK1 = arlk1;};
    void setArlequink2(double &arlk2) {arlequinK2 = arlk2;};
    void setArlequinEpsilon(double &arlqEps) {arlequinEpsilon = arlqEps;};
    void setGlueZoneThickness(double &thick) {glueZoneThickness = thick;};

    void setVelocityInf(double *velInf){
        for (int i = 0; i < DIM; i++){
            velocityInf[i] = velInf[i];
        };
    };
    
    void setNumTimeSteps(int &numTS) {numTimeSteps = numTS;};
    void setFreqPrint(int &pf) {printFreq = pf;};

    double& getTimeStep() {return timeStepSize;}
    double& getDensity() {return density;}
    double& getViscosity() {return viscosity;}
    double& getAlphaM() {return alpha_m;}
    double& getAlphaF() {return alpha_f;}
    double& getGamma() {return gamma;}
    double& getFieldForce(int dir) {return fieldForce[dir];}
    double& getArlequinK1() {return arlequinK1;};
    double& getArlequinK2() {return arlequinK2;};
    double& getArlequinEpsilon() {return arlequinEpsilon;};
    double& getGlueZoneThickness () {return glueZoneThickness;};
    double& getVelocityInf(int dir) {return velocityInf[dir];}
    double& getSpectralRadius() {return spectralRadius;};
    int& getNumTimeSteps() {return numTimeSteps;};
    int& getFreqPrint() {return printFreq;};


private:
    double viscosity;
    double density;
    double timeStepSize;
    double spectralRadius;
    double alpha_m;
    double alpha_f;
    double gamma;
    double fieldForce [DIM];
    double arlequinK1;
    double arlequinK2;
    double glueZoneThickness;
    double arlequinEpsilon;
    double velocityInf[DIM];
    int    numTimeSteps;
    int    printFreq;

    
};

#endif
