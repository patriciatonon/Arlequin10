#include "FluidParameters.h"

template<>
void FluidParameters<2>::setFieldForce(double *ff){
    fieldForce[0] =  ff[0];
    fieldForce[1] =  ff[1];
};

template<>
void FluidParameters<3>::setFieldForce(double *ff){
    fieldForce[0] =  ff[0];
    fieldForce[1] =  ff[1];
    fieldForce[2] =  ff[2];
};

