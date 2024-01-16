#ifndef CASES_H_
#define CASES_H_

#include <petsc.h>

PetscReal ExactSolution(PetscReal x) {
    return PetscSinReal(4*PETSC_PI*x);
}

PetscReal ExactGradSolution(PetscReal x) {
    return  4*PETSC_PI*PetscCosReal(4*PETSC_PI*x);
}

PetscReal RHSfunct_one(PetscReal x) {
    return (-16*PETSC_PI*PETSC_PI +1)*PetscSinReal(4*PETSC_PI*x);
}


void Coefficients(PetscReal coeffs[], PetscReal x) {
    coeffs[0] = -1; //Coefficient p
    coeffs[1] = 0; //Coefficient q
    coeffs[2] = 1; //Coefficient r
}

//Dirichlet boundary condition
PetscReal gD_one(PetscReal x) {
    return ExactSolution(x);
}


PetscReal gN(PetscReal x) {
    return 0.0;
}

#endif
