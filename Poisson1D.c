static char help[] = " 1D Finite Element Programm for the equation"
                "-p u'' + q u' + r u = f,x in (0,1)"
                "assume uniform unform mesh in [0,1]";

#include <stdio.h>
# include <math.h>
#include "Quadrature1D.h"
#include "UserParameters.h"
#include <petsc.h>

typedef struct {
    PetscInt    nParts, // number pf partitions. Assume uniform mesh in 1D.
                nNodes, // number of nodes
                dof, // degree of freedom
                quadOrder; // order of the quadrature
    // boundary data
    PetscBool   leftD,rightD,leftN,rghtN;
    // functions
    PetscReal   (*RHSfunct)(PetscReal); // right hand side
    PetscReal   (*exactSol)(PetscReal); // exact solution
    PetscReal   (*exactGradSol)(PetscReal); // exact solution
    PetscReal   (*gDfunct)(PetscReal);
    PetscReal   (*gNfunct)(PetscReal);
} paramsCtx;


extern PetscErrorCode FormRHS(Vec, paramsCtx* );
extern PetscErrorCode FormMatrix(Mat, paramsCtx* );
extern PetscErrorCode Errors(Vec, paramsCtx* );

int main(int argc,char **args){
    
    paramsCtx   user;
    Vec         fForce, uExact;
    Mat         K;
    KSP         ksp;
    
    PetscCall(PetscInitialize(&argc,&args,NULL,help));
    
    // - - DEFAULT VALUES
    user.quadOrder = 4; // order of the quadrature (2..6)
    user.nParts = 5; // number of partitions
    user.dof = 2; // quadratic
    
    // - - - - OPTIONS - - - -
    PetscOptionsBegin(PETSC_COMM_WORLD, "fem_", "options for 1D fem on Poisson equaiton", "");
    PetscCall(PetscOptionsInt("-parts", "number of partitions", "Poisson1D.c",user.nParts,&(user.nParts),NULL));
    PetscCall(PetscOptionsInt("-quadOrder", "quadrature degree (2, . . ,6)", "Poisson1D.c",user.quadOrder,&(user.quadOrder),NULL));
    PetscOptionsEnd();
    // - - - - END OPTIONS - - - -
    
    // Some parameter calculations and initializations
    user.nNodes = user.dof*user.nParts + 1; // number of nodes
    user.RHSfunct = &RHSfunct_one;
    user.exactSol = &ExactSolution;
    user.exactGradSol = &ExactGradSolution;
    user.gDfunct = &gD_one;
    user.leftD = PETSC_TRUE;
    user.rightD = PETSC_TRUE;
    user.leftN = PETSC_FALSE;
    user.rghtN = PETSC_FALSE;
    
   
    // - - - - RHS - - - -
    PetscCall(VecCreate(PETSC_COMM_WORLD,&fForce)); //fForce is the vector on RHS
    PetscCall(VecSetSizes(fForce,PETSC_DECIDE,user.nNodes));
    PetscCall(VecSetFromOptions(fForce));
    PetscCall(VecDuplicate(fForce,&uExact)); // duplicate vector for exact solution
    PetscCall(VecSet(fForce,0.0)); // initialize the vector to 0
    PetscCall(FormRHS(fForce,&user)); // calculate fForce
    
    // - - - - MATRIX - - - -
    PetscCall(MatCreate(PETSC_COMM_WORLD,&K)); // stiffness matrix
    PetscCall(MatSetSizes(K,PETSC_DECIDE,PETSC_DECIDE,user.nNodes,user.nNodes));
    PetscCall(MatSetFromOptions(K));
    PetscCall(MatSetOption(K,MAT_SYMMETRIC,PETSC_TRUE)); // CHECK with boundary cinditions
    PetscCall(MatSetUp(K));
    PetscCall(FormMatrix(K,&user));
    PetscCall(MatAssemblyBegin(K,MAT_FINAL_ASSEMBLY));
    PetscCall(MatAssemblyEnd(K,MAT_FINAL_ASSEMBLY));
    
    
    //  - - - - SOLVE linear system - - - -
    PetscCall(KSPCreate(PETSC_COMM_WORLD,&ksp));
    PetscCall(KSPSetOperators(ksp,K,K));
    PetscCall(KSPSetFromOptions(ksp));
    PetscCall(KSPSolve(ksp,fForce,uExact));
    
    // - - - VIEW - - -
    //PetscCall(VecView(uExact,PETSC_VIEWER_STDOUT_WORLD));
    
    // - - - L2 & H1 errors
    Errors(uExact,&user);
    
    PetscCall(VecDestroy(&fForce));
    PetscCall(VecDestroy(&uExact));
    PetscCall(MatDestroy(&K));
    PetscCall(PetscFinalize());
    return 0;
}

// = = = = = = FUNCTIONS = = = = = =

// - - - - START SHAPE FUNCTIONS - - - -
void QuadShape(PetscReal Basis[], PetscReal x, PetscReal x1, PetscReal x2) {
    PetscReal x3 = 0.5*(x1+x2); // middle element
    PetscReal h = x2 - x1;
    Basis[0] = 2.0*(x3-x)*(x2-x)/(h*h);
    Basis[1] = 4.0*(x-x1)*(x2-x)/(h*h);
    Basis[2] = 2.0*(x-x3)*(x-x1)/(h*h);
}

void QuadShapeGrad(PetscReal gradBasis[], PetscReal x, PetscReal x1, PetscReal x2) {
    PetscReal x3 = 0.5*(x1+x2); // middle element
    PetscReal h = x2 - x1;
    gradBasis[0] = 2*(2*x-x2-x3)/(h*h);
    gradBasis[1] = 4*(-2*x+x1+x2)/(h*h);
    gradBasis[2] = 2*(2*x-x1-x3)/(h*h);
}
// - - - - END SHAPE FUNCTIONS - - - -


PetscErrorCode FormRHS(Vec F, paramsCtx* ctx) {
    PetscReal   *auxF;
    PetscInt    leftNd;
    PetscReal   scale, shift, xGlobal, left, right;
    PetscReal   Basis[ctx->dof+1];
    Quad1D quad = GausseQuadrature[ctx->quadOrder-2]; // quadrature defined in "Quadrature1D.h".
    
    PetscCall(VecGetArray(F,&auxF)); // Alternatively 'VecSetValues' can be used
    
    scale = 0.5/ctx->nParts; // for quadrature
    // loop over the each element
    for (PetscInt ele = 0; ele < ctx->nParts; ele++){
        leftNd = ele*ctx->dof; // left nodal index of the element 'ele'
        left = ((PetscReal) ele)/ctx->nParts; // left coordinate of the element 'ele'
        right = ((PetscReal) ele+1)/ctx->nParts; // right coordinate of the element 'ele'
        shift = 0.5*(left + right);
        // quadrature
        for (PetscInt q = 0; q < ctx->quadOrder; q++ ){
            xGlobal = scale * quad.qPts[q] + shift;
            QuadShape(Basis,xGlobal,left,right); // quadratic shape functions
            // create the vector
            for (PetscInt n = 0; n < ctx->dof+1; n++){
                auxF[leftNd+n] += scale*quad.weghts[q]*Basis[n]*ctx->RHSfunct(xGlobal);
            } //end of n loop
        } // end of q loop
    } // end of ele loop
    auxF[0] = 0;
    auxF[ctx->nNodes-1] = 0;
    PetscCall(VecRestoreArray(F,&auxF));
    return 0;
}


// - - - - LHS
// Assemble local matrix
PetscErrorCode FormLocalMatrix(PetscReal localK[], PetscInt ele, PetscInt localNd[], paramsCtx* ctx){
    PetscInt    leftNd;
    PetscInt    ct = 0; // count for the local stiffness array (not the matrix)
    PetscReal   scale, shift, xGlobal, left, right;
    PetscReal   coeffs[3] = {0,0,0}; // vector of 3 coefficients
    PetscReal   Basis[ctx->dof+1];
    PetscReal   gradBasis[ctx->dof+1];
    Quad1D      quad = GausseQuadrature[ctx->quadOrder-2]; // quadrature defined in "Quadrature1D.h".
    
    leftNd = ele*ctx->dof; // left nodal index of the element 'ele'
    left = ((PetscReal) ele)/ctx->nParts; // left coordinate of the element 'ele'
    right = ((PetscReal) ele+1)/ctx->nParts; // right coordinate of the element 'ele'
    shift = 0.5*(left + right);
    scale = 0.5/ctx->nParts; // for quadrature
    for (PetscInt i = 0; i < ctx->dof + 1; i++){localNd[i] = leftNd +i;}
    for (PetscInt i = 0; i < (ctx->dof+1)*(ctx->dof+1); i++){localK[i] = 0;}// reset localK[dof+1] = 0
    
    // quadrature
    for (PetscInt q = 0; q < ctx->quadOrder; q++ ){
        xGlobal = scale * quad.qPts[q] + shift;
        Coefficients(coeffs,xGlobal); // evaluate coefficient functions at quadrature points
        QuadShape(Basis,xGlobal,left,right); // quadratic shape functions
        QuadShapeGrad(gradBasis,xGlobal,left,right); // gradiants of quadratic shape functions
        // create local matrix
        for (PetscInt r = 0; r < ctx->dof+1; r++){
            for (PetscInt c = 0; c < ctx->dof+1; c++){
                localK[ct] += scale*quad.weghts[q]*(coeffs[0]*gradBasis[r]*gradBasis[c] + coeffs[2]*Basis[r]*Basis[c]);
                ct++;
            } // end of the loop r
        } // end of the loop c
        ct = 0;
    } // end of the loop q
    return 0;
}

// Assemble gloablal matrix
PetscErrorCode FormMatrix(Mat K, paramsCtx* ctx){
    
    PetscInt    *nnz; //number of nonzeros in a row
    PetscReal   diagSum = 0.0; // sum of the diagonal elements
    PetscReal   localK[(ctx->dof + 1)*(ctx->dof + 1)]; // local stiffness matrix as an array

    PetscInt localNd[ctx->dof+1]; // local nodal indices
    
    // - - Finding number of nonzeros on each row
    // First and last and every other rows consists of 'dof+1' number of nonzeros
    // Other rows consists of '2*dof+1' number of nonzeros
    PetscCall(PetscMalloc1(ctx->nNodes,&nnz));
    for (PetscInt i = 0; i < ctx->nNodes-1; i+=2 ){
        nnz[i] = 2*ctx->dof + 1;
        nnz[i+1] = ctx->dof + 1;
    }
    nnz[0] = ctx->dof + 1;
    nnz[ctx->nNodes-1] = ctx->dof + 1;
    
    PetscCall(MatSeqAIJSetPreallocation(K,-1,nnz));
    PetscCall(PetscFree(nnz));
    // - - End of finding number of nonzeros on each row
    
    // loop over the internal elements
    for (PetscInt ele = 1; ele < ctx->nParts-1; ele++){
        PetscCall(FormLocalMatrix(localK,ele,localNd,ctx));
        // sum of diagonal entries for preconitioning
        for (PetscInt i = 0;i<ctx->dof + 1;i++){diagSum += localK[i*(ctx->dof+1)];}
        // create global matrix for the nonboundary elements
        PetscCall(MatSetValues(K,ctx->dof+1,localNd,ctx->dof+1,localNd,localK,ADD_VALUES));
    }
    // boundary elements
    // first element ele = 0
    PetscCall(FormLocalMatrix(localK,0,localNd,ctx));
    for (PetscInt i = 1;i<ctx->dof + 1;i++){
        localK[i*(ctx->dof+1)] = 0;
        localK[i] = 0;
    }
    localK[0] = diagSum/(ctx->nNodes-1); // average of diagonal elements
    PetscCall(MatSetValues(K,ctx->dof+1,localNd,ctx->dof+1,localNd,localK,ADD_VALUES));
    // last element  ele = ctx->nParts-1
    PetscCall(FormLocalMatrix(localK,ctx->nParts-1,localNd,ctx));
    for (PetscInt i = 0;i<ctx->dof;i++){
        localK[i*(ctx->dof+1)+ctx->dof] = 0;
        localK[(ctx->dof+1)*(ctx->dof+1)-i-2] = 0;
    }
    localK[(ctx->dof+1)*(ctx->dof+1)-1] = diagSum/(ctx->nNodes-1); // average of diagonal elements
    PetscCall(MatSetValues(K,ctx->dof+1,localNd,ctx->dof+1,localNd,localK,ADD_VALUES));
    
    return 0;
}


PetscErrorCode Errors(Vec uh, paramsCtx* ctx){
    PetscInt    leftNd;
    PetscInt    ix[ctx->dof+1];
    PetscReal   scale, shift, xGlobal, left, right;
    PetscReal   localAppx, localGradAppx;

    PetscReal   Basis[ctx->dof+1],localUh[ctx->dof+1];
    PetscReal   gradBasis[ctx->dof+1];
    PetscReal   l2Error = 0;
    PetscReal   h1Error = 0;
    
    Quad1D      quad = GausseQuadrature[ctx->quadOrder-2]; // quadrature defined in "Quadrature1D.h".
    
    for (PetscInt ele = 0; ele < ctx->nParts; ele++){
        leftNd = ele*ctx->dof; // left nodal index of the element 'ele'
        left = ((PetscReal) ele)/ctx->nParts; // left coordinate of the element 'ele'
        right = ((PetscReal) ele+1)/ctx->nParts; // right coordinate of the element 'ele'
        shift = 0.5*(left + right);
        scale = 0.5/ctx->nParts; // for quadrature
        // find approximated value at q in the element ele
        for (PetscInt i = 0; i < ctx->dof +1; i++){ix[i] = i+leftNd; } // indeces of element ele
        PetscCall(VecGetValues(uh, 3, ix, localUh));
        // quadrature
        for (PetscInt q = 0; q < ctx->quadOrder; q++ ){
            xGlobal = scale * quad.qPts[q] + shift;
            QuadShape(Basis,xGlobal,left,right); // quadratic shape functions
            QuadShapeGrad(gradBasis,xGlobal,left,right); // gradiants of quadratic shape functions
            localAppx = 0;
            localGradAppx = 0;
            for (PetscInt i = 0; i < ctx->dof +1; i++){
                localAppx += Basis[i]*localUh[i];
                localGradAppx += gradBasis[i]*localUh[i];
            } // end i
            l2Error += scale*quad.weghts[q]*(localAppx - ctx->exactSol(xGlobal))*(localAppx - ctx->exactSol(xGlobal));
            h1Error += scale*quad.weghts[q]*(localGradAppx - ctx->exactGradSol(xGlobal))*(localGradAppx - ctx->exactGradSol(xGlobal));
        } // end q
    }
    printf("\n L2 error: %f H1 error: %f \n", sqrt(l2Error), sqrt(h1Error));
    return 0;
}
