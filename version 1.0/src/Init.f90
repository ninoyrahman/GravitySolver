!****************************************************************
! Init.f90
!
! Created on: Jul 23, 2014
! Author: ninoy
!*****************************************************************

subroutine Init
#include <finclude/petscdef.h>
    use petscksp
    use MD_Parameter
    use MD_Quantity
    use MD_IO
    implicit none
    PetscErrorCode ierr

    call MatCreate(PETSC_COMM_WORLD,A,ierr)
    call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,inc*jnc*knc,inc*jnc*knc,ierr)
    call MatSetType(A,MATSEQAIJ,ierr)
    call MatSeqAIJSetPreallocation(A,(2*swidth+1)*(2*swidth+1)*(2*swidth+1),PETSC_NULL_INTEGER,ierr)
    call MatSetUp(A,ierr)

    call VecCreate(PETSC_COMM_WORLD,b,ierr)
    call VecSetType(b,VECSEQ,ierr)
    call VecSetSizes(b,PETSC_DECIDE,inc*jnc*knc,ierr)

    call VecCreate(PETSC_COMM_WORLD,phi,ierr)
    call VecSetType(phi,VECSEQ,ierr)
    call VecSetSizes(phi,PETSC_DECIDE,inc*jnc*knc,ierr)

    call ComputeTransverseGradient
    call CalculateStencil(imin,jmin,kmin)
    call ComputeMatrix

    call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
    call KSPSetOperators(ksp,A,A,DIFFERENT_NONZERO_PATTERN,ierr)
    call KSPSetType(ksp,KSPCG,ierr)
    call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)
    call KSPGetPC(ksp,pc,ierr)
    call PCSetType(pc,PCILU,ierr)
    call KSPSetTolerances(ksp,1.e-6,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION, &
        PETSC_DEFAULT_INTEGER,ierr)

    if(scenerio .eq. 0) then
        call IC_Vaccum
    else if(scenerio .eq. 1) then
        call IC_CentredSphere
    else if(scenerio .eq. 2) then
        call IC_OffCentredSphere
    else if(scenerio .eq. 3) then
        call IC_UniformDensity
    else if(scenerio .eq. 4) then
        call IC_SinusoidalDensity
    else if(scenerio .eq. 5) then
        call IC_Poly6Density
    end if

    call ViewInfo

end subroutine Init
