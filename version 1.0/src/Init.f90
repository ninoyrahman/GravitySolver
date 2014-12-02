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
    use MD_CalculateCordinate
    use MD_IO
    implicit none
    PetscErrorCode ierr

    call CalculateCordinate

    call MatCreate(PETSC_COMM_WORLD,A,ierr)
    call MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,inc*jnc*knc,inc*jnc*knc,ierr)
    call MatSetType(A,MATSEQAIJ,ierr)
    call MatSeqAIJSetPreallocation(A,61,PETSC_NULL_INTEGER,ierr)
    call MatSetUp(A,ierr)

    call VecCreate(PETSC_COMM_WORLD,b,ierr)
    call VecSetType(b,VECSEQ,ierr)
    call VecSetSizes(b,PETSC_DECIDE,inc*jnc*knc,ierr)

    call VecCreate(PETSC_COMM_WORLD,phi,ierr)
    call VecSetType(phi,VECSEQ,ierr)
    call VecSetSizes(phi,PETSC_DECIDE,inc*jnc*knc,ierr)

    if(scenerio .eq. 0) then
        call IC_Vaccum
    else if(scenerio .eq. 1) then
        call IC_CentredSphere
    else if(scenerio .eq. 2) then
        call IC_OffCentredGaussian
    else if(scenerio .eq. 3) then
        call IC_UniformDensity
    else if(scenerio .eq. 4) then
        call IC_SinusoidalDensity
    else if(scenerio .eq. 5) then
        call IC_Poly6Density
    else if(scenerio .eq. 6) then
        call IC_GaussianDensity
    else if(scenerio .eq. 7) then
        call IC_Ploytrop1
    else if(scenerio .eq. 8) then
        call IC_Ploytrop5
    else if(scenerio .eq. 9) then
        call IC_CondensedSphere
    else if(scenerio .eq. 10) then
        call IC_OffCentredCondensed
    end if

    call ComputeTransverseGradient
    call CalculateStencil
    call ComputeMatrix

    call KSPCreate(PETSC_COMM_WORLD,ksp,ierr)
    call KSPSetOperators(ksp,A,A,SAME_PRECONDITIONER,ierr)
    call KSPSetType(ksp,KSPCG,ierr)
    call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)
    call KSPGetPC(ksp,pc,ierr)
    call PCSetType(pc,PCSOR,ierr)
    call KSPSetTolerances(ksp,1.0d-17,PETSC_DEFAULT_DOUBLE_PRECISION, &
        PETSC_DEFAULT_DOUBLE_PRECISION, PETSC_DEFAULT_INTEGER,ierr)

    call ViewInfo

end subroutine Init
