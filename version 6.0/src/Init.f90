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
    use MD_PETScQuantity
    use MD_IO
    use MD_CalculateCordinate
    use MD_Helper
    implicit none

    PetscErrorCode ierr
    PetscInt low, high

    if(rank .eq. master) then
        print *, ''
        print *, '***************'
        print *, 'Initialization'
        print *, '***************'
        print *, ''
    endif

    h = domainsize/dble(inc)

    call DomainDecomposition

    call ComputeIndex

    call CalculateCordinate

    call VecCreate(cart_comm,b,ierr)
    call VecSetType(b,VECMPI,ierr)
    call VecSetSizes(b,ilnc*jlnc*klnc,inc*jnc*knc,ierr)

    call VecDuplicate(b, phi, ierr)


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
        call IC_CondensedSphere
    else if(scenerio .eq. 8) then
        call IC_OffCentredCondensed
    else if(scenerio .eq. 9) then
        call IC_Ploytrop1
    else if(scenerio .eq. 10) then
        call IC_Ploytrop5
    end if

    if(coordinate .eq. 0) then
        call ComputeCartesianMapping
    endif

    call ComputeTransverseGradient
    call ComputeStencil
    call ComputeMatrix

    call KSPCreate(cart_comm,ksp,ierr)
    call KSPSetOperators(ksp,A,A,SAME_PRECONDITIONER,ierr)
    call KSPSetType(ksp,KSPCG,ierr)
    call KSPSetInitialGuessNonzero(ksp,PETSC_TRUE,ierr)
    call KSPGetPC(ksp,pc,ierr)
    call PCSetType(pc,PCSOR,ierr)
    call KSPSetTolerances(ksp,1.0d-17,PETSC_DEFAULT_DOUBLE_PRECISION,PETSC_DEFAULT_DOUBLE_PRECISION, &
        PETSC_DEFAULT_INTEGER,ierr)

    if(rank .eq. master) then
        call ViewInfo
    end if

end subroutine Init
