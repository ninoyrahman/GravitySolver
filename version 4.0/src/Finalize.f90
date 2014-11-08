!****************************************************************
! Finalize.f90
!
! Created on: Jul 26, 2014
! Author: ninoy
!*****************************************************************

subroutine Finalize
#include <finclude/petscdef.h>
    use petscksp
    use MD_Quantity
    implicit none

    PetscErrorCode ierr

    call MatDestroy(A,ierr)
    call VecDestroy(b,ierr)
    call VecDestroy(phi,ierr)
    call KSPDestroy(ksp,ierr)

end subroutine Finalize
