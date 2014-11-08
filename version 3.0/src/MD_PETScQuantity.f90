!****************************************************************
! MD_PETScQuantity.f90
!
! Created on: Oct 28, 2014
! Author: ninoy
!*****************************************************************

MODULE MD_PETScQuantity
#include <finclude/petscdef.h>
    USE petscksp
    USE petscdm
    IMPLICIT NONE
    SAVE

    ! PETSc objects for A*phi=b, phi = gravitational potential
    Mat A
    Vec b
    Vec phi
    KSP ksp
    PC pc
    DM da

    ! MPI communicator
    MPI_Comm :: cart_comm


END MODULE MD_PETScQuantity

