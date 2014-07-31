 !****************************************************************  
 ! MD_Quantity.f90
 ! 
 ! Created on: Jul 23, 2014
 ! Author: ninoy
 !***************************************************************** 

MODULE MD_Quantity
#include <finclude/petscdef.h>
    USE petscksp
    USE MD_Definition
    USE MD_Parameter
    IMPLICIT NONE
    SAVE

    ! PETSc objects for A*phi=b, phi = gravitational potential
    Mat A
    Vec b
    Vec phi
    KSP ksp
    PC pc

    ! Stencil
    real(kind=double), dimension(-swidth:swidth,-swidth:swidth,-swidth:swidth) :: S

    ! Density
    real(kind=double), dimension(imin:imax,jmin:jmax,kmin:kmax) :: rho

END MODULE MD_Quantity
