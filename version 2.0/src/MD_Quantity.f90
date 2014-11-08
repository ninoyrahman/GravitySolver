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
    real(kind=double), dimension(0:inc*jnc*knc-1,-swidth:swidth,-swidth:swidth,-swidth:swidth) :: S

    ! Density
    real(kind=double), dimension(imin:imax,jmin:jmax,kmin:kmax) :: rho
    real(kind=double), dimension(imin:imax,jmin:jmax,kmin:kmax) :: rho4

    ! Error
    real(kind=double), dimension(imin:imax,jmin:jmax,kmin:kmax) :: error

    ! Coordinate
    real(kind=double), dimension(imin-ngh:imax+ngh) :: x
    real(kind=double), dimension(jmin-ngh:jmax+ngh) :: y
    real(kind=double), dimension(kmin-ngh:kmax+ngh) :: z

END MODULE MD_Quantity
