 !****************************************************************  
 ! MD_GeometricQuantity.f90
 ! 
 ! Created on: Jul 31, 2014
 ! Author: ninoy
 !***************************************************************** 

MODULE MD_GeometricQuantity
    USE MD_Parameter
    IMPLICIT NONE

     ! b matrix
    real(kind=double), dimension(imin-1:imax,jmin-1:jmax,kmin-1:kmax) :: bxx, bxy, bxz
    real(kind=double), dimension(imin-1:imax,jmin-1:jmax,kmin-1:kmax) :: byx, byy, byz
    real(kind=double), dimension(imin-1:imax,jmin-1:jmax,kmin-1:kmax) :: bzx, bzy, bzz

    ! G matrix
    real(kind=double), dimension(imin-1:imax,jmin-1:jmax,kmin-1:kmax) :: Gxy, Gxz
    real(kind=double), dimension(imin-1:imax,jmin-1:jmax,kmin-1:kmax) :: Gyx, Gyz
    real(kind=double), dimension(imin-1:imax,jmin-1:jmax,kmin-1:kmax) :: Gzx, Gzy


END MODULE MD_GeometricQuantity


