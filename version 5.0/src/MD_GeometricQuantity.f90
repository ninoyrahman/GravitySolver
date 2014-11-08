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
    real(kind=double), dimension(:,:,:), allocatable :: bxx, bxy, bxz
    real(kind=double), dimension(:,:,:), allocatable :: byx, byy, byz
    real(kind=double), dimension(:,:,:), allocatable :: bzx, bzy, bzz

    ! G matrix
    real(kind=double), dimension(:,:,:), allocatable :: Gxy, Gxz
    real(kind=double), dimension(:,:,:), allocatable :: Gyx, Gyz
    real(kind=double), dimension(:,:,:), allocatable :: Gzx, Gzy

END MODULE MD_GeometricQuantity


