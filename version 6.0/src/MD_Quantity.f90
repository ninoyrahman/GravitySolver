 !****************************************************************  
 ! MD_Quantity.f90
 ! 
 ! Created on: Jul 23, 2014
 ! Author: ninoy
 !***************************************************************** 

MODULE MD_Quantity
    USE MD_Definition
    USE MD_Parameter
    IMPLICIT NONE
    SAVE

    ! Stencil
    real(kind=double), dimension(:,:,:,:), allocatable :: S

    ! Density
    real(kind=double), dimension(:,:,:), allocatable :: rho
    real(kind=double), dimension(:,:,:), allocatable :: rho4

    ! Global and Local Index
    integer, dimension(:,:,:), allocatable :: globalindex
    integer, dimension(:,:,:), allocatable :: localindex

    ! Coordinate
    real(kind=double), dimension(:), allocatable :: x
    real(kind=double), dimension(:), allocatable :: y
    real(kind=double), dimension(:), allocatable :: z

END MODULE MD_Quantity
