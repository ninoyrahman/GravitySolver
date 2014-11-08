 !****************************************************************  
 ! MD_Helper.f90
 ! 
 ! Created on: Jul 24, 2014
 ! Author: ninoy
 !***************************************************************** 

MODULE MD_Helper
#include <finclude/petscdef.h>
    USE petscksp
    USE MD_Parameter
    IMPLICIT NONE

CONTAINS

    !***************************************************************
    ! returns 1 if cell(i,j,k) is a boundary cell
    ! returns 0 if cell(i,j,k) is a inner cell
    !***************************************************************
    integer function IsBoundary(i,j,k)
        implicit none
        integer, intent(in) :: i,j,k
        if(i .lt. imin .or. i .gt. imax .or. &
            j .lt. jmin .or. j .gt. jmax .or. &
            k .lt. kmin .or. k .gt. kmax) then
            IsBoundary = 1
        else
            IsBoundary = 0
        end if
    end function IsBoundary

    !***************************************************************
    ! this funtion convert multidimensional index(i,j,k)
    ! of inner cell to linear index starting from zero.
    ! cell index: (2,2,3)-->1
    ! with 2 boundary cell in each dimension
    !***************************************************************
    integer function GetIndex(i,j,k)
        implicit none
        integer, intent(in) :: i,j,k
        GetIndex = (i-imin) + inc*(j-jmin) + inc*jnc*(k-kmin)
    end function GetIndex

    !***************************************************************
    ! this funtion convert multidimensional index(row,column)
    ! of matrix entry to linear index starting from zero.
    ! cell index: (1,1)-->3
    ! with 2 column in each row
    !***************************************************************
    integer function GetLinearIndex(row,column)
        implicit none
        integer, intent(in) :: row,column
        GetLinearIndex = column + row*inc*jnc*knc
    end function GetLinearIndex

END MODULE MD_Helper

