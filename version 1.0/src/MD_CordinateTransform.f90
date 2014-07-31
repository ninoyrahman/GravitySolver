!****************************************************************
! MD_CordinateTransform.f90
!
! Created on: Jul 27, 2014
! Author: ninoy
!*****************************************************************
MODULE MD_CordinateTransform
    USE MD_Definition
    USE MD_Parameter
    IMPLICIT NONE

CONTAINS

    ! uniform cartesian coordinate
    subroutine GetCartesianCordinate(i,j,k,x)
        implicit none
        integer, intent(in) :: i,j,k
        real(kind=double), dimension(0:2), intent(out) :: x

        x(0) = (i-imin + 0.5)*h - 0.5*inc*h
        x(1) = (j-jmin + 0.5)*h - 0.5*jnc*h
        x(2) = (k-kmin + 0.5)*h - 0.5*knc*h

    end subroutine GetCartesianCordinate

    ! nonuniformly scaled cartesian coordinate
    subroutine NUSCCordinate(i,j,k,x)
        implicit none
        integer, intent(in) :: i,j,k
        real(kind=double), dimension(0:2), intent(out) :: x

        call GetCartesianCordinate(i,j,k,x)
        x(0) = x(0)
        x(1) = 2.0*x(1)
        x(2) = 2.0*x(2)

    end subroutine NUSCCordinate


    subroutine GetCordinate(i,j,k,x)
        implicit none
        integer, intent(in) :: i,j,k
        real(kind=double), dimension(0:2), intent(out) :: x

        if(coordinate .eq. 0) then
            call GetCartesianCordinate(i,j,k,x)
        else if(coordinate .eq. 1) then
            call NUSCCordinate(i,j,k,x)
        else
            print *, 'wrong coordinate'
        end if

    end subroutine GetCordinate

END MODULE MD_CordinateTransform
