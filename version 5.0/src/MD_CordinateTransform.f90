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
    subroutine GetCartesianCordinate(i,j,k,xcord)
        implicit none
        integer, intent(in) :: i,j,k
        real(kind=double), dimension(0:2), intent(out) :: xcord

        xcord(0) = (dble(i-imin) + 0.5_rp)*h - 0.5_rp*dble(inc)*h
        xcord(1) = (dble(j-jmin) + 0.5_rp)*h - 0.5_rp*dble(jnc)*h
        xcord(2) = (dble(k-kmin) + 0.5_rp)*h - 0.5_rp*dble(knc)*h

    end subroutine GetCartesianCordinate

    ! nonuniformly scaled cartesian coordinate
    subroutine NUSCCordinate(i,j,k,xcord)
        implicit none
        integer, intent(in) :: i,j,k
        real(kind=double), dimension(0:2), intent(out) :: xcord

        call GetCartesianCordinate(i,j,k,xcord)
        xcord(0) = xcord(0)
        xcord(1) = 2.0_rp*xcord(1)
        xcord(2) = 2.0_rp*xcord(2)

    end subroutine NUSCCordinate


    subroutine GetCordinate(i,j,k,xcord)
        implicit none
        integer, intent(in) :: i,j,k
        real(kind=double), dimension(0:2), intent(out) :: xcord

        if(coordinate .eq. 0) then
            call GetCartesianCordinate(i,j,k,xcord)
        else if(coordinate .eq. 1) then
            call NUSCCordinate(i,j,k,xcord)
        else
            print *, 'wrong coordinate'
        end if

    end subroutine GetCordinate

END MODULE MD_CordinateTransform
