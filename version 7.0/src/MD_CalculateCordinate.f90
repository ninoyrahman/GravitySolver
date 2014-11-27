!****************************************************************
! MD_CalculateCordinate.f90
!
! Created on: Sep 29, 2014
! Author: ninoy
!*****************************************************************

MODULE MD_CalculateCordinate
    USE MD_Definition
    USE MD_Parameter
    USE MD_Quantity
    IMPLICIT NONE

CONTAINS

    ! uniform cartesian coordinate
    subroutine CalculateCartesianCordinate
        implicit none

        integer :: i,j,k

        !$OMP PARALLEL PRIVATE(i)
        !$OMP DO SCHEDULE(STATIC)
        do i = imin-ngh, imax+ngh
            x(i) = (dble(i-ngh) + 0.5_rp)*h - 0.5_rp*dble(inc)*h + dble(igmin-ngh)*h
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        !$OMP PARALLEL PRIVATE(j)
        !$OMP DO SCHEDULE(STATIC)
        do j = jmin-ngh, jmax+ngh
            y(j) = (dble(j-ngh) + 0.5_rp)*h - 0.5_rp*dble(jnc)*h  + dble(jgmin-ngh)*h
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        !$OMP PARALLEL PRIVATE(k)
        !$OMP DO SCHEDULE(STATIC)
        do k = kmin-ngh, kmax+ngh
            z(k) = (dble(k-ngh) + 0.5_rp)*h - 0.5_rp*dble(knc)*h  + dble(kgmin-ngh)*h
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine CalculateCartesianCordinate

    ! nonuniformly scaled cartesian coordinate
    subroutine CalculateNUSCCordinate
        implicit none
        integer :: i,j,k

        call CalculateCartesianCordinate

        !$OMP PARALLEL PRIVATE(i)
        !$OMP DO SCHEDULE(STATIC)
        do i = imin-ngh, imax+ngh
            x(i) = x(i)
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        !$OMP PARALLEL PRIVATE(j)
        !$OMP DO SCHEDULE(STATIC)
        do j = jmin-ngh, jmax+ngh
            y(j) = 2.0_rp*y(j)
        end do
        !$OMP END DO
        !$OMP END PARALLEL

        !$OMP PARALLEL PRIVATE(k)
        !$OMP DO SCHEDULE(STATIC)
        do k = kmin-ngh, kmax+ngh
            z(k) = 2.0_rp*z(k)
        end do
        !$OMP END DO
        !$OMP END PARALLEL

    end subroutine CalculateNUSCCordinate


    subroutine CalculateCordinate
        implicit none

        if(coordinate .eq. 0) then
            call CalculateCartesianCordinate
        else if(coordinate .eq. 1) then
            call CalculateNUSCCordinate
        else
            print *, 'wrong coordinate'
        end if

    end subroutine CalculateCordinate


END MODULE MD_CalculateCordinate

