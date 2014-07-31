!****************************************************************
! IC_Vaccum.f90
!
! Created on: Jul 26, 2014
! Author: ninoy
!*****************************************************************

subroutine IC_Vaccum
    use MD_Parameter
    use MD_Quantity
    implicit none

    integer i,j,k

    print *, ''
    print *, '******************************'
    print *, 'Initializing Scenerio: Vaccum'
    print *, '******************************'
    print *, ''

    xc(0) = 0.0
    xc(1) = 0.0
    xc(2) = 0.0
    M = 0.0


    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do i = imin, imax
        do j = jmin, jmax
            do k = kmin, kmax
                rho(i,j,k) = 0.0
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


end subroutine IC_Vaccum
