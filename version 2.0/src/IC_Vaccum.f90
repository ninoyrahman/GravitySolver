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

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M = 0.0_rp


    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                rho(i,j,k) = 0.0_rp
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


end subroutine IC_Vaccum
