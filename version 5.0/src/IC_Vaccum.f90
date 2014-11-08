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

    if(rank .eq. master) then
        print *, ''
        print *, '******************************'
        print *, 'Initializing Scenerio: Vaccum'
        print *, '******************************'
        print *, ''
    endif

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M = 0.0_rp


    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin-1, kmax+1
        do j = jmin-1, jmax+1
            do i = imin-1, imax+1
                rho(i,j,k) = 0.0_rp
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


end subroutine IC_Vaccum
