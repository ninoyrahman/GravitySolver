!****************************************************************
! IC_OffCentredGaussian.f90
!
! Created on: Jul 31, 2014
! Author: ninoy
!*****************************************************************

subroutine IC_OffCentredGaussian
    use MD_Definition
    use MD_Parameter
    use MD_Quantity
    implicit none

    integer :: i,j,k
    real(kind=double) :: rl, rr, left, right, alpha, ixc, jxc, kxc


    print *, ''
    print *, '******************************************'
    print *, 'Initializing Scenerio: Off-Centred Density'
    print *, '******************************************'
    print *, ''

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M = 0.0_rp
    alpha = 64.0_rp
    ixc = 0.0_rp
    jxc = 0.0_rp
    kxc = 0.0_rp
    left = 1.0_rp
    right = 1.0_rp

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,rl,rr)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin-1, kmax+1
        do j = jmin-1, jmax+1
            do i = imin-1, imax+1
                rl = (x(i) + 0.5_rp)**2 + (y(j) - 0.0_rp)**2 + (z(k) - 0.0_rp)**2
                rr = (x(i) - 0.5_rp)**2 + (y(j) - 0.0_rp)**2 + (z(k) - 0.0_rp)**2
                rho(i,j,k) = left*exp(-alpha*rl) + right*exp(-alpha*rr)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+:M,ixc,jxc,kxc) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                M = M + rho(i,j,k)
                ixc = ixc + x(i)*rho(i,j,k)
                jxc = jxc + y(j)*rho(i,j,k)
                kxc = kxc + z(k)*rho(i,j,k)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    xc(0) = ixc/M
    xc(1) = jxc/M
    xc(2) = kxc/M

    M = M*h**3

    print *, xc(0),xc(1),xc(2)

end subroutine IC_OffCentredGaussian
