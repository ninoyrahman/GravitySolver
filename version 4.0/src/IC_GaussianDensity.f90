!****************************************************************
! IC_GaussianDensity.f90
!
! Created on: Jul 31, 2014
! Author: ninoy
!*****************************************************************

subroutine IC_GaussianDensity
    use MD_Definition
    use MD_Parameter
    use MD_Quantity
    implicit none

    integer :: i,j,k
    real(kind=double) :: r2, c, alpha

    print *, ''
    print *, '****************************************'
    print *, 'Initializing Scenerio: Gaussian Density'
    print *, '****************************************'
    print *, ''

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M = 0.0_rp
    alpha = 64.0_rp
    c= ((pi/alpha)**(1.5_rp))/(8.0_rp*h**3)

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,r2)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin-1, kmax+1
        do j = jmin-1, jmax+1
            do i = imin-1, imax+1
                r2 = x(i)**2 + y(j)**2 + z(k)**2
                rho(i,j,k) = exp(-alpha*r2)
!                rho(i,j,k) = c*(erf(sqrt(alpha)*(x(i) + 0.5_rp*h)) - erf(sqrt(alpha)*(x(i) - 0.5_rp*h)))* &
!                               (erf(sqrt(alpha)*(y(j) + 0.5_rp*h)) - erf(sqrt(alpha)*(y(j) - 0.5_rp*h)))* &
!                               (erf(sqrt(alpha)*(z(k) + 0.5_rp*h)) - erf(sqrt(alpha)*(z(k) - 0.5_rp*h)))
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+:M) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                M = M + rho(i,j,k)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


    M = M*h**3

end subroutine IC_GaussianDensity
