 !****************************************************************  
 ! IC_SinusoidalDensity.f90
 ! 
 ! Created on: Jul 30, 2014
 ! Author: ninoy
 !***************************************************************** 
subroutine IC_SinusoidalDensity
    use MD_Definition
    use MD_Parameter
    use MD_Quantity
    implicit none

    integer :: i,j,k
!    real(kind=double), dimension(0:2) :: xcord
    real(kind=double) :: r2, c

    print *, ''
    print *, '******************************************'
    print *, 'Initializing Scenerio: Sinusoidal Density'
    print *, '******************************************'
    print *, ''

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M = 0.0_rp
    R = 0.5_rp
!    c = (pi/(48.0_rp*R**3))*(2.0_rp*R/(pi*h))
    c = (pi/(48.0_rp*R**3))

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,r2)
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+:M) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
!                call GetCordinate(i,j,k,xcord)
                r2 = x(i)**2 + y(j)**2 + z(k)**2
!                if(r2 .le. R**2) then
                    rho(i,j,k) = c*(cos(pi*x(i)/(2.0_rp*R)) + cos(pi*y(j)/(2.0_rp*R)) + cos(pi*z(k)/(2.0_rp*R)))
!                    rho(i,j,k) = c*(sin(pi*(x(i)+0.5_rp*h)/(2.0_rp*R)) - sin(pi*(x(i)-0.5_rp*h)/(2.0_rp*R)) + &
!                                    sin(pi*(y(j)+0.5_rp*h)/(2.0_rp*R)) - sin(pi*(y(j)-0.5_rp*h)/(2.0_rp*R)) + &
!                                    sin(pi*(z(k)+0.5_rp*h)/(2.0_rp*R)) - sin(pi*(z(k)-0.5_rp*h)/(2.0_rp*R)))
!                else
!                    rho(i,j,k) = 0.0_rp
!                end if
                M = M + rho(i,j,k)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    M = M*h**3

end subroutine IC_SinusoidalDensity
