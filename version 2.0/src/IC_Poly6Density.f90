 !****************************************************************  
 ! IC_Poly6Density.f90
 ! 
 ! Created on: Jul 30, 2014
 ! Author: ninoy
 !***************************************************************** 
subroutine IC_Poly6Density
    use MD_Definition
    use MD_Parameter
    use MD_Quantity
    implicit none

    integer :: i,j,k
!    real(kind=double), dimension(0:2) :: xcord
    real(kind=double) :: c

    print *, ''
    print *, '*************************************'
    print *, 'Initializing Scenerio: ploy6 Density'
    print *, '*************************************'
    print *, ''

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M = 0.0_rp
!    c = 30.0_rp/(4.0_rp*pi*G*5.0_rp*h)
    c = 30.0_rp/(4.0_rp*pi*G)

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
!                call GetCordinate(i,j,k,xcord)
                rho(i,j,k) = c*(x(i)**4 + y(j)**4 + z(k)**4)
!                rho(i,j,k) = c*((x(i) + 0.5_rp*h)**5 - (x(i) - 0.5_rp*h)**5 + &
!                                (y(j) + 0.5_rp*h)**5 - (y(j) - 0.5_rp*h)**5 + &
!                                (z(k) + 0.5_rp*h)**5 - (z(k) - 0.5_rp*h)**5)
                M = M + rho(i,j,k)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    M = M*h**3


end subroutine IC_Poly6Density
