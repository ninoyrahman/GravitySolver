 !****************************************************************  
 ! IC_InverseSquareDensity.f90
 ! 
 ! Created on: Aug 4, 2014
 ! Author: ninoy
 !***************************************************************** 

subroutine IC_InverseSquareDensity
    use MD_Parameter
    use MD_Quantity
    implicit none

    integer :: i,j,k,counter
    real(kind=double) :: r2
!    real(kind=double), dimension(0:2) :: xcord

    print *, ''
    print *, '*********************************************'
    print *, 'Initializing Scenerio: Inverse Square Density'
    print *, '*********************************************'
    print *, ''

    xc(0) = 0.0
    xc(1) = 0.0
    xc(2) = 0.0
    M = 0.0
    R = 1.0

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,r2)
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+:M) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
!                call GetCordinate(i,j,k,xcord)
                r2 = x(i)**2 + y(j)**2 + z(k)**2
                rho(i,j,k) = R**2/r2
                M = M + rho(i,j,k)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    M = M*h**3

end subroutine IC_InverseSquareDensity
