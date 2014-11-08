 !****************************************************************  
 ! IC_ExpDecayingDensity.f90
 ! 
 ! Created on: Aug 6, 2014
 ! Author: ninoy
 !***************************************************************** 
subroutine IC_ExpDecayingDensity
    use MD_Definition
    use MD_Parameter
    use MD_Quantity
    implicit none

    integer :: i,j,k
!    real(kind=double), dimension(0:2) :: xcord
    real(kind=double) :: r2

    print *, ''
    print *, '*****************************************************'
    print *, 'Initializing Scenerio: Exponentially Decaying Density'
    print *, '*****************************************************'
    print *, ''

    xc(0) = 0.0
    xc(1) = 0.0
    xc(2) = 0.0
    M = 0.0
    R = 100

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,r2)
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+:M) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
!                call GetCordinate(i,j,k,xcord)
                r2 = x(i)**2 + y(j)**2 + z(k)**2
                if(r2 .le. 0.1*0.1) then
                    rho(i,j,k) = exp(-R*(abs(x(i)) + abs(y(j)) + abs(z(k))))
                else
                    rho(i,j,k) = 0.0
                end if
                M = M + rho(i,j,k)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    M = M*h**3

end subroutine IC_ExpDecayingDensity
