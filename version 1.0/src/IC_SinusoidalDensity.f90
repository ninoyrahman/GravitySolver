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
    use MD_CordinateTransform
    implicit none

    integer :: i,j,k
    real(kind=double), dimension(0:2) :: x

    print *, ''
    print *, '******************************************'
    print *, 'Initializing Scenerio: Sinusoidal Density'
    print *, '******************************************'
    print *, ''

    xc(0) = 0.0
    xc(1) = 0.0
    xc(2) = 0.0
    M = 0.0

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,x)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do i = imin, imax
        do j = jmin, jmax
            do k = kmin, kmax
                call GetCordinate(i,j,k,x)
                rho(i,j,k) = -(sin(x(0))+sin(x(1))+sin(x(2)))/(4*pi*G)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


end subroutine IC_SinusoidalDensity
