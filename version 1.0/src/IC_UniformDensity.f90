 !****************************************************************  
 ! IC_UniformDensity.f90
 ! 
 ! Created on: Jul 29, 2014
 ! Author: ninoy
 !***************************************************************** 

subroutine IC_UniformDensity
    use MD_Definition
    use MD_Parameter
    use MD_Quantity
    implicit none

    integer :: i,j,k
    real(kind=double) :: c

    print *, ''
    print *, '***************************************'
    print *, 'Initializing Scenerio: Uniform Density'
    print *, '***************************************'
    print *, ''

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M = 1.0_rp
    R = 3.0_rp
    c = (M/((4.0_rp/3.0_rp)*pi*R**3))

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+:M) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                rho(i,j,k) = c
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


end subroutine IC_UniformDensity
