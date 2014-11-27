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
    real(kind=double) :: unirho

    if(rank .eq. master) then
        print *, ''
        print *, '***************************************'
        print *, 'Initializing Scenerio: Uniform Density'
        print *, '***************************************'
        print *, ''
    endif

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M = 1.0_rp
    R = 1.5_rp*domainsize
    unirho = M/((4.0_rp/3.0_rp)*pi*R**3)

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin-1, kmax+1
        do j = jmin-1, jmax+1
            do i = imin-1, imax+1
                rho(i,j,k) = unirho
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL


end subroutine IC_UniformDensity
