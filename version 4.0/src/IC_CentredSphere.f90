 !****************************************************************  
 ! IC_CentredSphere.f90
 ! 
 ! Created on: Jul 26, 2014
 ! Author: ninoy
 !***************************************************************** 

subroutine IC_CentredSphere
    use MD_Parameter
    use MD_Quantity
    implicit none

    integer :: i,j,k,counter
    real(kind=double) :: r2,unirho
!    real(kind=double), dimension(0:2) :: xcord

    print *, ''
    print *, '**************************************'
    print *, 'Initializing Scenerio: Centred Sphere'
    print *, '**************************************'
    print *, ''

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M = 0.0_rp
    R = 0.5_rp
    unirho = 1.0_rp/((4.0_rp/3.0_rp)*pi*R**3)

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,r2)
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+:M) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
!                call GetCordinate(i,j,k,xcord)
                r2 = x(i)**2 + y(j)**2 + z(k)**2
                if(r2 .le. R**2) then
                    rho(i,j,k) = unirho
                else
                    rho(i,j,k) = 0.0_rp
                end if
                M = M + rho(i,j,k)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    M = M*h**3

end subroutine IC_CentredSphere
