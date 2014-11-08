 !****************************************************************  
 ! IC_OffCentredSphere.f90
 ! 
 ! Created on: Jul 27, 2014
 ! Author: ninoy
 !***************************************************************** 
subroutine IC_OffCentredSphere
    use MD_Parameter
    use MD_Quantity
    implicit none

    integer :: i,j,k,counter
    real(kind=double) :: r2,unirho
!    real(kind=double), dimension(0:2) :: xcord
    real(kind=double), dimension(0:2) :: xtmp

    print *, ''
    print *, '******************************************'
    print *, 'Initializing Scenerio: Off-Centred Sphere'
    print *, '******************************************'
    print *, ''

    xc(0) = 0.1_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M = 0.0_rp
    R = 0.2_rp
    unirho = 1.0_rp/((4.0_rp/3.0_rp)*pi*R**3)
    xtmp(0) = 0.0_rp
    xtmp(1) = 0.0_rp
    xtmp(2) = 0.0_rp

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,r2)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3) REDUCTION(+:M)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
!                call GetCordinate(i,j,k,xcord)
                r2 = (x(i)-xc(0))**2+(y(j)-xc(1))**2+(z(k)-xc(2))**2
                if(r2 .le. R**2) then
                    rho(i,j,k) = unirho
                    xtmp(0) = xtmp(0) + x(i)*rho(i,j,k)
                    xtmp(1) = xtmp(1) + y(j)*rho(i,j,k)
                    xtmp(2) = xtmp(2) + z(k)*rho(i,j,k)
                else
                    rho(i,j,k) = 0.0_rp
                end if
                M = M + rho(i,j,k)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    xc(0) = xtmp(0)/M
    xc(1) = xtmp(1)/M
    xc(2) = xtmp(2)/M

    M = M*h**3

end subroutine IC_OffCentredSphere
