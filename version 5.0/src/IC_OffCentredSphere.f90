 !****************************************************************  
 ! IC_OffCentredSphere.f90
 ! 
 ! Created on: Jul 27, 2014
 ! Author: ninoy
 !***************************************************************** 
subroutine IC_OffCentredSphere
    use MD_Parameter
    use MD_Quantity
    use MD_PETScQuantity
    implicit none

    integer :: i,j,k,counter
    real(kind=double) :: r2,unirho
    real(kind=double), dimension(0:2) :: xtmp
    integer :: ierr

    if(rank .eq. master) then
        print *, ''
        print *, '******************************************'
        print *, 'Initializing Scenerio: Off-Centred Sphere'
        print *, '******************************************'
        print *, ''
    endif

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
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin-1, kmax+1
        do j = jmin-1, jmax+1
            do i = imin-1, imax+1
                r2 = (x(i)-xc(0))**2+(y(j)-xc(1))**2+(z(k)-xc(2))**2
                if(r2 .le. R**2) then
                    rho(i,j,k) = unirho
                else
                    rho(i,j,k) = 0.0_rp
                end if
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

                xtmp(0) = xtmp(0) + x(i)*rho(i,j,k)
                xtmp(1) = xtmp(1) + y(j)*rho(i,j,k)
                xtmp(2) = xtmp(2) + z(k)*rho(i,j,k)

                M = M + rho(i,j,k)

            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    M = M*h**3

    call MPI_Barrier(cart_comm, ierr)
    call MPI_Allreduce(xtmp, xtmp, 3, MPI_DOUBLE_PRECISION, MPI_SUM, cart_comm, ierr)
    call MPI_Allreduce(M, M, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cart_comm, ierr)

    xc(0) = xtmp(0)*h**3/M
    xc(1) = xtmp(1)*h**3/M
    xc(2) = xtmp(2)*h**3/M

end subroutine IC_OffCentredSphere
