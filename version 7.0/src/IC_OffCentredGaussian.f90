 !****************************************************************  
 ! IC_OffCentredSphere.f90
 ! 
 ! Created on: Jul 27, 2014
 ! Author: ninoy
 !***************************************************************** 
subroutine IC_OffCentredGaussian
    use MD_Parameter
    use MD_Quantity
    use MD_PETScQuantity
    implicit none

    integer :: i,j,k
    real(kind=double), dimension(0:2) :: xtmp
    real(kind=double) :: rl, rr, left, right, alpha, buffer
    integer :: ierr

    if(rank .eq. master) then
        print *, ''
        print *, '*******************************************'
        print *, 'Initializing Scenerio: Off-Centred Gaussian'
        print *, '*******************************************'
        print *, ''
    endif

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M = 0.0_rp
    alpha = 64.0_rp
    xtmp(0) = 0.0_rp
    xtmp(1) = 0.0_rp
    xtmp(2) = 0.0_rp
    left = 1.0_rp
    right = 1.0_rp


    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,rl,rr)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin-1, kmax+1
        do j = jmin-1, jmax+1
            do i = imin-1, imax+1
                rl = (x(i) + 0.5_rp)**2 + (y(j) - 0.0_rp)**2 + (z(k) - 0.0_rp)**2
                rr = (x(i) - 0.5_rp)**2 + (y(j) - 0.0_rp)**2 + (z(k) - 0.0_rp)**2
                rho(i,j,k) = left*exp(-alpha*rl) + right*exp(-alpha*rr)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+:M,xtmp) COLLAPSE(3)
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

    buffer = M*h**3

    call MPI_Barrier(cart_comm, ierr)
    call MPI_Allreduce(xtmp, xc, 3, MPI_DOUBLE_PRECISION, MPI_SUM, cart_comm, ierr)
    call MPI_Allreduce(buffer, M, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cart_comm, ierr)

    xc(0) = xc(0)*h**3/M
    xc(1) = xc(1)*h**3/M
    xc(2) = xc(2)*h**3/M

    if(rank .eq. master) print *, xc(0), xc(1), xc(2)

end subroutine IC_OffCentredGaussian
