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
    use MD_PETScQuantity
    implicit none

    integer :: i,j,k
    real(kind=double) :: c, buffer
    integer :: ierr

    if(rank .eq. master) then
        print *, ''
        print *, '*************************************'
        print *, 'Initializing Scenerio: ploy6 Density'
        print *, '*************************************'
        print *, ''
    endif

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M = 0.0_rp
!    c = 30.0_rp/(4.0_rp*pi*G*5.0_rp*h)
    c = 30.0_rp/(4.0_rp*pi*G)

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin-1, kmax+1
        do j = jmin-1, jmax+1
            do i = imin-1, imax+1
                rho(i,j,k) = c*(x(i)**4 + y(j)**4 + z(k)**4)
!                rho(i,j,k) = c*((x(i) + 0.5_rp*h)**5 - (x(i) - 0.5_rp*h)**5 + &
!                                (y(j) + 0.5_rp*h)**5 - (y(j) - 0.5_rp*h)**5 + &
!                                (z(k) + 0.5_rp*h)**5 - (z(k) - 0.5_rp*h)**5)
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
                M = M + rho(i,j,k)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    buffer = M*h**3

    call MPI_Barrier(cart_comm, ierr)
    call MPI_Allreduce(buffer, M, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cart_comm, ierr)

end subroutine IC_Poly6Density
