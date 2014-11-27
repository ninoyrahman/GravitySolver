 !****************************************************************  
 ! IC_CentredSphere.f90
 ! 
 ! Created on: Jul 26, 2014
 ! Author: ninoy
 !***************************************************************** 

subroutine IC_CentredSphere
    use MD_Parameter
    use MD_Quantity
    use MD_PETScQuantity
    implicit none

    integer :: i,j,k,counter
    real(kind=double) :: r2,unirho,buffer
    integer :: ierr

    if(rank .eq. master) then
        print *, ''
        print *, '**************************************'
        print *, 'Initializing Scenerio: Centred Sphere'
        print *, '**************************************'
        print *, ''
    endif

    xc(0) = 0.0_rp
    xc(1) = 0.0_rp
    xc(2) = 0.0_rp
    M = 0.0_rp
    R = 0.75_rp
    unirho = 1.0_rp/((4.0_rp/3.0_rp)*pi*R**3)

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,r2)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin-1, kmax+1
        do j = jmin-1, jmax+1
            do i = imin-1, imax+1
                r2 = x(i)**2 + y(j)**2 + z(k)**2
                if(r2 .le. R**2) then
                    rho(i,j,k) = unirho
                else
                    rho(i,j,k) = 1.0d-17
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
                M = M + rho(i,j,k)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    buffer = M*h**3

    call MPI_Barrier(cart_comm, ierr)
    call MPI_Allreduce(buffer, M, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cart_comm, ierr)

end subroutine IC_CentredSphere
