!****************************************************************
! Test_Accuracy.f90
!
! Created on: Jul 27, 2014
! Author: ninoy
!*****************************************************************

subroutine Test_Accuracy
#include <finclude/petscdef.h>
    use petscksp
    use MD_Parameter
    use MD_Quantity
    use MD_PETScQuantity
    use MD_Helper
    use MD_BoundaryCondition
    implicit none

    Vec pphi
    Vec err
    PetscErrorCode ierr
    character(len=100) :: format
    character(len=100) :: rankstr
    integer :: iteration, counter
    integer :: i, j, k, row, ierror
    real(kind=double) :: l1, l2, linf, r2
    PetscScalar, dimension(:), allocatable :: exacti
    PetscInt, dimension(:), allocatable :: rowi
    PetscScalar one

    if(rank .eq. master) then
        print *, ''
        print *, '*****************'
        print *, 'Testing Accuracy'
        print *, '*****************'
        print *, ''
    endif

    one = 1.0_rp

    call Init
    call Solver(iteration)
    call Output_3D(0)

    if(rank .eq. master) then
        print *, ''
        print *, '*********************'
        print *, 'Calculating Accuracy'
        print *, '*********************'
        print *, ''
    endif

    call VecDuplicate(phi, err, ierr)

    if(.not. allocated(exacti)) allocate(exacti(0:ilnc*jlnc*klnc-1))
    if(.not. allocated(rowi)) allocate(rowi(0:ilnc*jlnc*klnc-1))

    ! looping over inner cell
    !$OMP PARALLEL PRIVATE(i,j,k,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                row = GetLocalIndex(i,j,k)
                exacti(row) = GetBoundaryValue(x(i),y(j),z(k))
                rowi(row) = GetIndex(i,j,k)
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    call VecSetValues(err,ilnc*jlnc*klnc,rowi,exacti,INSERT_VALUES,ierr)
    call VecAssemblyBegin(err,ierr)
    call VecAssemblyEnd(err,ierr)

    if(allocated(exacti)) deallocate(exacti)
    if(allocated(rowi)) deallocate(rowi)

    call VecAXPBY(err,one,-one,phi,ierr)

    call VecNorm(err,NORM_1,l1,ierr)
    call VecNorm(err,NORM_2,l2,ierr)
    call VecNorm(err,NORM_INFINITY,linf,ierr)

    call MPI_Barrier(cart_comm, ierr)
    call MPI_Allreduce(l1, l1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cart_comm, ierr)
    call MPI_Allreduce(l2, l2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cart_comm, ierr)
    call MPI_Allreduce(linf, linf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, cart_comm, ierr)

    l1 = l1/(inc*jnc*knc)
    l2 = sqrt(l2**2/(inc*jnc*knc))

    if(rank .eq. master) then

        format = "(A12, I5, I5, I5)"
        write (*,format) 'resolution=',inc,jnc,knc
        format = "(A9, E10.3)"
        write (*,format) 'normL1=',l1
        write (*,format) 'normL2=',l2
        write (*,format) 'normInf=',linf

        open(unit=1, file='data/Accuracy.log',status='UNKNOWN',action= 'WRITE',&
            position='append',iostat= ierror)
        format = "(I3,A2,I3,A2,I3,A2,E10.3,A2,E10.3,A2,E10.3,A2,E10.3)"
        write(1,format) inc,',',jnc,',',knc,',',M,',',l1,',',l2,',',linf
        close(unit=1)

    end if

end subroutine Test_Accuracy
