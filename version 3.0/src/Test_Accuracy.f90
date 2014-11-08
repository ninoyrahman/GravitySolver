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

    character(len=100) :: format
    character(len=100) :: rankstr
    integer :: iteration, counter
    integer :: i,j,k,row,ierror,lrow
    real(kind=double) :: l1, l2, linf, r2, exact, err
    PetscScalar, pointer :: phii(:)
    PetscErrorCode ierr

    if(rank .eq. master) then
        print *, ''
        print *, '*****************'
        print *, 'Testing Accuracy'
        print *, '*****************'
        print *, ''
    endif

    coordinate = 0
    scenerio = 4

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

    l1 = 0.0_rp
    l2 = 0.0_rp
    linf = 0.0_rp

    call VecGetArrayF90(phi,phii,ierr)

    !$OMP PARALLEL PRIVATE(i,j,k,row,exact,err)
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+:l1,l2) REDUCTION(MAX:linf) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                row = GetLocalIndex(i,j,k)
                exact = GetBoundaryValue(x(i),y(j),z(k))
                err = abs(PetscRealPart(phii(row+1)) - exact)
                l1 = l1 + err
                l2 = l2 + err*err
                if(linf .lt. err) then
                    linf = err
                end if
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    call VecRestoreArrayF90(phi,phii,ierr)

    call MPI_Barrier(cart_comm, ierr)
    call MPI_Allreduce(l1, l1, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cart_comm, ierr)
    call MPI_Allreduce(l2, l2, 1, MPI_DOUBLE_PRECISION, MPI_SUM, cart_comm, ierr)
    call MPI_Allreduce(linf, linf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, cart_comm, ierr)

    if(rank .eq. master) then

        format = "(A12, I5, I5, I5)"
        write (*,format) 'resolution=',inc,jnc,knc
        format = "(A9, E10.3)"
        write (*,format) 'normL1=',l1/(inc*jnc*knc)
        write (*,format) 'normL2=',sqrt(l2/(inc*jnc*knc))
        write (*,format) 'normInf=',linf

    end if

end subroutine Test_Accuracy
