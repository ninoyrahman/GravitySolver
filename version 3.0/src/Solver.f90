!****************************************************************
! Solver.f90
!
! Created on: Jul 23, 2014
! Author: ninoy
!*****************************************************************

SUBROUTINE Solver(iteration)
#include <finclude/petscdef.h>
    use petscksp
    use MD_Definition
    use MD_Parameter
    use MD_Quantity
    use MD_PETScQuantity
    use MD_Helper
    use MD_BoundaryCondition
    implicit none

    Vec res
    PetscScalar, pointer :: resi(:)
    PetscErrorCode ierr
    KSPConvergedReason reason
    integer :: i,j,k,row
    integer, intent(out) :: iteration
    real(kind=double) :: normInf
    character(len=100) :: format
    character(len=100) :: rankstr

    call CalculateAvgDensity
    call ComputeRHS

    if(rank .eq. master) then
        print *, ''
        print *, '************************************************'
        print *, 'Solving Linear System of Equation for Potential'
        print *, '************************************************'
        print *, ''
    endif

    call KSPSolve(ksp,b,phi,ierr)

    format = "(A20, I3)"

    call KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD,ierr)
    if(rank .eq. master) then
        call KSPGetConvergedReason(ksp,reason,ierr)
        if(reason .lt. 0) then
            print *,'Diverged!!!, ', reason
        else
            call KSPGetIterationNumber(ksp,iteration,ierr)
            print format,'iteration number=', iteration
        end if
    end if

    call VecDuplicate(phi, res, ierr)
    call MatMult(A, phi, res, ierr)
    call VecAXPBY(res, 1.0_rp, -1.0_rp, b, ierr)
    call VecAssemblyBegin(res,ierr)
    call VecAssemblyEnd(res,ierr)

    normInf = 0.0

    call VecGetArrayF90(res,resi,ierr)

    !$OMP PARALLEL PRIVATE(i,j,k,row)
    !$OMP DO SCHEDULE(STATIC) REDUCTION(MAX:normInf) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                row = GetLocalIndex(i,j,k)
                if(normInf .lt. abs(PetscRealPart(resi(row+1)))) then
                    normInf = abs(PetscRealPart(resi(row+1)))
                end if
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    call VecRestoreArrayF90(res,resi,ierr)

    call VecDestroy(res, ierr)

    call MPI_Barrier(cart_comm, ierr)
    call MPI_Allreduce(normInf, normInf, 1, MPI_DOUBLE_PRECISION, MPI_MAX, cart_comm, ierr)

    if(rank .eq. master) then
        format = "(A20, E10.3)"
        print format, 'normInf(b-Ax)=', normInf
    end if

END SUBROUTINE Solver
