!****************************************************************
! Test_SolverConvergence.f90
!
! Created on: Jul 31, 2014
! Author: ninoy
!*****************************************************************
subroutine Test_SolverConvergence
#include <finclude/petscdef.h>
    use petscksp
    use MD_Parameter
    use MD_PETScQuantity
    use MD_Quantity
    use MD_Helper
    use MD_BoundaryCondition
    implicit none

    character(len=100) :: format
    integer :: iteration0,iterationp,ierror
    integer :: i,j,k,row
    real(kind=double) :: maxvalue, perturbation
    PetscScalar, pointer :: phii(:)
    PetscOffset offset
    PetscErrorCode ierr
    KSPConvergedReason reason

    if(rank .eq. master) then
        print *, ''
        print *, '***********************'
        print *, 'Solver Convergence Test'
        print *, '***********************'
        print *, ''
    end if

    perturbation = 1.0d-6

    call Init
    call Solver(iteration0)

    ! setting initial guess
    call VecNorm(phi,NORM_INFINITY,maxvalue,ierr)
    call VecGetArrayF90(phi,phii,ierr)
    !$OMP PARALLEL PRIVATE(i,j,k,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                row = localindex(i,j,k)
                phii(1 + row) = phii(1 + row) + dble(rand()-0.5_rp)*maxvalue*perturbation
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    call VecRestoreArrayF90(phi,phii,ierr)
    call VecAssemblyBegin(phi,ierr)
    call VecAssemblyEnd(phi,ierr)

    call KSPSolve(ksp,b,phi,ierr)
    call KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD,ierr)
    format = "(A20, I3)"
    if(rank .eq. master) then
        call KSPGetConvergedReason(ksp,reason,ierr)
        if(reason .lt. 0) then
            print *,'Diverged!!!, ', reason
        else
            call KSPGetIterationNumber(ksp,iterationp,ierr)
            print format,'iteration number=', iterationp
        end if
    end if

    if(rank .eq. master) then
        open(unit=1, file='data/SolverConvergence.log',status='UNKNOWN',action= 'WRITE',&
            position='append',iostat= ierror)
        format = "(I5,A2,I5,A2,I5,A2,E10.3,A2,I5,A2,I5,A2,I5)"
        write(1,format) inc,',',jnc,',',knc,',',perturbation,',',iteration0,',',iterationp,',',iteration0-iterationp
        close(unit=1)
    end if

end subroutine Test_SolverConvergence
