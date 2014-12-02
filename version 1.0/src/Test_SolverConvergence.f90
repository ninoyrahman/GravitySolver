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
    use MD_Quantity
    use MD_Helper
    use MD_BoundaryCondition
    implicit none

    character(len=100) :: format
    integer :: iteration0,iterationp,ierror
    integer :: i,j,k,row
    real(kind=double) :: maxvalue, perturbation
    PetscScalar phii(1)
    PetscOffset offset
    PetscErrorCode ierr
    KSPConvergedReason reason

    print *, ''
    print *, '***********************'
    print *, 'Solver Convergence Test'
    print *, '***********************'
    print *, ''


    coordinate = 0
    scenerio = 6
    perturbation = 1.0d-1

    call Init
    call Solver(iteration0)

    ! setting initial guess
    call VecNorm(phi,NORM_INFINITY,maxvalue,ierr)
    call VecGetArray(phi,phii,offset,ierr)
    !$OMP PARALLEL PRIVATE(i,j,k,row)
    !$OMP DO SCHEDULE(STATIC) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                row = GetIndex(i,j,k)
                phii(offset + 1 + row) = phii(offset + 1 + row) + dble(rand())*maxvalue*perturbation
            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL
    call VecRestoreArray(phi,phii,offset,ierr)
    call VecAssemblyBegin(phi,ierr)
    call VecAssemblyEnd(phi,ierr)

    call KSPSolve(ksp,b,phi,ierr)
    format = "(A20, I3)"
    call KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD,ierr)
    call KSPGetConvergedReason(ksp,reason,ierr)
    if(reason .lt. 0) then
        print *,'Diverged!!!, ', reason
    else
        call KSPGetIterationNumber(ksp,iterationp,ierr)
        print format,'iteration number=', iterationp
    end if

    open(unit=1, file='data/SolverConvergence.dat',status='UNKNOWN',action= 'WRITE',&
        position='append',iostat= ierror)
    format = "(I5,A2,I5,A2,I5,A2,E10.3,A2,I5,A2,I5,A2,I5)"
    write(1,format) inc,',',jnc,',',knc,',',perturbation,',',iteration0,',',iterationp,',',iteration0-iterationp
    close(unit=1)

end subroutine Test_SolverConvergence
