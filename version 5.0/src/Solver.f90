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
    PetscScalar, pointer :: phii(:)
    PetscErrorCode ierr
    KSPConvergedReason reason
    integer :: i,j,k,row,ierror
    integer, intent(out) :: iteration
    real(kind=double) :: normInf
    character(len=100) :: format
    character(len=100) :: rankstr

    call ComputeAvgDensity
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

!    call VecView(phi,PETSC_VIEWER_STDOUT_WORLD,ierr)

END SUBROUTINE Solver
