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
    use MD_Helper
    implicit none

    character(len=100) :: format
    integer, intent(out) :: iteration
    PetscViewer viewer
    PetscErrorCode ierr
    KSPConvergedReason reason

    integer :: i,j,k,row
    real(kind=double) :: normInf
    Vec y
    PetscScalar bi(1)
    PetscScalar yi(1)
    PetscOffset offset
    PetscOffset yoffset


    print *, ''
    print *, '************************************************'
    print *, 'Solving Linear System of Equation for Potential'
    print *, '************************************************'
    print *, ''

    call ComputeRHS

    call KSPSolve(ksp,b,phi,ierr)

    format = "(A20, I3)"

    call KSPView(ksp, PETSC_VIEWER_STDOUT_WORLD,ierr)
    call KSPGetConvergedReason(ksp,reason,ierr)
    if(reason .lt. 0) then
        print *,'Diverged!!!, ', reason
    else
        call KSPGetIterationNumber(ksp,iteration,ierr)
        print format,'iteration number=', iteration
    end if

    call VecCreate(PETSC_COMM_WORLD,y,ierr)
    call VecSetType(y,VECSEQ,ierr)
    call VecSetSizes(y,PETSC_DECIDE,inc*jnc*knc,ierr)
    call MatMult(A,phi,y,ierr)


    call VecGetArray(b,bi,offset,ierr)
    call VecGetArray(y,yi,yoffset,ierr)
    normInf = 0.0
    do i = imin, imax
        do j = jmin, jmax
            do k = kmin, kmax
                row = GetIndex(i,j,k)
                if(normInf .lt. abs(PetscRealPart(bi(offset + 1 + row))-PetscRealPart(yi(yoffset + 1 + row)))) then
                    normInf = abs(PetscRealPart(bi(offset + 1 + row))-PetscRealPart(yi(yoffset + 1 + row)))
                end if

            end do
        end do
    end do

    format = "(A20, E10.3)"
    print format, 'normInf(b-Ax)=', normInf

    call VecRestoreArray(b,bi,offset,ierr)
    call VecRestoreArray(y,yi,yoffset,ierr)
    call VecDestroy(y,ierr)



END SUBROUTINE Solver
