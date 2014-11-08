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

    integer :: i,j,k,row,ierror
    real(kind=double) :: normInf
    Vec res
    PetscScalar bi(1)
    PetscScalar resi(1)
    PetscOffset offset
    PetscOffset resoffset
    PetscScalar phia

    print *, ''
    print *, '************************************************'
    print *, 'Solving Linear System of Equation for Potential'
    print *, '************************************************'
    print *, ''

    call CalculateAvgDensity
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

    call VecCreate(PETSC_COMM_WORLD,res,ierr)
    call VecSetType(res,VECSEQ,ierr)
    call VecSetSizes(res,PETSC_DECIDE,inc*jnc*knc,ierr)
    call MatMult(A,phi,res,ierr)

    open(unit=5000, file='phi.dat',status='UNKNOWN', &
        action= 'WRITE', position='rewind',iostat= ierror)

    call VecGetArray(b,bi,offset,ierr)
    call VecGetArray(res,resi,resoffset,ierr)
    normInf = 0.0
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax
                row = GetIndex(i,j,k)
                if(normInf .lt. abs(PetscRealPart(bi(offset + 1 + row))-PetscRealPart(resi(resoffset + 1 + row)))) then
                    normInf = abs(PetscRealPart(bi(offset + 1 + row))-PetscRealPart(resi(resoffset + 1 + row)))
                end if
                call VecGetValues(phi, 1, row, phia, ierr)
                write(5000,*) row, '    ', PetscRealPart(phia)
            end do
        end do
    end do

    close(unit=5000)

    format = "(A20, E10.3)"
    print format, 'normInf(b-Ax)=', normInf

    call VecRestoreArray(b,bi,offset,ierr)
    call VecRestoreArray(res,resi,resoffset,ierr)
    call VecDestroy(res,ierr)

!    call VecView(phi,PETSC_VIEWER_STDOUT_WORLD,ierr)

END SUBROUTINE Solver
