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
    use MD_Helper
    use MD_CordinateTransform
    use MD_BoundaryCondition
    implicit none

    character(len=100) :: format
    integer :: iteration, counter
    integer :: i,j,k,row,ierror
    real(kind=double) :: normL1, normL2, normInf
    real(kind=double), dimension(0:2) :: x
    PetscScalar phii(1)
    PetscOffset offset
    PetscErrorCode ierr

    coordinate = 0
    scenerio = 4

    call Init
    call Solver(iteration)
    call Output(0)

    normL1 = 0.0
    normL2 = 0.0
    normInf = 0.0

    call VecGetArray(phi,phii,offset,ierr)

    counter = 0
    i = imin
    !$OMP PARALLEL PRIVATE(j,k,row,x)
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+:normL1,normL2) COLLAPSE(2)
    do j = jmin, jmax
        do k = kmin, kmax
!            if(rho(i,j,k) .eq. 0.0) then
                row = GetIndex(i,j,k)
                call GetCordinate(i,j,k,x)
!                print *, PetscRealPart(phii(offset + 1 + row)), GetBoundaryValue(x)
                normL1 = normL1 + abs(PetscRealPart(phii(offset + 1 + row)) - GetBoundaryValue(x))
                normL2 = normL2 + (PetscRealPart(phii(offset + 1 + row)) - GetBoundaryValue(x))* &
                    (PetscRealPart(phii(offset + 1 + row)) - GetBoundaryValue(x))
                counter = counter + 1
                if(normInf .lt. abs(PetscRealPart(phii(offset + 1 + row)) - GetBoundaryValue(x))) then
                    normInf = abs(PetscRealPart(phii(offset + 1 + row)) - GetBoundaryValue(x))
                end if
!            end if
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    call VecRestoreArray(phi,phii,offset,ierr)

    normL1 = normL1/counter
    normL2 = sqrt(normL2)/counter

    open(unit=1, file='data/Accuracy.dat',status='UNKNOWN',action= 'WRITE',&
        position='append',iostat= ierror)

    format = "(I3,A2,I3,A2,I3,A2,E10.3,A2,E10.3,A2,E10.3)"
    write(1,format) inc,',',jnc,',',knc,',',normL1,',',normL2,',',normInf

    close(unit=1)

    print *, ''
    print *, '**************'
    print *, 'Accuracy Test'
    print *, '**************'
    print *, ''

    format = "(A12, I5, I5, I5)"
    write (*,format) 'resolution=',inc,jnc,knc
    format = "(A8, E10.3)"
    write (*,format) 'normL1=',normL1
    write (*,format) 'normL2=',normL2
    write (*,format) 'normInf=',normInf

end subroutine Test_Accuracy
