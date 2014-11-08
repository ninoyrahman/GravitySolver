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
    use MD_BoundaryCondition
    implicit none

    character(len=100) :: format
    integer :: iteration, counter
    integer :: i,j,k,row,ierror
    real(kind=double) :: normL1, normL2, normInf, r2, exact
    PetscScalar phii(1)
    PetscOffset offset
    PetscErrorCode ierr

    coordinate = 0
    scenerio = 4

    call Init
    call Solver(iteration)
    call Output_3D(0)

    normL1 = 0.0_rp
    normL2 = 0.0_rp
    normInf = 0.0_rp
    counter = 0


    call VecGetArray(phi,phii,offset,ierr)

    !$OMP PARALLEL PRIVATE(i,j,k,row,r2,exact)
    !$OMP DO SCHEDULE(STATIC) REDUCTION(+:normL1,normL2,counter) REDUCTION(MAX:normInf) COLLAPSE(3)
    do k = kmin, kmax
        do j = jmin, jmax
            do i = imin, imax

                row = GetIndex(i,j,k)
                r2 = x(i)**2 + y(j)**2 + z(k)**2
                error(i,j,k) = 0.0_rp

!                if(rho(i,j,k) .lt. 1.0e-10) then
                    exact = GetBoundaryValue(x(i),y(j),z(k))
                    counter = counter + 1
                    error(i,j,k) = abs(PetscRealPart(phii(offset + 1 + row)) - exact)!/abs(exact)
!                else
!                    exact = GetR2Solution(x(i),y(j),z(k)) !GetSinSolution(x(i),y(j),z(k))
!                    counter = counter + 1
!                    error(i,j,k) = abs(PetscRealPart(phii(offset + 1 + row)) - exact)!/abs(exact)
!                end if

                normL1 = normL1 + error(i,j,k)
                normL2 = normL2 + error(i,j,k)*error(i,j,k)
                if(normInf .lt. error(i,j,k)) then
                    normInf = error(i,j,k)
                end if

            end do
        end do
    end do
    !$OMP END DO
    !$OMP END PARALLEL

    normL1 = normL1/counter
    normL2 = sqrt(normL2/counter)

    call VecRestoreArray(phi,phii,offset,ierr)

    open(unit=1, file='data/Accuracy.dat',status='UNKNOWN',action= 'WRITE',&
        position='append',iostat= ierror)

    format = "(I3,A2,I3,A2,I3,A2,E10.3,A2,E10.3,A2,E10.3,A2,E10.3)"
    write(1,format) inc,',',jnc,',',knc,',',M,',',normL1,',',normL2,',',normInf

    close(unit=1)

    call OutputError(0)
    call Output_1D

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
