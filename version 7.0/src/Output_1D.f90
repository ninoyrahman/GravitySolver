 !****************************************************************  
 ! Output_1D.f90
 ! 
 ! Created on: Aug 13, 2014
 ! Author: ninoy
 !***************************************************************** 

subroutine Output_1D(ig,jg,kg)
#include <finclude/petscdef.h>
    use petscksp
    use MD_Parameter
    use MD_Quantity
    use MD_PETScQuantity
    use MD_Helper
    implicit none

    integer, intent(in) :: ig,jg,kg
    integer :: i,j,k,row,ierror
    character(len=100) :: format
    character(len=100) :: rankstr
    PetscScalar, pointer :: phii(:)
    PetscErrorCode ierr

    if(rank .eq. master) then
        format = "(A10, I3)"
        print *, ''
        print *, '*******************'
        print *, 'writing 1D output '
        print *, '*******************'
        print *, ''
    endif

    format = '(I5.5)'
    write (rankstr, format) rank

    call VecGetArrayF90(phi,phii,ierr)

    format = "(E20.3E5, A4, E20.3E5)"

    if((jg .ge. jgmin .and. jg .le. jgmax) .and. (kg .ge. kgmin .and. kg .le. kgmax)) then
        open(unit=1000, file='data/rho_x.'//trim(rankstr)//'.dat',status='UNKNOWN',action= 'WRITE', &
            position='rewind',iostat= ierror)
!        open(unit=1003, file='data/phi_x.'//trim(rankstr)//'.dat',status='UNKNOWN',action= 'WRITE', &
!            position='rewind',iostat= ierror)
        j = jg-jgmin+ngh
        k = kg-kgmin+ngh
        do i = imin, imax
!            row = localindex(i,j,k)
            write(1000,format) x(i), '    ', rho(i,j,k)
!            write(1003,format) x(i), '    ', PetscRealPart(phii(1 + row))
        end do
        close(unit=1000)
!        close(unit=1003)
    end if

    if((ig .ge. igmin .and. ig .le. igmax) .and. (kg .ge. kgmin .and. kg .le. kgmax)) then
        i = ig-igmin+ngh
        k = kg-kgmin+ngh
        open(unit=1001, file='data/rho_y.'//trim(rankstr)//'.dat',status='UNKNOWN',action= 'WRITE', &
            position='rewind',iostat= ierror)
!        open(unit=1004, file='data/phi_y.'//trim(rankstr)//'.dat',status='UNKNOWN',action= 'WRITE', &
!            position='rewind',iostat= ierror)
        do j = jmin, jmax
!            row = localindex(i,j,k)
            write(1001,format) y(j), '    ', rho(i,j,k)
!            write(1004,format) y(j), '    ', PetscRealPart(phii(1 + row))
        end do
        close(unit=1001)
!        close(unit=1004)
    end if

    if((ig .ge. igmin .and. ig .le. igmax) .and. (jg .ge. jgmin .and. jg .le. jgmax)) then
        i = ig-igmin+ngh
        j = jg-jgmin+ngh
        open(unit=1002, file='data/rho_z.'//trim(rankstr)//'.dat',status='UNKNOWN',action= 'WRITE', &
            position='rewind',iostat= ierror)
!        open(unit=1005, file='data/phi_z.'//trim(rankstr)//'.dat',status='UNKNOWN',action= 'WRITE', &
!            position='rewind',iostat= ierror)
        do k = kmin, kmax
!            row = localindex(i,j,k)
            write(1002,format) z(k), '    ', rho(i,j,k)
!            write(1005,format) z(k), '    ', PetscRealPart(phii(1 + row))
        end do
        close(unit=1002)
!        close(unit=1005)
    end if

    call VecRestoreArrayF90(phi,phii,ierr)

end subroutine Output_1D
