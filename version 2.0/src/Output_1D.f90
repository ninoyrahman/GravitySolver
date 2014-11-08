 !****************************************************************  
 ! Output_1D.f90
 ! 
 ! Created on: Aug 13, 2014
 ! Author: ninoy
 !***************************************************************** 

subroutine Output_1D
#include <finclude/petscdef.h>
    use MD_Parameter
    use MD_Quantity
    use MD_Helper
    implicit none

    integer :: i,j,k,row
    integer :: ihalf, jhalf, khalf, ierror
!    real(kind=double), dimension(0:2) :: xcord
    PetscScalar phii(1)
    PetscOffset offset
    PetscErrorCode ierr
    character(len=100) :: format

    ihalf = (imin + imax)/2
    jhalf = (jmin + jmax)/2
    khalf = (kmin + kmax)/2

    open(unit=1000, file='data/rho_x.dat',status='UNKNOWN',action= 'WRITE', &
        position='rewind',iostat= ierror)

    open(unit=1001, file='data/rho_y.dat',status='UNKNOWN',action= 'WRITE', &
        position='rewind',iostat= ierror)

    open(unit=1002, file='data/rho_z.dat',status='UNKNOWN',action= 'WRITE', &
        position='rewind',iostat= ierror)

    open(unit=1003, file='data/phi_x.dat',status='UNKNOWN',action= 'WRITE', &
        position='rewind',iostat= ierror)

    open(unit=1004, file='data/phi_y.dat',status='UNKNOWN',action= 'WRITE', &
        position='rewind',iostat= ierror)

    open(unit=1005, file='data/phi_z.dat',status='UNKNOWN',action= 'WRITE', &
        position='rewind',iostat= ierror)

    open(unit=1006, file='data/error_x.dat',status='UNKNOWN',action= 'WRITE', &
        position='rewind',iostat= ierror)

    open(unit=1007, file='data/error_y.dat',status='UNKNOWN',action= 'WRITE',&
        position='rewind',iostat= ierror)

    open(unit=1008, file='data/error_z.dat',status='UNKNOWN',action= 'WRITE',&
        position='rewind',iostat= ierror)


    call VecGetArray(phi,phii,offset,ierr)

    format = "(E20.3E5, A4, E20.3E5)"

    j = jhalf
    k = khalf
    do i = imin, imax
!        call GetCordinate(i,j,k,xcord)
        row = GetIndex(i,j,k)
        write(1000,*) x(i), '    ', rho(i,j,k)
        write(1003,*) x(i), '    ', PetscRealPart(phii(offset + 1 + row))
        write(1006,*) x(i), '    ', error(i,j,k)
    end do

    i = jhalf
    k = khalf
    do j = jmin, jmax
!        call GetCordinate(i,j,k,xcord)
        row = GetIndex(i,j,k)
        write(1001,format) y(j), '    ', rho(i,j,k)
        write(1004,format) y(j), '    ', PetscRealPart(phii(offset + 1 + row))
        write(1007,format) y(j), '    ', error(i,j,k)
    end do

    i = ihalf
    j = jhalf
    do k = kmin, kmax
!        call GetCordinate(i,j,k,xcord)
        row = GetIndex(i,j,k)
        write(1002,format) z(k), '    ', rho(i,j,k)
        write(1005,format) z(k), '    ', PetscRealPart(phii(offset + 1 + row))
        write(1008,format) z(k), '    ', error(i,j,k)
    end do

    call VecRestoreArray(phi,phii,offset,ierr)

    do i = 1000, 1008
        close(unit=i)
    end do

end subroutine Output_1D
